#!/usr/bin/env python3
"""
build_plasmid_db.py

1. Assembles each sample with Unicycler (hybrid: long + short reads)
2. Extracts plasmid contigs (circular + size filter)
3. Builds combined plasmid FASTA database

Read path conventions:
  Long:        {LONG_READS_DIR}/{run_accession}.fastq
  Short paired:{SHORT_READS_DIR}/{illumina_run_accession}_1.fastq + _2.fastq
  Short single:{SHORT_READS_DIR}/{illumina_run_accession}.fastq
"""

import gzip
import re
import subprocess
import argparse
from pathlib import Path
import pandas as pd

LONG_READS_DIR  = Path("/wynton/home/rotation/valeria-se/Mueller_Lab/long_reads")
SHORT_READS_DIR = Path("/wynton/home/rotation/valeria-se/Mueller_Lab/short_reads")
DEFAULT_CHROM_CUTOFF = 2_000_000


def open_fasta(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header, chunks = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def resolve_short_reads(short_acc: str) -> tuple[str, list[Path]] | tuple[None, None]:
    """
    Returns (mode, paths) where mode is 'paired' or 'unpaired'.
    Checks for paired (_1/_2) first, then falls back to single-end.
    Returns (None, None) if no files found.
    """
    r1 = SHORT_READS_DIR / f"{short_acc}_1.fastq"
    r2 = SHORT_READS_DIR / f"{short_acc}_2.fastq"
    unpaired = SHORT_READS_DIR / f"{short_acc}.fastq"

    if r1.exists() and r2.exists():
        return "paired", [r1, r2]
    elif r1.exists() and not r2.exists():
        print(f"  [WARN] Found {r1.name} but not {r2.name} — treating as unpaired")
        return "unpaired", [r1]
    elif unpaired.exists():
        return "unpaired", [unpaired]
    else:
        return None, None


def run_unicycler(long_reads: Path, short_mode: str, short_paths: list[Path],
                  sample: str, out_root: Path, threads: int, mode: str) -> Path | None:
    out_dir  = out_root / sample
    assembly = out_dir / "assembly.fasta"
    out_dir.mkdir(parents=True, exist_ok=True)

    if assembly.exists():
        print(f"  [SKIP] assembly already exists: {assembly}")
        return assembly

    cmd = ["unicycler"]

    if short_mode == "paired":
        cmd += ["-1", str(short_paths[0]), "-2", str(short_paths[1])]
        print(f"  Short reads: paired ({short_paths[0].name}, {short_paths[1].name})")
    else:
        cmd += ["-s", str(short_paths[0])]
        print(f"  Short reads: unpaired ({short_paths[0].name})")

    cmd += [
        "-l", str(long_reads),
        "-o", str(out_dir),
        "--threads", str(threads),
        "--mode", mode,
    ]

    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not assembly.exists():
        print(f"  [ERROR] Unicycler failed for {sample}:\n{result.stderr[-800:]}")
        return None

    return assembly


def is_circular(header: str) -> bool:
    return bool(re.search(r"circular=true", header.lower()))


def classify_contig(header: str, seq: str, size_cutoff: int) -> tuple[bool, str]:
    length   = len(seq)
    circular = is_circular(header)
    print(f"    CONTIG | len={length:>10,} | circular={str(circular):<5} | header='{header[:80]}'")
    if length >= size_cutoff:
        return False, f"chromosomal size ({length:,} bp)"
    if not circular:
        return False, "not circular"
    return True, f"plasmid — size={length:,}, circular"


def build_plasmid_db(sample_df, out_root, threads, size_cutoff, output, mode):
    total = kept = failed = missing = 0
    out_root.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as out_fh:
        for _, row in sample_df.iterrows():
            sample    = row["sample_accession"]
            long_acc  = row["run_accession"]
            short_acc = row["illumina_run_accession"]

            print(f"\n── {sample} (long={long_acc}, short={short_acc}) ──")

            # Resolve long reads
            long_reads = LONG_READS_DIR / f"{long_acc}.fastq"
            if not long_reads.exists():
                print(f"  [WARN] Missing long reads: {long_reads}")
                missing += 1
                continue

            # Resolve short reads — paired or unpaired
            short_mode, short_paths = resolve_short_reads(short_acc)
            if short_mode is None:
                print(f"  [WARN] No short reads found for {short_acc} in {SHORT_READS_DIR}")
                missing += 1
                continue

            assembly = run_unicycler(long_reads, short_mode, short_paths,
                                     sample, out_root, threads, mode)
            if assembly is None:
                failed += 1
                continue

            sample_kept = 0
            for header, seq in open_fasta(assembly):
                total += 1
                keep, reason = classify_contig(header, seq, size_cutoff)
                if keep:
                    out_fh.write(f">{sample}__{header}\n{seq}\n")
                    sample_kept += 1
                    kept += 1

            print(f"  → Kept {sample_kept} plasmid contig(s)")

    print(f"\n── Summary ──────────────────────────────────────────")
    print(f"  Missing input files    : {missing}")
    print(f"  Unicycler failures     : {failed}")
    print(f"  Total contigs seen     : {total}")
    print(f"  Plasmid contigs kept   : {kept}")
    print(f"  Output                 : {output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads",     type=int, default=8)
    parser.add_argument("--size-cutoff", type=int, default=DEFAULT_CHROM_CUTOFF)
    parser.add_argument("--outdir",      default="unicycler_assemblies")
    parser.add_argument("--output",      default="plasmid_db.fa")
    parser.add_argument("--mode",        default="normal",
                        choices=["conservative", "normal", "bold"])
    parser.add_argument("--csv",         default="unicycler.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.csv).dropna(subset=["run_accession", "illumina_run_accession"])
    print(f"Loaded {len(df)} samples from {args.csv}")

    build_plasmid_db(
        sample_df   = df,
        out_root    = Path(args.outdir),
        threads     = args.threads,
        size_cutoff = args.size_cutoff,
        output      = Path(args.output),
        mode        = args.mode,
    )


if __name__ == "__main__":
    main()