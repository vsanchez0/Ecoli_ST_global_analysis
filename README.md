# Ecoli_ST_global_analysis

### Overview
This repository corresponds to our project investigating how antimicrobial resistance (AMR) gene content shapes global dispersion patterns of *Escherichia coli* sequence types (STs).

### Research Aims
We seek to determine whether *E. coli* STs exhibit consistent, predictable relationships among:
* Genomic traits (AMR gene content, virulence factors)
* Global dispersion patterns (geographic range, introduction events)
* AMR gene gain/loss dynamics over time and across regions

### Methods
**1. Genomic Data Collection**
* Collection of publicly available *E. coli* genome sequences with associated metadata
- [x] Select four STs with around 10,000 samples per ST
- [x] Figure out which sequences have corresponding granular time and location data (months, states)
- [x] Figure out how many of these samples have long read data
- [x] Download assemblies

**2. AMR Gene Profiling**
* AMR gene detection using CARD/RGI
- [x] Annotate assemblies

**3. Phylogenetic Reconstruction**
* [x] Generate core-genome alignment (Gubbins)
* [ ] Temporal signal check and outlier detection (TreeTime)
* [ ] Split ST131 into sub-clades sharing MRCA within ~50 years; rerun alignment per clade
* [ ] Phylogenetic inference per ST with subsampling (BEAST2)
* [ ] Ancestral state reconstruction of AMR gene gain/loss events (DTA)

**4. BEAUti / BEAST XML Setup**
* [ ] Prepare per-ST alignment FASTA with tip dates embedded in sequence names (`isolate|YYYY-MM-DD`)
* [ ] Prepare binary AMR trait file (isolate × gene presence/absence matrix; exclude genes at 0% or 100% prevalence)
* [ ] Prepare discrete location trait file (HHS region per isolate)
* [ ] Import alignment + traits into BEAUti using template XML; confirm DTA model setup
* [ ] Run BEAST DTA for AMR genes and HHS migration

**5. Visualization & Interpretation**
* [ ] Flag long-read (Oxford Nanopore / PacBio) isolates in metadata for tree annotation
* [ ] Visualize flagged long-read isolates on trees
* [ ] Plot AMR gene gain/loss events mapped onto time-scaled tree
