import baltic as bt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import re
import tempfile
import os
import math
from datetime import datetime


GENES = [
    {
        "key":   "sul1",
        "label": "sul1",
        "file":  "st131_clades/sul1/sul1_clade_1.nwk",
    },
    {
        "key":   "sul2",
        "label": "sul2",
        "file":  "st131_clades/sul2/sul2_clade_1.nwk",
    },
    {
        "key":   "Escherichia_coli_mdfA",
        "label": "mdfA",
        "file":  "st131_clades/mdfA/Escherichia_coli_mdfA_clade_1.nwk",
    },
]

OUTPUT = "figures/test.png"

GENE_PALETTE = [
    "#D4A020",   # Kyburg Gold
    "#18A3AC",   # UCSF Teal
    "#6B2D90",   # Hutch Purple
]

PRESENT_COLOR = "#C0362C"
ABSENT_COLOR  = "#FFFFFF"
BORDER_COLOR  = "#CCCCCC"
MISSING_COLOR = "#888886"

GAIN_MARKER = "^"
LOSS_MARKER = "v"
MARKER_SIZE = 60


assert 1 <= len(GENES) <= 3, "GENES must contain 1–3 entries."

for i, gene in enumerate(GENES):
    gene["color"] = GENE_PALETTE[i]


def parse_translate_block(text):
    """Return {number_str: taxon_name} from a NEXUS Translate block."""
    translate_map = {}
    tb = re.search(r'Translate\s+(.*?);', text, re.DOTALL | re.IGNORECASE)
    if tb:
        for m in re.finditer(r'(\d+)\s+(\S+)', tb.group(1)):
            translate_map[m.group(1)] = m.group(2).rstrip(',')
    return translate_map


def extract_last_tree_newick(filepath):
    """
    Return (translate_map, newick_string) for the last sampled tree.
 
    Handles three formats:
      - BEAST .trees nexus  (Translate block + 'tree STATE_...' lines)
      - Generic NEXUS       ('tree ...' lines, no Translate block)
      - Plain Newick        (one or more bare Newick strings; last line used)
    """
    with open(filepath) as f:
        text = f.read()
 
    nexus_trees = re.findall(
        r'^tree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)?(.+)$',
        text, re.MULTILINE | re.IGNORECASE,
    )
    if nexus_trees:
        translate_map = parse_translate_block(text)
        return translate_map, nexus_trees[-1].strip()
 
    newick_lines = [
        ln.strip() for ln in text.splitlines()
        if ln.strip().startswith('(')
    ]
    if newick_lines:
        return {}, newick_lines[-1]
 
    raise ValueError(
        f"No recognisable tree found in {filepath}.\n"
        "Expected a BEAST/NEXUS file with 'tree ... = ...' lines, "
        "or a plain Newick file whose tree lines begin with '('."
    )


def load_baltic_tree(newick_string, translate_map):
    """Write newick to a temp file, load with baltic, remap tip names."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write(newick_string + '\n')
        tmppath = f.name
    try:
        tree = bt.loadNewick(tmppath, absoluteTime=False)
    finally:
        os.unlink(tmppath)
    for node in tree.Objects:
        if node.branchType == 'leaf' and node.name in translate_map:
            node.name = translate_map[node.name]
    return tree


def build_state_map(tree, gene_key):
    """
    Read gene state from node.traits (populated by baltic during newick parsing).
    Returns {id(node): int}.
    """
    state_map = {}
    for node in tree.Objects:
        traits = getattr(node, 'traits', {}) or {}
        if gene_key in traits:
            try:
                state_map[id(node)] = int(traits[gene_key])
            except (ValueError, TypeError):
                pass
    return state_map


def rekey_to_reference_tree(ref_tree, src_tree, src_map):
    """
    Translate a state map keyed by id(src_tree node) to id(ref_tree node),
    using positional correspondence (both trees must share the same topology).
    """
    out = {}
    for n_ref, n_src in zip(ref_tree.Objects, src_tree.Objects):
        if id(n_src) in src_map:
            out[id(n_ref)] = src_map[id(n_src)]
    return out


def get_transitions(node_state_map, tree_objects):
    """
    Detect parent→child state changes on every branch.
    parent=0, child=1 → gain;  parent=1, child=0 → loss.
    Marker is placed at the branch midpoint.
    """
    transitions = []
    for node in tree_objects:
        if node.parent is None:
            continue
        p = node_state_map.get(id(node.parent))
        c = node_state_map.get(id(node))
        if p is None or c is None or p == c:
            continue
        transitions.append({
            "x":    (node.parent.x + node.x) / 2.0,
            "y":    node.y,
            "type": "gain" if c == 1 else "loss",
        })
    return transitions


def scatter_markers_with_offset(ax, all_markers, y_step=0.28):
    """
    Plot transition markers, offsetting vertically when multiple markers
    fall on the same branch (same x, y position).

    all_markers: list of (x, y, marker_symbol, color)
    """
    from collections import defaultdict

    groups = defaultdict(list)
    for item in all_markers:
        key = (round(item[0], 4), round(item[1], 4))
        groups[key].append(item)

    for items in groups.values():
        n = len(items)
        offsets = [(i - (n - 1) / 2.0) * y_step for i in range(n)]
        for (x, y, marker, color), dy in zip(items, offsets):
            ax.scatter(x, y + dy,
                       marker=marker, color=color,
                       s=MARKER_SIZE, zorder=200,
                       linewidths=0.5, edgecolors="white")


def date_to_decimal(date_str):
    """'2020-09-13' → 2020.70"""
    d = datetime.strptime(date_str.strip(), "%Y-%m-%d")
    start = datetime(d.year, 1, 1)
    end   = datetime(d.year + 1, 1, 1)
    return d.year + (d - start).days / (end - start).days

print(f"Loading {len(GENES)} gene(s) …")

# Load the reference tree from the first gene's file
translate_map, ref_newick = extract_last_tree_newick(GENES[0]["file"])
ref_tree = load_baltic_tree(ref_newick, translate_map)

state_maps      = {}
gene_transitions = {}

for gene in GENES:
    label = gene["label"]
    fpath = gene["file"]
    key   = gene["key"]

    _, newick = extract_last_tree_newick(fpath)
    src_tree  = load_baltic_tree(newick, translate_map)
    raw_map   = build_state_map(src_tree, key)
    state_maps[label] = rekey_to_reference_tree(ref_tree, src_tree, raw_map)

    gene_transitions[label] = get_transitions(state_maps[label], ref_tree.Objects)


tips      = [n for n in ref_tree.Objects if n.branchType == 'leaf']
tip_names = [t.name for t in tips]

tip_decimal_dates = {}
for name in tip_names:
    if '|' in name:
        try:
            tip_decimal_dates[name] = date_to_decimal(name.split('|', 1)[1])
        except ValueError:
            pass

x_max  = max(t.x for t in tips)
x_min  = min(n.x for n in ref_tree.Objects if n.x is not None)
x_span = x_max - x_min

if tip_decimal_dates:
    max_date = max(tip_decimal_dates.values())
    def x_to_year(x):
        return max_date - (x_max - x)
else:
    max_date = None
    def x_to_year(x):
        return x


n_genes  = len(GENES)
cell_w   = x_span * 0.045
col_gap  = x_span * 0.010
heat_gap = x_span * 0.025   # gap between rightmost tip and first column
half_h   = 0.42

# x position of each heatmap column
col_x = [
    x_max + heat_gap + i * (cell_w + col_gap)
    for i in range(n_genes)
]
right_edge = col_x[-1] + cell_w

# Widen figure slightly for 3 genes
fig_w = 7.0 + n_genes * 0.4


y_top = max(t.y for t in tips)

fig, ax = plt.subplots(figsize=(fig_w, 6), facecolor="w")
ref_tree.plotTree(ax, colour="#052049")

all_markers = []
for gene in GENES:
    label = gene["label"]
    color = gene["color"]
    for t in gene_transitions[label]:
        symbol = GAIN_MARKER if t["type"] == "gain" else LOSS_MARKER
        all_markers.append((t["x"], t["y"], symbol, color))

scatter_markers_with_offset(ax, all_markers)

for i, gene in enumerate(GENES):
    label = gene["label"]
    cx    = col_x[i]

    ax.text(cx + cell_w / 2, y_top + 1.5,
            label, rotation=45, ha="left", va="bottom",
            fontsize=8, color="#052049", fontweight="bold")

    for leaf in tips:
        state = state_maps[label].get(id(leaf))
        fill  = (MISSING_COLOR if state is None
                 else PRESENT_COLOR if state == 1
                 else ABSENT_COLOR)
        ax.add_patch(mpatches.Rectangle(
            (cx, leaf.y - half_h), cell_w, half_h * 2,
            linewidth=0.3, edgecolor=BORDER_COLOR, facecolor=fill,
        ))

year_min = x_to_year(x_min)
year_max = x_to_year(x_max)
span     = year_max - year_min
step     = max(1, round(span / 5))
tick_years = list(range(math.ceil(year_min), math.floor(year_max) + 1, step))
if not tick_years:
    tick_years = [round(year_min), round(year_max)]

if max_date is not None:
    tick_xs = [x_max - (max_date - y) for y in tick_years]
else:
    tick_xs = tick_years

ax.set_xticks(tick_xs)
ax.set_xticklabels([str(int(y)) for y in tick_years], fontsize=8, color="#052049")
ax.tick_params(axis='x', length=3, color="#AAAAAA")

ax.set_xlim(x_min - x_span * 0.06, right_edge + cell_w * 3.5)
ax.set_ylim(-1, y_top + 3.5)
ax.set_yticks([])
for loc in ["left", "right", "top"]:
    ax.spines[loc].set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.spines["bottom"].set_color("#AAAAAA")

heatmap_handles = [
    mpatches.Patch(facecolor=PRESENT_COLOR, edgecolor=BORDER_COLOR, label="Present"),
    mpatches.Patch(facecolor=ABSENT_COLOR,  edgecolor=BORDER_COLOR, label="Absent"),
    mpatches.Patch(facecolor=MISSING_COLOR, edgecolor=BORDER_COLOR, label="No data"),
]

transition_handles = []
for gene in GENES:
    label = gene["label"]
    color = gene["color"]
    transition_handles += [
        Line2D([0], [0], marker=GAIN_MARKER, color='w',
               markerfacecolor=color, markersize=8,
               label=f"{label} gain", linewidth=0),
        Line2D([0], [0], marker=LOSS_MARKER, color='w',
               markerfacecolor=color, markersize=8,
               label=f"{label} loss", linewidth=0),
    ]

leg1 = ax.legend(handles=heatmap_handles,
                 title="Gene at tips",
                 bbox_to_anchor=(1.01, 0.75), loc="upper left",
                 frameon=False, fontsize=7, title_fontsize=8,
                 borderaxespad=0)
ax.add_artist(leg1)

ax.legend(handles=transition_handles,
          title="Inferred transitions",
          bbox_to_anchor=(1.01, 0.55), loc="upper left",
          frameon=False, fontsize=7, title_fontsize=8,
          borderaxespad=0)

plt.tight_layout()
plt.savefig(OUTPUT, dpi=300, bbox_inches="tight")
print(f"Saved → {OUTPUT}")
plt.show()
