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
* [ ] Core-genome alignment and phylogenetic inference per ST (Gubbins)
* [ ] Ancestral state reconstruction of AMR gene gain/loss events (DTA)
