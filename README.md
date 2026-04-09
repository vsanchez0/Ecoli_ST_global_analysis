# Ecoli_ST_global_analysis

### Overview
This repository corresponds to our project investigating how antimicrobial resistance (AMR) gene contenet shapes global dispersion patterns of *Escherichia coli* sequence types (STs).

### Research Aims
We seek to determine whether *E. coli* STs exhibit consistent, predictable relationships among:
* Genomic traits (AMR gene content, virulence factors)
* Global dispersion patterns (geographic range, introduction events)
* AMR gene gain/loss dynamics over time and across regions

### Methods
**1. Genomic Data Collection**
* Collection of publicly available *E. coli* genome sequences with associated metadata
- [ ] Select four STs with around 10,000 samples per ST
- [ ] Figure out which sequences have corresponding granular time and location data (months, states)
- [ ] Figure out how many of these samples have long read data
- [ ] Download assemblies

**2. AMR Gene Profiling**
* AMR gene detection using CARD/RGI
- [ ] Annotate assemblies

**3. Phylogenetic Reconstruction**
* Core-genome alignmnet and phylogenetic inference per ST (Gubbins)
* Ancestral state reconstruction of AMR gene gain/loss events (DTA)
