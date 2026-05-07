#!/usr/bin/env python3
"""Generate a BEAST 2 XML for joint AMR ancestral-state + DTA analysis."""

import argparse
import os
import re
import sys


def normalize_gene_id(name, maxlen=30):
    s = re.sub(r"[()']", '', name)
    s = re.sub(r'[^A-Za-z0-9]+', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    return s[:maxlen]


def parse_fasta(path):
    records = []
    header, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(parts)))
                header = line[1:]
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        records.append((header, ''.join(parts)))
    return records


def parse_gene_map(path):
    gene_map = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    gene_map[int(parts[0])] = parts[1]
                except ValueError:
                    pass  # skip header row
    return gene_map


def parse_locations(path):
    locs = {}
    first = True
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if first:
                first = False
                # skip header if first column doesn't look like a taxon name
                parts = line.split('\t')
                if '|' not in parts[0]:
                    continue
            parts = line.split('\t')
            if len(parts) >= 2:
                locs[parts[0]] = parts[1]
    return locs


def extract_date(taxon):
    if '|' in taxon:
        return taxon.split('|', 1)[1]
    return None


def build_xml(aln_name, aln_records, gene_map, amr_records, locations,
              trait_name, chain_length, log_every, tree_log_every):

    taxa = [h.split()[0] for h, _ in aln_records]

    # --- Per-gene data (filter invariant) ---
    gene_order = []
    all_gene_data = {}
    for col_idx in sorted(gene_map.keys()):
        raw_name = gene_map[col_idx]
        gene_id  = normalize_gene_id(raw_name)
        taxon_vals = {}
        for header, seq in amr_records:
            taxon = header.split()[0]
            if col_idx < len(seq) and seq[col_idx] in ('0', '1'):
                taxon_vals[taxon] = seq[col_idx]
        unique = {v for v in taxon_vals.values()}
        if len(unique) < 2:
            print(f'  Skipping {raw_name} (invariant: {unique or "no data"})', file=sys.stderr)
            continue
        all_gene_data[gene_id] = taxon_vals
        gene_order.append(gene_id)

    # --- Location setup ---
    loc_values = sorted(set(locations.values()))
    n_locs = len(loc_values)
    # Dimension for rateIndicator and relativeGeoRates: n*(n-1)/2 (upper triangle of rate matrix)
    bssvs_dim = n_locs * (n_locs - 1) // 2
    code_map = ','.join(f'{v}={i}' for i, v in enumerate(loc_values))
    code_map += ',? = ' + ' '.join(str(i) for i in range(n_locs)) + ' '
    loc_trait_str = ',\n\n'.join(f'{t}={locations[t]}' for t in taxa if t in locations)

    # --- Date trait string ---
    dates_str = ','.join(f'{t}={extract_date(t)}' for t in taxa if extract_date(t))

    L = []  # output lines

    # =========================================================================
    # Header
    # =========================================================================
    L.append('<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
             '<beast beautitemplate=\'Standard\' beautistatus=\'\' '
             'namespace="beast.core:beast.evolution.alignment:'
             'beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:'
             'beast.evolution.operators:beast.evolution.sitemodel:'
             'beast.evolution.substitutionmodel:beast.base.evolution.alignment:'
             'beast.pkgmgmt:beast.base.core:beast.base.inference:'
             'beast.base.evolution.tree.coalescent:beast.pkgmgmt:beast.base.core:'
             'beast.base.inference.util:beast.evolution.nuc:'
             'beast.base.evolution.operator:beast.base.inference.operator:'
             'beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:'
             'beast.base.evolution.likelihood" '
             'required="BEAST.base v2.7.8:BICEPS v1.1.2:MM v1.2.1:BEAST_CLASSIC v1.6.4" '
             'version="2.7">')
    L.append('')

    # =========================================================================
    # Main alignment
    # =========================================================================
    L.append(f'    <data')
    L.append(f'id="{aln_name}"')
    L.append(f'spec="Alignment"')
    L.append(f'name="alignment">')
    for header, seq in aln_records:
        taxon = header.split()[0]
        L.append(f'        <sequence id="seq_{taxon}1" spec="Sequence" '
                 f'taxon="{taxon}" totalcount="4" value="{seq}"/>')
    L.append('    </data>')
    L.append('')

    # =========================================================================
    # AMR binary data blocks
    # =========================================================================
    for gene_id in gene_order:
        taxon_vals = all_gene_data[gene_id]
        L.append(f'    <data id="amr_{gene_id}" spec="Alignment" dataType="binary">')
        for header, _ in amr_records:
            taxon = header.split()[0]
            val = taxon_vals.get(taxon, '?')
            L.append(f'        <sequence id="seq_{taxon}_{gene_id}" spec="Sequence" '
                     f'taxon="{taxon}" totalcount="2" value="{val}"/>')
        L.append('    </data>')
        L.append('')

    # =========================================================================
    # Map elements
    # =========================================================================
    for name, cls in [
        ('Uniform',            'beast.base.inference.distribution.Uniform'),
        ('Exponential',        'beast.base.inference.distribution.Exponential'),
        ('LogNormal',          'beast.base.inference.distribution.LogNormalDistributionModel'),
        ('Normal',             'beast.base.inference.distribution.Normal'),
        ('Beta',               'beast.base.inference.distribution.Beta'),
        ('Gamma',              'beast.base.inference.distribution.Gamma'),
        ('LaplaceDistribution','beast.base.inference.distribution.LaplaceDistribution'),
        ('prior',              'beast.base.inference.distribution.Prior'),
        ('InverseGamma',       'beast.base.inference.distribution.InverseGamma'),
        ('OneOnX',             'beast.base.inference.distribution.OneOnX'),
    ]:
        L.append(f'    <map name="{name}" >{cls}</map>')
    L.append('')

    # =========================================================================
    # MCMC run
    # =========================================================================
    L.append(f'    <run id="mcmc" spec="MCMC" chainLength="{chain_length}">')

    # ----- State -----
    L.append(f'        <state id="state" spec="State" storeEvery="5000">')
    L.append(f'            <tree id="Tree.t:{aln_name}" spec="beast.base.evolution.tree.Tree" name="stateNode">')
    L.append(f'                <trait id="dateTrait.t:{aln_name}" spec="beast.base.evolution.tree.TraitSet" '
             f'dateFormat="yyyy-M-dd" traitname="date" value="{dates_str}">')
    L.append(f'                    <taxa id="TaxonSet.{aln_name}" spec="TaxonSet">')
    L.append(f'                        <alignment idref="{aln_name}"/>')
    L.append(f'                    </taxa>')
    L.append(f'                </trait>')
    L.append(f'                <taxonset idref="TaxonSet.{aln_name}"/>')
    L.append(f'            </tree>')
    # Explicit dimension = n*(n-1)/2 for the BSSVS rate indicator and rate vectors
    L.append(f'            <stateNode id="rateIndicator.s:{trait_name}" spec="parameter.BooleanParameter" dimension="{bssvs_dim}">true</stateNode>')
    L.append(f'            <parameter id="relativeGeoRates.s:{trait_name}" spec="parameter.RealParameter" dimension="{bssvs_dim}" name="stateNode">1.0</parameter>')
    L.append(f'            <parameter id="traitClockRate.c:{trait_name}" spec="parameter.RealParameter" name="stateNode">1.0</parameter>')
    L.append(f'            <parameter id="clockRate.c:{aln_name}" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>')
    for rate in ('rateAC', 'rateAG', 'rateAT', 'rateCG', 'rateGT'):
        L.append(f'            <parameter id="{rate}.s:{aln_name}" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>')
    L.append(f'            <parameter id="freqParameter.s:{aln_name}" spec="parameter.RealParameter" '
             f'dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>')
    L.append(f'            <parameter id="gammaShape.s:{aln_name}" spec="parameter.RealParameter" lower="0.1" name="stateNode">1.0</parameter>')
    L.append(f'            <parameter id="BICEPSPopSize.t:{aln_name}" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>')
    L.append(f'        </state>')

    # ----- Init -----
    L.append(f'        <init id="RandomTree.t:{aln_name}" spec="RandomTree" estimate="false" '
             f'initial="@Tree.t:{aln_name}" taxa="@{aln_name}">')
    L.append(f'            <populationModel id="ConstantPopulation0.t:{aln_name}" spec="ConstantPopulation">')
    L.append(f'                <parameter id="randomPopSize.t:{aln_name}" spec="parameter.RealParameter" name="popSize">1.0</parameter>')
    L.append(f'            </populationModel>')
    L.append(f'        </init>')

    # ----- Posterior -----
    L.append(f'        <distribution id="posterior" spec="CompoundDistribution">')
    L.append(f'            <distribution id="prior" spec="CompoundDistribution">')

    # BICEPS
    L.append(f'                <distribution id="BICEPS.t:{aln_name}" spec="biceps.BICEPS" '
             f'linkedMean="true" ploidy="1.0" populationMean="@BICEPSPopSize.t:{aln_name}">')
    L.append(f'                    <parameter id="RealParameter.21" spec="parameter.RealParameter" '
             f'lower="0.0" name="populationShape" upper="0.0">3.0</parameter>')
    L.append(f'                    <groupSizes id="BICEPSGroupSizes.t:{aln_name}" spec="parameter.IntegerParameter" '
             f'dimension="10" estimate="false">1</groupSizes>')
    L.append(f'                    <treeIntervals id="BICEPSTreeIntervals.t:{aln_name}" '
             f'spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree.t:{aln_name}"/>')
    L.append(f'                </distribution>')

    # BICEPSPopSize prior
    L.append(f'                <prior id="BICEPSPopSizePrior.t:{aln_name}" name="distribution" x="@BICEPSPopSize.t:{aln_name}">')
    L.append(f'                    <LogNormal id="LogNormalDistributionModel.0" meanInRealSpace="true" name="distr">')
    L.append(f'                        <M id="Function$Constant.0" spec="Function$Constant" value="1.0"/>')
    L.append(f'                        <S id="Function$Constant.1" spec="Function$Constant" value="1.0"/>')
    L.append(f'                    </LogNormal>')
    L.append(f'                </prior>')

    # Clock prior
    L.append(f'                <prior id="ClockPrior.c:{aln_name}" name="distribution" x="@clockRate.c:{aln_name}">')
    L.append(f'                    <Uniform id="Uniform.0" name="distr" upper="Infinity"/>')
    L.append(f'                </prior>')

    # Freq prior
    L.append(f'                <prior id="FrequenciesPrior.s:{aln_name}" name="distribution" x="@freqParameter.s:{aln_name}">')
    L.append(f'                    <distr id="Dirichlet.0" spec="distribution.Dirichlet">')
    L.append(f'                        <parameter id="RealParameter.20" spec="parameter.RealParameter" '
             f'dimension="4" estimate="false" name="alpha">4.0 4.0 4.0 4.0</parameter>')
    L.append(f'                    </distr>')
    L.append(f'                </prior>')

    # GammaShape prior
    L.append(f'                <prior id="GammaShapePrior.s:{aln_name}" name="distribution" x="@gammaShape.s:{aln_name}">')
    L.append(f'                    <Exponential id="Exponential.0" name="distr">')
    L.append(f'                        <parameter id="RealParameter.0" spec="parameter.RealParameter" estimate="false" name="mean">1.0</parameter>')
    L.append(f'                    </Exponential>')
    L.append(f'                </prior>')

    # nonZeroRatePrior (BSSVS Poisson)
    L.append(f'                <prior id="nonZeroRatePrior.s:{trait_name}" name="distribution">')
    L.append(f'                    <x id="nonZeroRates.s:{trait_name}" spec="beast.base.evolution.Sum">')
    L.append(f'                        <arg idref="rateIndicator.s:{trait_name}"/>')
    L.append(f'                    </x>')
    L.append(f'                    <distr id="Poisson.0" spec="distribution.Poisson" offset="1.0">')
    L.append(f'                        <parameter id="RealParameter.7" spec="parameter.RealParameter" estimate="false" name="lambda">0.693</parameter>')
    L.append(f'                    </distr>')
    L.append(f'                </prior>')

    # GTR rate priors
    for rate, gid, rp_a, rp_b, alpha, beta in [
        ('rateAC', 'Gamma.2', 'RealParameter.8',  'RealParameter.9',  '0.05', '10.0'),
        ('rateAG', 'Gamma.3', 'RealParameter.10', 'RealParameter.11', '0.05', '20.0'),
        ('rateAT', 'Gamma.4', 'RealParameter.12', 'RealParameter.13', '0.05', '10.0'),
        ('rateCG', 'Gamma.5', 'RealParameter.14', 'RealParameter.15', '0.05', '10.0'),
        ('rateGT', 'Gamma.7', 'RealParameter.18', 'RealParameter.19', '0.05', '10.0'),
    ]:
        suffix = rate[4:].upper()
        L.append(f'                <prior id="Rate{suffix}Prior.s:{aln_name}" name="distribution" x="@{rate}.s:{aln_name}">')
        L.append(f'                    <Gamma id="{gid}" name="distr">')
        L.append(f'                        <parameter id="{rp_a}" spec="parameter.RealParameter" estimate="false" name="alpha">{alpha}</parameter>')
        L.append(f'                        <parameter id="{rp_b}" spec="parameter.RealParameter" estimate="false" name="beta">{beta}</parameter>')
        L.append(f'                    </Gamma>')
        L.append(f'                </prior>')

    # relativeGeoRates prior
    L.append(f'                <prior id="relativeGeoRatesPrior.s:{trait_name}" name="distribution" x="@relativeGeoRates.s:{trait_name}">')
    L.append(f'                    <Gamma id="Gamma.1" name="distr">')
    L.append(f'                        <parameter id="RealParameter.5" spec="parameter.RealParameter" estimate="false" name="alpha">1.0</parameter>')
    L.append(f'                        <parameter id="RealParameter.6" spec="parameter.RealParameter" estimate="false" name="beta">1.0</parameter>')
    L.append(f'                    </Gamma>')
    L.append(f'                </prior>')

    # geoclockPrior
    L.append(f'                <prior id="geoclockPrior.c:{trait_name}" name="distribution" x="@traitClockRate.c:{trait_name}">')
    L.append(f'                    <Gamma id="Gamma.0" name="distr">')
    L.append(f'                        <parameter id="RealParameter.3" spec="parameter.RealParameter" estimate="false" name="alpha">0.001</parameter>')
    L.append(f'                        <parameter id="RealParameter.4" spec="parameter.RealParameter" estimate="false" name="beta">1000.0</parameter>')
    L.append(f'                    </Gamma>')
    L.append(f'                </prior>')

    L.append(f'            </distribution>')  # close prior

    # ----- Likelihood -----
    L.append(f'            <distribution id="likelihood" spec="CompoundDistribution" useThreads="true">')

    # GTR tree likelihood
    L.append(f'                <distribution id="treeLikelihood.{aln_name}" '
             f'spec="ThreadedTreeLikelihood" data="@{aln_name}" tree="@Tree.t:{aln_name}">')
    L.append(f'                    <siteModel id="SiteModel.s:{aln_name}" spec="SiteModel" '
             f'gammaCategoryCount="4" shape="@gammaShape.s:{aln_name}">')
    L.append(f'                        <parameter id="mutationRate.s:{aln_name}" spec="parameter.RealParameter" '
             f'estimate="false" lower="0.0" name="mutationRate">1.0</parameter>')
    L.append(f'                        <parameter id="proportionInvariant.s:{aln_name}" spec="parameter.RealParameter" '
             f'estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>')
    L.append(f'                        <substModel id="gtr.s:{aln_name}" spec="GTR" '
             f'rateAC="@rateAC.s:{aln_name}" rateAG="@rateAG.s:{aln_name}" '
             f'rateAT="@rateAT.s:{aln_name}" rateCG="@rateCG.s:{aln_name}" rateGT="@rateGT.s:{aln_name}">')
    L.append(f'                            <parameter id="rateCT.s:{aln_name}" spec="parameter.RealParameter" '
             f'estimate="false" lower="0.0" name="rateCT">1.0</parameter>')
    L.append(f'                            <frequencies id="estimatedFreqs.s:{aln_name}" spec="Frequencies" '
             f'frequencies="@freqParameter.s:{aln_name}"/>')
    L.append(f'                        </substModel>')
    L.append(f'                    </siteModel>')
    L.append(f'                    <branchRateModel id="StrictClock.c:{aln_name}" '
             f'spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:{aln_name}"/>')
    L.append(f'                </distribution>')

    # Per-gene AMR likelihoods (LewisMK binary)
    for gene_id in gene_order:
        L.append(f'                <distribution id="traitedtreeLikelihood.{gene_id}" '
                 f'spec="beastclassic.evolution.likelihood.AncestralStateTreeLikelihood" '
                 f'tag="{gene_id}" tree="@Tree.t:{aln_name}">')
        L.append(f'                    <data idref="amr_{gene_id}"/>')
        L.append(f'                    <siteModel id="SiteModel.s:{gene_id}" spec="SiteModel">')
        L.append(f'                        <parameter id="mutationRate.s:{gene_id}" spec="parameter.RealParameter" '
                 f'estimate="false" lower="0.0" name="mutationRate">1.0</parameter>')
        L.append(f'                        <substModel id="LewisMK.s:{gene_id}" '
                 f'spec="morphmodels.evolution.substitutionmodel.LewisMK" stateNumber="2"/>')
        L.append(f'                    </siteModel>')
        L.append(f'                    <branchRateModel id="StrictClockModel.c:{gene_id}" '
                 f'spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="1.0"/>')
        L.append(f'                </distribution>')

    # DTA (location) likelihood
    L.append(f'                <distribution id="traitedtreeLikelihood.{trait_name}" '
             f'spec="beastclassic.evolution.likelihood.AncestralStateTreeLikelihood" '
             f'tag="location" tree="@Tree.t:{aln_name}">')
    L.append(f'                    <data')
    L.append(f'id="{trait_name}"')
    L.append(f'spec="beastclassic.evolution.alignment.AlignmentFromTrait">')
    L.append(f'                        <traitSet id="traitSet.{trait_name}" '
             f'spec="beast.base.evolution.tree.TraitSet" '
             f'taxa="@TaxonSet.{aln_name}" traitname="discrete">{loc_trait_str}</traitSet>')
    L.append(f'                        <userDataType id="traitDataType.{trait_name}" '
             f'spec="beast.base.evolution.datatype.UserDataType" '
             f'codeMap="{code_map}" codelength="-1" states="{n_locs}"/>')
    L.append(f'                    </data>')
    L.append(f'                    <siteModel id="geoSiteModel.s:{trait_name}" spec="SiteModel" gammaCategoryCount="1">')
    L.append(f'                        <parameter id="mutationRate.s:{trait_name}" spec="parameter.RealParameter" '
             f'estimate="false" name="mutationRate">1.0</parameter>')
    L.append(f'                        <parameter id="gammaShape.s:{trait_name}" spec="parameter.RealParameter" '
             f'estimate="false" name="shape">1.0</parameter>')
    L.append(f'                        <parameter id="proportionInvariant.s:{trait_name}" spec="parameter.RealParameter" '
             f'estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>')
    L.append(f'                        <substModel id="svs.s:{trait_name}" '
             f'spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel" '
             f'rateIndicator="@rateIndicator.s:{trait_name}" rates="@relativeGeoRates.s:{trait_name}">')
    uniform_freq = ' '.join([f'{1.0/n_locs:.10f}'] * n_locs)
    L.append(f'                            <frequencies id="traitfreqs.s:{trait_name}" spec="Frequencies">')
    L.append(f'                                <parameter id="traitfrequencies.s:{trait_name}" '
             f'spec="parameter.RealParameter" dimension="{n_locs}" name="frequencies">{uniform_freq}</parameter>')
    L.append(f'                            </frequencies>')
    L.append(f'                        </substModel>')
    L.append(f'                    </siteModel>')
    L.append(f'                    <branchRateModel id="StrictClockModel.c:{trait_name}" '
             f'spec="beast.base.evolution.branchratemodel.StrictClockModel" '
             f'clock.rate="@traitClockRate.c:{trait_name}"/>')
    L.append(f'                </distribution>')

    L.append(f'            </distribution>')  # close likelihood
    L.append(f'            <distribution id="fossilCalibrations" spec="CompoundDistribution"/>')
    L.append(f'        </distribution>')  # close posterior

    # =========================================================================
    # Operators
    # =========================================================================
    # Geo operators
    L.append(f'        <operator id="georateScaler.s:{trait_name}" spec="ScaleOperator" '
             f'parameter="@relativeGeoRates.s:{trait_name}" scaleAllIndependently="true" scaleFactor="0.99" weight="30.0"/>')
    L.append(f'        <operator id="indicatorFlip.s:{trait_name}" spec="operator.BitFlipOperator" '
             f'parameter="@rateIndicator.s:{trait_name}" weight="30.0"/>')
    L.append(f'        <operator id="geoMuScaler.c:{trait_name}" spec="ScaleOperator" '
             f'parameter="@traitClockRate.c:{trait_name}" scaleFactor="0.9" weight="3.0"/>')
    L.append(f'        <operator id="BSSVSoperator.c:{trait_name}" '
             f'spec="beastclassic.evolution.operators.BitFlipBSSVSOperator" '
             f'indicator="@rateIndicator.s:{trait_name}" mu="@traitClockRate.c:{trait_name}" weight="30.0"/>')

    # AVMN + clock up-down
    L.append(f'        <operator id="StrictClockRateScaler.c:{aln_name}" spec="AdaptableOperatorSampler" weight="1.5">')
    L.append(f'            <parameter idref="clockRate.c:{aln_name}"/>')
    L.append(f'            <operator id="AVMNOperator.{aln_name}" '
             f'spec="kernel.AdaptableVarianceMultivariateNormalOperator" '
             f'allowNonsense="true" beta="0.05" burnin="400" initial="800" weight="0.1">')
    L.append(f'                <transformations id="AVMNSumTransform.{aln_name}" '
             f'spec="operator.kernel.Transform$LogConstrainedSumTransform">')
    L.append(f'                    <f idref="freqParameter.s:{aln_name}"/>')
    L.append(f'                </transformations>')
    L.append(f'                <transformations id="AVMNLogTransform.{aln_name}" '
             f'spec="operator.kernel.Transform$LogTransform">')
    for p in ('clockRate', 'rateAC', 'rateAG', 'rateAT', 'rateCG', 'rateGT', 'gammaShape'):
        suffix = 'c' if p == 'clockRate' else 's'
        L.append(f'                    <f idref="{p}.{suffix}:{aln_name}"/>')
    L.append(f'                </transformations>')
    L.append(f'                <transformations id="AVMNNoTransform.{aln_name}" '
             f'spec="operator.kernel.Transform$NoTransform">')
    L.append(f'                    <f idref="Tree.t:{aln_name}"/>')
    L.append(f'                </transformations>')
    L.append(f'            </operator>')
    L.append(f'            <operator id="StrictClockRateScalerX.c:{aln_name}" '
             f'spec="kernel.BactrianScaleOperator" parameter="@clockRate.c:{aln_name}" upper="10.0" weight="3.0"/>')
    L.append(f'        </operator>')

    L.append(f'        <operator id="strictClockUpDownOperator.c:{aln_name}" spec="AdaptableOperatorSampler" weight="1.5">')
    L.append(f'            <parameter idref="clockRate.c:{aln_name}"/>')
    L.append(f'            <tree idref="Tree.t:{aln_name}"/>')
    L.append(f'            <operator idref="AVMNOperator.{aln_name}"/>')
    L.append(f'            <operator id="strictClockUpDownOperatorX.c:{aln_name}" '
             f'spec="operator.kernel.BactrianUpDownOperator" scaleFactor="0.75" weight="3.0">')
    L.append(f'                <up idref="clockRate.c:{aln_name}"/>')
    L.append(f'                <down idref="Tree.t:{aln_name}"/>')
    L.append(f'            </operator>')
    L.append(f'        </operator>')

    # GTR rate operators
    for rate, suffix in [('rateAC','AC'),('rateAG','AG'),('rateAT','AT'),('rateCG','CG'),('rateGT','GT')]:
        L.append(f'        <operator id="Rate{suffix}Scaler.s:{aln_name}" spec="AdaptableOperatorSampler" weight="0.05">')
        L.append(f'            <parameter idref="{rate}.s:{aln_name}"/>')
        L.append(f'            <operator idref="AVMNOperator.{aln_name}"/>')
        L.append(f'            <operator id="Rate{suffix}ScalerX.s:{aln_name}" '
                 f'spec="kernel.BactrianScaleOperator" parameter="@{rate}.s:{aln_name}" scaleFactor="0.5" upper="10.0" weight="0.1"/>')
        L.append(f'        </operator>')

    # Freq exchanger
    L.append(f'        <operator id="FrequenciesExchanger.s:{aln_name}" spec="AdaptableOperatorSampler" weight="0.05">')
    L.append(f'            <parameter idref="freqParameter.s:{aln_name}"/>')
    L.append(f'            <operator idref="AVMNOperator.{aln_name}"/>')
    L.append(f'            <operator id="FrequenciesExchangerX.s:{aln_name}" '
             f'spec="operator.kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1">')
    L.append(f'                <parameter idref="freqParameter.s:{aln_name}"/>')
    L.append(f'            </operator>')
    L.append(f'        </operator>')

    # GammaShape scaler
    L.append(f'        <operator id="gammaShapeScaler.s:{aln_name}" spec="AdaptableOperatorSampler" weight="0.05">')
    L.append(f'            <parameter idref="gammaShape.s:{aln_name}"/>')
    L.append(f'            <operator idref="AVMNOperator.{aln_name}"/>')
    L.append(f'            <operator id="gammaShapeScalerX.s:{aln_name}" '
             f'spec="kernel.BactrianScaleOperator" parameter="@gammaShape.s:{aln_name}" scaleFactor="0.5" upper="10.0" weight="0.1"/>')
    L.append(f'        </operator>')

    # BICEPS tree operators
    L.append(f'        <operator id="BICEPSEpochTop.t:{aln_name}" spec="EpochFlexOperator" scaleFactor="0.1" tree="@Tree.t:{aln_name}" weight="2.0"/>')
    L.append(f'        <operator id="BICEPSEpochAll.t:{aln_name}" spec="EpochFlexOperator" fromOldestTipOnly="false" scaleFactor="0.1" tree="@Tree.t:{aln_name}" weight="2.0"/>')
    L.append(f'        <operator id="BICEPSTreeFlex.t:{aln_name}" spec="TreeStretchOperator" scaleFactor="0.01" tree="@Tree.t:{aln_name}" weight="2.0"/>')
    L.append(f'        <operator id="BICEPSTreeRootScaler.t:{aln_name}" spec="kernel.BactrianScaleOperator" rootOnly="true" scaleFactor="0.1" tree="@Tree.t:{aln_name}" upper="10.0" weight="3.0"/>')
    L.append(f'        <operator id="BICEPSUniformOperator.t:{aln_name}" spec="kernel.BactrianNodeOperator" tree="@Tree.t:{aln_name}" weight="30.0"/>')
    L.append(f'        <operator id="BICEPSSubtreeSlide.t:{aln_name}" spec="kernel.BactrianSubtreeSlide" tree="@Tree.t:{aln_name}" weight="15.0"/>')
    L.append(f'        <operator id="BICEPSNarrow.t:{aln_name}" spec="Exchange" tree="@Tree.t:{aln_name}" weight="15.0"/>')
    L.append(f'        <operator id="BICEPSWide.t:{aln_name}" spec="Exchange" isNarrow="false" tree="@Tree.t:{aln_name}" weight="3.0"/>')
    L.append(f'        <operator id="BICEPSWilsonBalding.t:{aln_name}" spec="WilsonBalding" tree="@Tree.t:{aln_name}" weight="3.0"/>')
    L.append(f'        <operator id="BICEPSPopSizesScaler.t:{aln_name}" spec="kernel.BactrianScaleOperator" parameter="@BICEPSPopSize.t:{aln_name}" scaleFactor="0.1" upper="10.0" weight="5.0"/>')

    # =========================================================================
    # Loggers
    # =========================================================================
    L.append(f'        <logger id="tracelog" spec="Logger" fileName="$(filebase).log" '
             f'logEvery="{log_every}" model="@posterior" sanitiseHeaders="true" sort="smart">')
    L.append(f'            <log idref="posterior"/>')
    L.append(f'            <log idref="likelihood"/>')
    L.append(f'            <log idref="prior"/>')
    L.append(f'            <log idref="treeLikelihood.{aln_name}"/>')
    L.append(f'            <log id="TreeHeight.t:{aln_name}" spec="beast.base.evolution.tree.TreeStatLogger" tree="@Tree.t:{aln_name}"/>')
    L.append(f'            <log idref="rateIndicator.s:{trait_name}"/>')
    L.append(f'            <log idref="relativeGeoRates.s:{trait_name}"/>')
    L.append(f'            <log idref="traitClockRate.c:{trait_name}"/>')
    L.append(f'            <log id="geoSubstModelLogger.s:{trait_name}" '
             f'spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger" '
             f'dataType="@traitDataType.{trait_name}" model="@svs.s:{trait_name}"/>')
    L.append(f'            <log idref="clockRate.c:{aln_name}"/>')
    for rate in ('rateAC', 'rateAG', 'rateAT', 'rateCG', 'rateGT'):
        L.append(f'            <log idref="{rate}.s:{aln_name}"/>')
    L.append(f'            <log idref="freqParameter.s:{aln_name}"/>')
    L.append(f'            <log idref="gammaShape.s:{aln_name}"/>')
    L.append(f'            <log idref="BICEPS.t:{aln_name}"/>')
    L.append(f'            <log idref="BICEPSPopSize.t:{aln_name}"/>')
    L.append(f'        </logger>')

    L.append(f'        <logger id="screenlog" spec="Logger" logEvery="{log_every}">')
    L.append(f'            <log idref="posterior"/>')
    L.append(f'            <log idref="likelihood"/>')
    L.append(f'            <log idref="prior"/>')
    L.append(f'        </logger>')

    L.append(f'        <logger id="treelog.t:{aln_name}" spec="Logger" '
             f'fileName="$(filebase)-$(tree).trees" logEvery="{tree_log_every}" mode="tree">')
    L.append(f'            <log id="TreeWithMetaDataLogger.t:{aln_name}" '
             f'spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@Tree.t:{aln_name}"/>')
    L.append(f'        </logger>')

    # DTA tree-with-trait logger
    L.append(f'        <logger id="treeWithTraitLogger.{trait_name}" spec="Logger" '
             f'fileName="{trait_name}_tree_with_trait.trees" logEvery="{tree_log_every}" mode="tree">')
    L.append(f'            <log id="treeWithTraitLoggerItem.t:{aln_name}" '
             f'spec="beastclassic.evolution.tree.TreeWithTraitLogger" tree="@Tree.t:{aln_name}">')
    L.append(f'                <metadata idref="posterior"/>')
    L.append(f'                <metadata idref="traitedtreeLikelihood.{trait_name}"/>')
    L.append(f'            </log>')
    L.append(f'        </logger>')

    # Per-gene AMR tree loggers
    for gene_id in gene_order:
        L.append(f'        <logger id="treeLogger.{gene_id}" spec="Logger" '
                 f'fileName="{gene_id}.trees" logEvery="{tree_log_every}" mode="tree">')
        L.append(f'            <log id="treeWithTraitLoggerItem.{gene_id}" '
                 f'spec="beastclassic.evolution.tree.TreeWithTraitLogger" tree="@Tree.t:{aln_name}">')
        L.append(f'                <metadata idref="posterior"/>')
        L.append(f'                <metadata idref="traitedtreeLikelihood.{gene_id}"/>')
        L.append(f'            </log>')
        L.append(f'        </logger>')

    L.append(f'        <operatorschedule id="OperatorSchedule" spec="OperatorSchedule"/>')
    L.append(f'    </run>')
    L.append('</beast>')

    return '\n'.join(L)


def main():
    p = argparse.ArgumentParser(description='Generate BEAST 2 XML for joint AMR + DTA analysis.')
    p.add_argument('--alignment',      required=True, help='Main alignment FASTA')
    p.add_argument('--gene-map',       required=True, help='AMR gene map TSV (col_index\\tgene_name)')
    p.add_argument('--amr-traits',     required=True, help='AMR binary traits FASTA')
    p.add_argument('--locations',      required=True, help='Locations TSV (taxon\\tlocation)')
    p.add_argument('--output-dir',     default='xmls', help='Output directory (default: xmls)')
    p.add_argument('--location-trait', default='hhs',  help='Trait name for DTA (default: hhs)')
    p.add_argument('--chain-length',   type=int, default=10_000_000)
    p.add_argument('--log-every',      type=int, default=1_000)
    p.add_argument('--tree-log-every', type=int, default=10_000)
    args = p.parse_args()

    aln_name = os.path.splitext(os.path.basename(args.alignment))[0]

    print(f'Alignment:   {args.alignment}  →  id="{aln_name}"')
    aln_records = parse_fasta(args.alignment)
    print(f'  {len(aln_records)} taxa')

    print(f'Gene map:    {args.gene_map}')
    gene_map = parse_gene_map(args.gene_map)
    print(f'  {len(gene_map)} columns')

    print(f'AMR traits:  {args.amr_traits}')
    amr_records = parse_fasta(args.amr_traits)

    print(f'Locations:   {args.locations}')
    locations = parse_locations(args.locations)
    loc_vals = sorted(set(locations.values()))
    print(f'  {len(loc_vals)} unique locations: {loc_vals}')

    os.makedirs(args.output_dir, exist_ok=True)

    xml = build_xml(
        aln_name=aln_name,
        aln_records=aln_records,
        gene_map=gene_map,
        amr_records=amr_records,
        locations=locations,
        trait_name=args.location_trait,
        chain_length=args.chain_length,
        log_every=args.log_every,
        tree_log_every=args.tree_log_every,
    )

    out_path = os.path.join(args.output_dir, f'{aln_name}.xml')
    with open(out_path, 'w') as f:
        f.write(xml)
    print(f'Written → {out_path}')


if __name__ == '__main__':
    main()