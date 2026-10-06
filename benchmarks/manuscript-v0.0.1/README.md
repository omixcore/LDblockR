# Quantitative benchmark used in the LDblockR 0.0.1 manuscript

## Scope and environment

These are measured results from LDblockR 0.0.1, the version reported in the manuscript. They are historical, version-specific results and should not be read as a performance benchmark of the current development branch.

- Ubuntu 24.04.3 LTS, x86_64; AMD EPYC 9V74 allocation (9 vCPUs); 8 GiB container memory limit
- R 4.3.3; Ubuntu reference `libblas`
- Three independent R processes per configuration; values in the manuscript are medians
- Synthetic diploid dosages with allele frequencies drawn uniformly from 0.05–0.50 and 1% missing calls
- Dosage-`r²` calculation and heatmap rendering; GNU `/usr/bin/time -v` recorded process peak RSS

## Marker scaling

At 281 samples, median dosage-`r²` calculation time ranged from 0.010 s (100 markers) to 8.138 s (3,093 markers); median peak RSS ranged from 79.7 to 479.9 MiB.

## Sample scaling

At 500 markers, median dosage-`r²` calculation time ranged from 0.084 s (100 samples) to 0.696 s (1,000 samples); median peak RSS ranged from 132.3 to 134.0 MiB.

## Reproducibility files

- [Marker-scaling replicate measurements](marker_scaling.tsv)
- [Sample-scaling replicate measurements](sample_scaling.tsv)
- [Benchmark R script](benchmark_regions.R)

The TSV files retain all three replicates, including LD computation time, rendering time, output/object sizes, R/package versions, peak RSS, and process wall time. The benchmark covers dosage-`r²` and plotting only; it does not measure D′, block detection, GWAS, or performance against other software. Measurements are machine-dependent.
