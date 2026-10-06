# Mixed-model GWAS validation record

This package implements a P3D-style single-variance-component MLM independently
of GAPIT.  The fitted null model is

`y = X beta + Zu + e`, `u ~ N(0, K sigma_g^2)`, `e ~ N(0, I sigma_e^2)`.

## Computational contract

1. Phenotype rows are matched to genotype Taxa IDs. Repeated rows are retained
   by `replicate = "expand"` or averaged by `replicate = "mean"`.
2. A standardized marker relationship is computed once when `kinship = NULL`;
   an input TASSEL matrix is subset to shared Taxa IDs.
3. The symmetric kinship matrix is eigendecomposed once. A one-dimensional
   REML-like profile over `log(delta)`, `delta = sigma_e^2/sigma_g^2`, is
   optimized on `[-8, 8]`.
4. Each marker is mean-imputed only for the test matrix, whitened with the
   cached eigenvectors, residualized against the fixed effects, and tested by
   a one-degree-of-freedom GLS t statistic.
5. `PVE` is `max(0, 1 - SSE_full/SSE_null)` in the whitened model; `model_PVE`
   is `1/(1 + delta)` from the null model.

The fixed-variance scan is the same P3D principle used for efficient MLM GWAS;
it is not a copy of GAPIT source code and it does not claim a separate REML
optimization for every marker.

## Regional figure correspondence

Use `gwas_ld_region()` with the exact threshold used in the Manhattan/GWAS
analysis. It selects the strongest local cluster of over-threshold markers and
returns both a region string and the regional association table. Pass that
table to `plot_ld()`. The plotting code matches rows by marker ID and then by
`chr + position`; unmatched rows are dropped with a warning. With
`show_connectors = TRUE`, the same LD matrix index is used for the upper point,
connector line, and heatmap-baseline marker, preventing visual misalignment.

## Association diagnostics

`plot_manhattan()` and `plot_qq()` consume the same normalized table returned by
`gwas_mlm()`. The Manhattan function supports chromosome, significance, and PVE
color encodings, marker-level PVE point sizes, Bonferroni/suggestive lines, and
SNP labels. The Q-Q function uses beta order-statistic confidence limits and
reports the 1-df genomic inflation lambda. `save_gwas_plot()` is tested for PDF,
SVG, PNG, and TIFF output and does not depend on an external plotting package.

## Required audit files for submission

Save the exact `gwas_mlm()` call, package version, `sessionInfo()`, the output
TSV, `attr(result, "model")`, the input SHA-256 values in
`extdata/tassel/PROVENANCE.json`, and the complete installation/build/check
logs.  For an external LDBlockShow comparison, retain the 861-pair raw output,
software version, command line, parameters, and the row-wise comparison table;
do not retain only the summary count.
