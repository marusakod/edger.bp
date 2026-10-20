# Using edger.bp

`edger.bp` runs edgeR-compatible analyses on genes-by-samples count
matrices stored with BPCells. It is intended for pseudobulk matrices
where the number of genes is large enough that row-blocked GLM fitting
and file-backed storage are useful.

The package does not define new statistical methods. It mirrors
supported edgeR workflows for BPCells-backed matrices, so statistical
choices should follow the edgeR documentation.

## Data Contract

Inputs should satisfy the usual edgeR count-matrix requirements:

- rows are genes or features;
- columns are samples or pseudobulk sample-cluster columns;
- entries are non-negative integer counts;
- row names and column names are present and unique;
- `nrow(design)` equals the number of matrix columns;
- optional `group` values are aligned to the matrix columns.

For large analyses, create the BPCells matrix before calling `edger.bp`.
Keep a column metadata table alongside it with any sample, cluster,
batch, or donor fields needed to build design matrices.

## A Small Runnable Example

The example below simulates a small two-group pseudobulk experiment and
writes the counts to a temporary BPCells directory.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`edger.bp`](https://marusakod.github.io/edger.bp/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`BPCells`](https://bnprks.github.io/BPCells)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`Matrix`](https://codeberg.org/jaganmn/Matrix)`)`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``7``)`\
\
`ngenes`` ``<-`` ``120L`\
`nsamples`` ``<-`` ``8L`\
`group`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"ctrl"``, ``"stim"``)``, each ``=`` ``4L``)``)`\
\
`base_mean`` ``<-`` `[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(``ngenes``, shape ``=`` ``2``, rate ``=`` ``0.15``)`\
`fold_change`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``1``, ``ngenes``)`\
`fold_change``[`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``12L``)``]`` ``<-`` ``2.5`\
\
`counts`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``nsamples``)``, ``function``(``j``)`` ``{`\
`    ``mu`` ``<-`` ``base_mean`` ``*`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``group``[[``j``]``]`` ``==`` ``"stim"``, ``fold_change``, ``1``)`\
`    `[`as.integer`](https://rdrr.io/r/base/integer.html)`(`[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(``ngenes``, mu ``=`` ``mu``, size ``=`` ``8``)``)`\
`}``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``ngenes``)``)`\
\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"gene"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``ngenes``)``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"sample"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``nsamples``)``)`\
\
`bp_dir`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``"edger-bp-counts-"``)`\
`sparse_counts`` ``<-`` `[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(``counts``, sparse ``=`` ``TRUE``)`\
`sparse_counts`` ``<-`` ``BPCells``::`[`convert_matrix_type`](https://bnprks.github.io/BPCells/reference/convert_matrix_type.html)`(``sparse_counts``, ``"uint32_t"``)`\
`bp_counts`` ``<-`` ``BPCells``::`[`write_matrix_dir`](https://bnprks.github.io/BPCells/reference/matrix_io.html)`(``sparse_counts``, dir ``=`` ``bp_dir``)`\
\
`design`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``group``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``design``)`\
`#> [1] "(Intercept)" "groupstim"`

The audited package versions can be checked explicitly:

\
[`bp_check_versions`](https://marusakod.github.io/edger.bp/reference/bp_check_versions.md)`(``error ``=`` ``FALSE``)`\
`#>   package installed expected version_supported   ok     stack`\
`#> 1       R     4.5.3    4.5.3              TRUE TRUE conda-r45`\
`#> 2   edgeR     4.8.2    4.8.2              TRUE TRUE conda-r45`\
`#> 3   limma    3.66.0   3.66.0              TRUE TRUE conda-r45`\
`#> 4 BPCells     0.3.1    0.3.1              TRUE TRUE conda-r45`

## Likelihood-Ratio Test

[`bp_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)
runs preparation, optional filtering, TMM normalization, dispersion
estimation, and blockwise GLM LRT testing. The returned object keeps
those stages separate.

\
`lrt_fit`` ``<-`` `[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"groupstim"``,`\
`    dispersion.method ``=`` ``"stream"``,`\
`    block.size ``=`` ``50L``,`\
`    workers ``=`` ``1L`\
`)`\
\
[`names`](https://rdrr.io/r/base/names.html)`(``lrt_fit``)`\
`#> [1] "prepare"    "dispersion" "lrt"`\
[`head`](https://rdrr.io/r/utils/head.html)`(``lrt_fit``$``lrt``$``table``[`[`order`](https://rdrr.io/r/base/order.html)`(``lrt_fit``$``lrt``$``table``$``PValue``)``, ``]``)`\
`#>             logFC   logCPM        LR      PValue`\
`#> gene50 -1.3440732 13.14280 10.087414 0.001492854`\
`#> gene75  1.2826286 13.48391  8.719293 0.003148602`\
`#> gene86  0.9841068 13.39220  5.240966 0.022061027`\
`#> gene49 -0.9028414 13.35884  4.728981 0.029658459`\
`#> gene29  0.8451257 13.37477  3.901494 0.048243195`\
`#> gene45  0.7526657 13.49282  3.164025 0.075277385`

Important components are:

- `prepare`: filtered counts, library sizes, normalization factors, and
  offsets;
- `dispersion`: common, trended, and tagwise dispersion estimates;
- `lrt$table`: the edgeR-style result table for the requested
  coefficient or contrast.

## QL F-Test and Treat

[`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)
and
[`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)
use edgeR’s modern quasi-likelihood path by default. They can also
accept a supplied scalar or gene-level dispersion vector when that
matches the intended edgeR analysis.

\
`qlf_fit`` ``<-`` `[`bp_qlf`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"groupstim"``,`\
`    dispersion.method ``=`` ``"stream"``,`\
`    block.size ``=`` ``50L``,`\
`    workers ``=`` ``1L`\
`)`\
\
[`head`](https://rdrr.io/r/utils/head.html)`(``qlf_fit``$``qlf``$``table``[`[`order`](https://rdrr.io/r/base/order.html)`(``qlf_fit``$``qlf``$``table``$``PValue``)``, ``]``)`\
`#>             logFC   logCPM        F      PValue`\
`#> gene50 -1.3439926 13.14280 9.692360 0.001932632`\
`#> gene75  1.2827450 13.48391 8.211283 0.004299308`\
`#> gene86  0.9840272 13.39220 4.971189 0.026118029`\
`#> gene49 -0.9028356 13.35884 4.499580 0.034285776`\
`#> gene29  0.8453132 13.37477 3.701896 0.054792438`\
`#> gene45  0.7526489 13.49282 2.979832 0.084789141`

[`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)
tests whether the effect exceeds a log2 fold-change threshold. The
default threshold follows edgeR’s
[`glmTreat()`](https://rdrr.io/pkg/edgeR/man/glmTreat.html) default,
`log2(1.2)`.

\
`treat_fit`` ``<-`` `[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"groupstim"``,`\
`    dispersion.method ``=`` ``"stream"``,`\
`    lfc ``=`` `[`log2`](https://rdrr.io/r/base/Log.html)`(``1.2``)``,`\
`    block.size ``=`` ``50L``,`\
`    workers ``=`` ``1L`\
`)`\
\
[`head`](https://rdrr.io/r/utils/head.html)`(``treat_fit``$``treat``$``table``[`[`order`](https://rdrr.io/r/base/order.html)`(``treat_fit``$``treat``$``table``$``PValue``)``, ``]``)`\
`#>             logFC unshrunk.logFC   logCPM      PValue`\
`#> gene50 -1.3439926     -1.3512817 13.14280 0.003282173`\
`#> gene75  1.2827450      1.2885989 13.48391 0.006672171`\
`#> gene86  0.9840272      0.9885210 13.39220 0.035175862`\
`#> gene49 -0.9028356     -0.9067826 13.35884 0.045901236`\
`#> gene29  0.8453132      0.8491799 13.37477 0.069344274`\
`#> gene39  0.8594424      0.8693028 12.20338 0.100362196`

For legacy edgeR quasi-likelihood comparisons, use `qlf.args` or
`treat.args`:

\
[`bp_qlf`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"groupstim"``,`\
`    qlf.args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``legacy ``=`` ``TRUE``)`\
`)`\
\
[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"groupstim"``,`\
`    treat.args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``legacy ``=`` ``TRUE``)`\
`)`

## One-vs-Rest Marker Tests

A common marker workflow builds one binary design per cluster. The
pseudobulk matrix still has genes as rows and pseudobulk columns as
samples; only the design changes.

\
`cluster`` ``<-`` ``col_data``$``cluster`\
`selected`` ``<-`` ``cluster`` ``==`` ``"cl01"`\
`design_marker`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``selected``)`\
\
`marker_cl01`` ``<-`` `[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design_marker``,`\
`    group ``=`` ``selected``,`\
`    coef ``=`` ``"selectedTRUE"``,`\
`    lfc ``=`` `[`log2`](https://rdrr.io/r/base/Log.html)`(``1.2``)``,`\
`    block.size ``=`` ``1000L``,`\
`    workers ``=`` ``4L`\
`)`\
\
`marker_table`` ``<-`` ``marker_cl01``$``treat``$``table`

This binary one-vs-rest design directly compares the selected cluster to
all other columns pooled together. A shared `~ 0 + cluster` design with
contrasts can be faster for many clusters, but it is a different
estimand and should be chosen deliberately.

## Storage Order and Parallelism

BPCells matrices have a physical storage order. Column-ordered storage
is useful for sample-wise operations such as library sizes and TMM
normalization. Row- ordered storage is useful for blockwise gene
operations. If you need a companion matrix with the other storage order,
use
[`bp_ensure_storage_order()`](https://marusakod.github.io/edger.bp/reference/bp_ensure_storage_order.md).

\
`bp_counts_row`` ``<-`` `[`bp_ensure_storage_order`](https://marusakod.github.io/edger.bp/reference/bp_ensure_storage_order.md)`(`\
`    ``bp_counts``,`\
`    order ``=`` ``"row"``,`\
`    outdir ``=`` ``"cache/counts_row"``,`\
`    tmpdir ``=`` ``"cache/tmp"`\
`)`

The main tuning arguments are:

- `block.size`: number of genes materialized per edgeR block;
- `workers`: number of row blocks processed in parallel where supported;
- Slurm resources: for full datasets, run package workflows in batch
  jobs and budget memory as roughly proportional to
  `workers * block.size * ncol(y)`.

Small interactive checks should use `workers = 1L` or `2L`. Full
pseudobulk marker runs should use Slurm and an explicit wall-time/QoS
choice appropriate to the cluster.

## Full-Scale Workflow

A typical full-scale analysis is:

1.  Build or load a genes-by-pseudobulk BPCells count matrix.
2.  Load column metadata and construct the design matrix.
3.  Run a small smoke test on a subset of genes or clusters.
4.  Run the full workflow with recorded `block.size`, `workers`, package
    versions, and Slurm resources.
5.  Save the result table plus enough metadata to reproduce the design
    and contrast.

For example:

\
`bp_counts`` ``<-`` ``BPCells``::`[`open_matrix_dir`](https://bnprks.github.io/BPCells/reference/matrix_io.html)`(``"cache/counts_col"``)`\
`col_data`` ``<-`` `[`read.delim`](https://rdrr.io/r/utils/read.table.html)`(``"cache/sample_metadata.tsv"``)`\
\
`group`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(``col_data``$``condition``)`\
`design`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``group``)`\
\
`fit`` ``<-`` `[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``"grouptreated"``,`\
`    block.size ``=`` ``1000L``,`\
`    workers ``=`` ``8L`\
`)`\
\
[`write.table`](https://rdrr.io/r/utils/write.table.html)`(`\
`    ``fit``$``treat``$``table``,`\
`    file ``=`` ``"results/edgeR_treat.tsv"``,`\
`    sep ``=`` ``"\t"``,`\
`    quote ``=`` ``FALSE``,`\
`    col.names ``=`` ``NA`\
`)`
