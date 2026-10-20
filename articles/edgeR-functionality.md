# edgeR functionality in edger.bp

`edger.bp` is a BPCells-backed implementation layer for selected edgeR
workflows. The package aims to reproduce edgeR results for supported
inputs while avoiding full dense materialization of large pseudobulk
matrices.

It does not replace the edgeR statistical documentation. Use edgeR’s
guidance to choose the model, contrast, test, and filtering strategy;
use `edger.bp` when the same workflow needs BPCells-backed storage and
blockwise execution.

## Function Mapping

`edger.bp` has two layers. The high-level workflow functions run
preparation, normalization, optional filtering, dispersion estimation
when needed, and the final GLM test. The lower-level GLM functions
assume those preparation steps have already been handled and run only
the blockwise fitting/testing stage.

| edgeR operation | High-level workflow | Lower-level function | Notes |
|----|----|----|----|
| [`filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html) | `bp_prepare_dge(filter = TRUE)` | [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md) | Runs on BPCells-backed counts. |
| `normLibSizes(method = "TMM")` | [`bp_prepare_dge()`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md) | [`bp_calc_norm_factors()`](https://marusakod.github.io/edger.bp/reference/bp_calc_norm_factors.md) | `TMM` and `none` are supported. |
| [`DGEList()`](https://rdrr.io/pkg/edgeR/man/DGEList.html) preparation | [`bp_prepare_dge()`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md) | none | Returns filtered counts, library sizes, normalization factors, and offsets. |
| [`estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html) | [`bp_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md), [`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md), or [`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md) when `dispersion = NULL` | [`bp_estimate_disp_stream()`](https://marusakod.github.io/edger.bp/reference/bp_estimate_disp_stream.md) | Main streamed estimator. |
| dense [`estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html) oracle | [`bp_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md), [`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md), or [`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md) with `dispersion.method = "dense"` | [`bp_estimate_disp_dense()`](https://marusakod.github.io/edger.bp/reference/bp_estimate_disp_dense.md) | Intended for small-data checks. |
| [`glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html) + [`glmLRT()`](https://rdrr.io/pkg/edgeR/man/glmLRT.html) | [`bp_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md) | [`bp_glm_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_glm_lrt.md) | Supports coefficient and contrast tests. |
| [`glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) + [`glmQLFTest()`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html) | [`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md) | [`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md) | Modern QL is the default; legacy QL is available. |
| [`glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) + [`glmTreat()`](https://rdrr.io/pkg/edgeR/man/glmTreat.html) | [`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md) | [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md) | Modern Treat is the default; legacy Treat is available. |

## Choosing a Statistical Test

The package provides three test families and seven test configurations.
The modern quasi-likelihood implementations are the defaults.

| Test configuration | Lower-level call | Null hypothesis | Reference distribution | Multiple coefficients? | Typical use |
|----|----|----|----|----|----|
| Likelihood-ratio test | [`bp_glm_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_glm_lrt.md) | Selected coefficient(s) or contrast(s) equal zero | Chi-squared | Yes | Fast testing when NB dispersions are treated as fixed. |
| Modern QLF | `bp_glm_qlf(..., legacy = FALSE)` | Selected coefficient(s) or contrast(s) equal zero | Modern quasi-likelihood F | Yes | Recommended general-purpose differential-expression test. |
| Legacy QLF | `bp_glm_qlf(..., legacy = TRUE)` | Selected coefficient(s) or contrast(s) equal zero | Legacy quasi-likelihood F | Yes | Compatibility with older edgeR QL workflows. |
| Modern Treat, interval null | `bp_glm_treat(..., legacy = FALSE, null = "interval")` | Absolute log2 fold change is no larger than `lfc` | Modern QL threshold test | One contrast | Recommended when a minimum biologically relevant effect is required. |
| Modern Treat, worst-case null | `bp_glm_treat(..., legacy = FALSE, null = "worst.case")` | Threshold null evaluated at its least-favourable boundary | Modern QL threshold test | One contrast | More conservative threshold formulation. |
| Legacy Treat, interval null | `bp_glm_treat(..., legacy = TRUE, null = "interval")` | Absolute log2 fold change is no larger than `lfc` | Legacy QL threshold test | One contrast | Threshold testing compatible with older edgeR QL workflows. |
| Legacy Treat, worst-case null | `bp_glm_treat(..., legacy = TRUE, null = "worst.case")` | Threshold null evaluated at its least-favourable boundary | Legacy QL threshold test | One contrast | Conservative threshold testing with legacy QL semantics. |

Each family has a lower-level fitting/testing function and an end-to-end
workflow:

| Test family | Lower-level function | End-to-end workflow |
|----|----|----|
| LRT | [`bp_glm_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_glm_lrt.md) | [`bp_lrt()`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md) |
| QLF | [`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md) | [`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md) |
| Treat | [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md) | [`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md) |

The lower-level functions require a design, dispersions, and offsets or
library sizes. The end-to-end workflows can also filter genes, calculate
TMM normalization factors, and estimate dispersions. Select legacy QL
semantics in the high-level workflows through `qlf.args` or
`treat.args`:

\
[`bp_qlf`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)`(`\
`    ``y``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``2``,`\
`    qlf.args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``legacy ``=`` ``TRUE``)`\
`)`\
\
[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``y``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``2``,`\
`    lfc ``=`` `[`log2`](https://rdrr.io/r/base/Log.html)`(``1.2``)``,`\
`    null ``=`` ``"interval"``,`\
`    treat.args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``legacy ``=`` ``TRUE``)`\
`)`

The default Treat threshold is `log2(1.2)`. Setting `lfc = 0` delegates
to the corresponding QLF test. Treat tests one coefficient or one
contrast; LRT and QLF also support joint tests of multiple coefficients
or contrast columns.

As a practical starting point, use modern QLF for a zero-effect primary
analysis, modern Treat when the scientific question includes a minimum
effect size, and LRT when speed or exploratory sensitivity is the
priority. Legacy QLF and Treat are primarily compatibility and
sensitivity-analysis paths.

In most analyses, start with a high-level workflow:

\
[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(``y``, design ``=`` ``design``, group ``=`` ``group``, coef ``=`` ``"groupB"``)`\
[`bp_qlf`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)`(``y``, design ``=`` ``design``, group ``=`` ``group``, coef ``=`` ``"groupB"``)`\
[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(``y``, design ``=`` ``design``, group ``=`` ``group``, coef ``=`` ``"groupB"``)`

These calls correspond to the usual native edgeR sequence of
constructing a `DGEList`, filtering, normalizing, estimating
dispersions, fitting a GLM, and testing a coefficient or contrast.

Use the lower-level functions when you want explicit control over
intermediate objects. Common reasons are to reuse one prepared matrix,
supply a dispersion vector estimated elsewhere, compare one stage
directly with native edgeR, or avoid repeating preparation inside a
loop.

For example, the high-level call:

\
[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(``y``, design ``=`` ``design``, group ``=`` ``group``, coef ``=`` ``"groupB"``)`

is roughly equivalent to:

\
`prep`` ``<-`` `[`bp_prepare_dge`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md)`(``y``, design ``=`` ``design``, group ``=`` ``group``)`\
`disp`` ``<-`` `[`bp_estimate_disp_stream`](https://marusakod.github.io/edger.bp/reference/bp_estimate_disp_stream.md)`(`\
`    ``prep``$``counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    lib.size ``=`` ``prep``$``effective.lib.size``,`\
`    offset ``=`` ``prep``$``offset`\
`)`\
\
`lrt`` ``<-`` `[`bp_glm_lrt`](https://marusakod.github.io/edger.bp/reference/bp_glm_lrt.md)`(`\
`    ``prep``$``counts``,`\
`    design ``=`` ``design``,`\
`    dispersion ``=`` ``disp``$``tagwise.dispersion``,`\
`    offset ``=`` ``prep``$``offset``,`\
`    coef ``=`` ``"groupB"`\
`)`

The same relationship holds for
[`bp_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)
versus
[`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md),
and
[`bp_treat()`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)
versus
[`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md).

## Supported Scope

The audited path covers non-negative count matrices with genes as rows
and samples or pseudobulk columns as columns. The current package
focuses on:

- TMM normalization or no normalization;
- edgeR-compatible expression filtering;
- common, trended, and tagwise dispersion estimation;
- coefficient and contrast tests;
- LRT, modern QLF, legacy QLF, modern Treat, and legacy Treat;
- scalar or gene-level supplied dispersions;
- sample-level offsets represented as a vector, or row-constant offset
  matrices.

The following edgeR features are not part of the supported path unless
they are separately audited and tested:

- robust dispersion or robust QL estimation;
- observation weights;
- gene-specific offset matrices;
- normalization methods other than `TMM` and `none`;
- exact compatibility with edgeR versions outside the audited dependency
  stacks.

Check the audited versions with:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`edger.bp`](https://marusakod.github.io/edger.bp/)`)`\
\
[`bp_check_versions`](https://marusakod.github.io/edger.bp/reference/bp_check_versions.md)`(``error ``=`` ``FALSE``)`\
`#>   package installed expected version_supported   ok     stack`\
`#> 1       R     4.5.3    4.5.3              TRUE TRUE conda-r45`\
`#> 2   edgeR     4.8.2    4.8.2              TRUE TRUE conda-r45`\
`#> 3   limma    3.66.0   3.66.0              TRUE TRUE conda-r45`\
`#> 4 BPCells     0.3.1    0.3.1              TRUE TRUE conda-r45`

## Native edgeR and edger.bp Side by Side

This small example runs the same two-group workflow with native edgeR
and `edger.bp`.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`edgeR`](https://bioinf.wehi.edu.au/edgeR/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`BPCells`](https://bnprks.github.io/BPCells)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`Matrix`](https://codeberg.org/jaganmn/Matrix)`)`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``48``)`\
\
`ngenes`` ``<-`` ``115L`\
`nsamples`` ``<-`` ``6L`\
`group`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"A"``, ``"B"``)``, each ``=`` ``3L``)``)`\
`design`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``group``)`\
\
`base_mean`` ``<-`` `[`rgamma`](https://rdrr.io/r/stats/GammaDist.html)`(``ngenes``, shape ``=`` ``2``, rate ``=`` ``0.2``)`\
`counts`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(``nsamples``)``, ``function``(``j``)`` ``{`\
`    `[`as.integer`](https://rdrr.io/r/base/integer.html)`(`[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(``ngenes``, mu ``=`` ``base_mean``, size ``=`` ``4``)``)`\
`}``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``ngenes``)``)`\
\
`counts``[`[`sample`](https://rdrr.io/r/base/sample.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``counts``)``, size ``=`` `[`floor`](https://rdrr.io/r/base/Round.html)`(`[`length`](https://rdrr.io/r/base/length.html)`(``counts``)`` ``*`` ``0.35``)``)``]`` ``<-`` ``0L`\
`counts``[``1``:``3``, ``]`` ``<-`` ``0L`\
\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"gene"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``ngenes``)``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``counts``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"sample"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(``nsamples``)``)`\
\
`bp_dir`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``"edger-bp-functionality-"``)`\
`bp_counts`` ``<-`` ``BPCells``::`[`write_matrix_dir`](https://bnprks.github.io/BPCells/reference/matrix_io.html)`(`\
`    ``BPCells``::`[`convert_matrix_type`](https://bnprks.github.io/BPCells/reference/convert_matrix_type.html)`(`[`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html)`(``counts``, sparse ``=`` ``TRUE``)``, ``"uint32_t"``)``,`\
`    dir ``=`` ``bp_dir`\
`)`

The native edgeR workflow is:

\
`dge`` ``<-`` `[`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html)`(``counts``)`\
`keep`` ``<-`` `[`filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)`(``dge``, design ``=`` ``design``, group ``=`` ``group``)`\
`dge`` ``<-`` ``dge``[``keep``, , keep.lib.sizes ``=`` ``FALSE``]`\
`dge`` ``<-`` `[`normLibSizes`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)`(``dge``)`\
`dge`` ``<-`` `[`estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)`(``dge``, design ``=`` ``design``)`\
\
`fit_lrt`` ``<-`` `[`glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html)`(``dge``, design ``=`` ``design``)`\
`edge_lrt`` ``<-`` `[`glmLRT`](https://rdrr.io/pkg/edgeR/man/glmLRT.html)`(``fit_lrt``, coef ``=`` ``2``)`\
\
`fit_ql`` ``<-`` `[`glmQLFit`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html)`(`\
`    ``dge``,`\
`    design ``=`` ``design``,`\
`    dispersion ``=`` ``dge``$``tagwise.dispersion``,`\
`    legacy ``=`` ``FALSE`\
`)`\
`edge_qlf`` ``<-`` `[`glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html)`(``fit_ql``, coef ``=`` ``2``)`\
`edge_treat`` ``<-`` `[`glmTreat`](https://rdrr.io/pkg/edgeR/man/glmTreat.html)`(``fit_ql``, coef ``=`` ``2``)`

The corresponding high-level `edger.bp` calls are:

\
`bp_lrt_fit`` ``<-`` `[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``2``,`\
`    block.size ``=`` ``10L``,`\
`    workers ``=`` ``1L`\
`)`\
\
`bp_qlf_fit`` ``<-`` `[`bp_qlf`](https://marusakod.github.io/edger.bp/reference/bp_qlf.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``2``,`\
`    block.size ``=`` ``10L``,`\
`    workers ``=`` ``1L`\
`)`\
\
`bp_treat_fit`` ``<-`` `[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(`\
`    ``bp_counts``,`\
`    design ``=`` ``design``,`\
`    group ``=`` ``group``,`\
`    coef ``=`` ``2``,`\
`    block.size ``=`` ``10L``,`\
`    workers ``=`` ``1L`\
`)`

A compact numerical comparison:

\
`max_abs_diff`` ``<-`` ``function``(``x``, ``y``)`` ``{`\
`    `[`max`](https://rdrr.io/r/base/Extremes.html)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(``x`` ``-`` ``y``)``, na.rm ``=`` ``TRUE``)`\
`}`\
\
`comparison`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`    quantity ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`        ``"filter keep vector"``,`\
`        ``"TMM norm factors"``,`\
`        ``"tagwise dispersion"``,`\
`        ``"LRT p-values"``,`\
`        ``"modern QLF p-values"``,`\
`        ``"modern Treat p-values"`\
`    ``)``,`\
`    max_abs_diff ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`        ``max_abs_diff``(`[`as.integer`](https://rdrr.io/r/base/integer.html)`(``bp_lrt_fit``$``prepare``$``keep``)``, `[`as.integer`](https://rdrr.io/r/base/integer.html)`(``keep``)``)``,`\
`        ``max_abs_diff``(``bp_lrt_fit``$``prepare``$``norm.factors``, ``dge``$``samples``$``norm.factors``)``,`\
`        ``max_abs_diff``(``bp_lrt_fit``$``dispersion``$``tagwise.dispersion``, ``dge``$``tagwise.dispersion``)``,`\
`        ``max_abs_diff``(``bp_lrt_fit``$``lrt``$``table``$``PValue``, ``edge_lrt``$``table``$``PValue``)``,`\
`        ``max_abs_diff``(``bp_qlf_fit``$``qlf``$``table``$``PValue``, ``edge_qlf``$``table``$``PValue``)``,`\
`        ``max_abs_diff``(``bp_treat_fit``$``treat``$``table``$``PValue``, ``edge_treat``$``table``$``PValue``)`\
`    ``)`\
`)`\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(`\
`    `[`transform`](https://rdrr.io/r/base/transform.html)`(``comparison``, max_abs_diff ``=`` `[`formatC`](https://rdrr.io/r/base/formatc.html)`(``max_abs_diff``, format ``=`` ``"e"``, digits ``=`` ``2``)``)`\
`)`

| quantity              | max_abs_diff |
|:----------------------|:-------------|
| filter keep vector    | 0.00e+00     |
| TMM norm factors      | 0.00e+00     |
| tagwise dispersion    | 0.00e+00     |
| LRT p-values          | 0.00e+00     |
| modern QLF p-values   | 1.52e-13     |
| modern Treat p-values | 1.16e-13     |

These comparisons are the package’s main correctness target: for
supported inputs, the BPCells-backed workflow should match the
corresponding native edgeR workflow to numerical tolerance.

## Coefficients, Contrasts, and Marker Designs

`coef` and `contrast` follow edgeR conventions. Coefficients can be
selected by index or by design column name:

\
[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(``bp_counts``, design ``=`` ``design``, group ``=`` ``group``, coef ``=`` ``"groupB"``)`\
\
`contrast`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"(Intercept)"`` ``=`` ``0``, ``"groupB"`` ``=`` ``1``)`\
[`bp_lrt`](https://marusakod.github.io/edger.bp/reference/bp_lrt.md)`(``bp_counts``, design ``=`` ``design``, group ``=`` ``group``, contrast ``=`` ``contrast``)`

For marker testing, a binary one-vs-rest design and a multi-cluster
contrast design answer related but different questions:

\
`# Binary one-vs-rest: selected cluster against all other columns pooled together.`\
`selected`` ``<-`` ``col_data``$``cluster`` ``==`` ``"cl01"`\
`design_binary`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``selected``)`\
[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(``bp_counts``, design ``=`` ``design_binary``, group ``=`` ``selected``, coef ``=`` ``"selectedTRUE"``)`\
\
`# Multi-cluster contrasts: one shared model with explicit cluster contrasts.`\
`cluster`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(``col_data``$``cluster``)`\
`design_cluster`` ``<-`` `[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)`(``~`` ``0`` ``+`` ``cluster``)`\
`contrast_cl01`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``-``1`` ``/`` ``(`[`nlevels`](https://rdrr.io/r/base/nlevels.html)`(``cluster``)`` ``-`` ``1``)``, `[`nlevels`](https://rdrr.io/r/base/nlevels.html)`(``cluster``)``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``contrast_cl01``)`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``design_cluster``)`\
`contrast_cl01``[``"clustercl01"``]`` ``<-`` ``1`\
[`bp_treat`](https://marusakod.github.io/edger.bp/reference/bp_treat.md)`(``bp_counts``, design ``=`` ``design_cluster``, group ``=`` ``cluster``, contrast ``=`` ``contrast_cl01``)`

The second pattern can be faster for many clusters because the QL fit
can be shared, but it changes the estimand. Choose the design for the
scientific question first, then tune the implementation.

## Runtime Controls

`edger.bp` materializes blocks of genes for edgeR GLM fitting. The main
controls are:

- `block.size`: larger blocks reduce overhead but increase memory per
  worker;
- `workers`: parallel row-block workers for supported stages;
- BPCells storage order: column-ordered matrices are efficient for
  sample-wise operations, while row-ordered companions can help
  gene-block operations.

For full pseudobulk datasets, use Slurm rather than an interactive
session and record the package versions, design matrix, contrast,
`block.size`, `workers`, and peak memory.
