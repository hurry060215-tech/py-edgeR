# pyedger

A **pure-Python port of [Bioconductor edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)** (Robinson, McCarthy & Smyth, *Bioinformatics* 2010) — negative-binomial models for differential expression of count data.

- **No `rpy2`**, no R install — the edgeR negative-binomial GLM workflow reimplemented in NumPy / SciPy
- The canonical pipelines: `DGEList → calcNormFactors → estimateDisp → glmQLFit → glmQLFTest` and the classic `exactTest`
- TMM normalization, common / trended / tagwise dispersion, GLM and quasi-likelihood F-tests, `filterByExpr`, `cpm` / `aveLogCPM`
- Both Python-style (`glm_fit`, `estimate_disp`, `top_tags`) and R-style (`glmFit`, `estimateDisp`, `topTags`) names exported

> This is a **standalone mirror** of the implementation developed in [`omicverse`](https://github.com/Starlitnightly/omicverse), where it powers the edgeR differential-expression backend of `ov.bulk` / `pyDEG`.

## Install

```bash
pip install pyedger
```

## Quick start

```python
import numpy as np
import pyedger

# counts: genes x samples raw count matrix; group: per-sample condition labels
dge = pyedger.DGEList(counts=counts, group=group)
dge = pyedger.calcNormFactors(dge)                       # TMM normalization
keep = pyedger.filterByExpr(dge, group=group)
dge = dge[keep]

# Quasi-likelihood F-test workflow (the recommended edgeR pipeline)
dge = pyedger.estimateDisp(dge, design)
fit = pyedger.glmQLFit(dge, design)
qlf = pyedger.glmQLFTest(fit, coef=1)
res = pyedger.topTags(qlf, n=np.inf)
res.head()
```

### Classic exact test

```python
dge = pyedger.estimateDisp(dge, design)
et = pyedger.exactTest(dge, pair=("control", "treated"))
pyedger.topTags(et)
```

## API

| Python | R counterpart |
|---|---|
| `DGEList` | `DGEList` |
| `calc_norm_factors` / `calcNormFactors` | `calcNormFactors` (TMM) |
| `filter_by_expr` / `filterByExpr` | `filterByExpr` |
| `estimate_disp` / `estimateDisp` | `estimateDisp` |
| `glm_fit` / `glmFit`, `glm_lrt` / `glmLRT` | `glmFit`, `glmLRT` |
| `glm_ql_fit` / `glmQLFit`, `glm_qlf_test` / `glmQLFTest` | `glmQLFit`, `glmQLFTest` |
| `exact_test` / `exactTest` | `exactTest` |
| `glm_treat` / `glmTreat` | `glmTreat` |
| `cpm`, `ave_log_cpm` / `aveLogCPM` | `cpm`, `aveLogCPM` |
| `top_tags` / `topTags` | `topTags` |
| `decide_tests_dge` / `decideTests` | `decideTestsDGE` |
| `DGEGLM`, `DGELRT`, `DGEExact`, `TestResults` | the corresponding S4 classes |

## Citation

> Robinson, M.D., McCarthy, D.J., Smyth, G.K. **edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.** *Bioinformatics* 26(1), 139–140 (2010).

…and acknowledge omicverse / this repo for the Python port.

## License

LGPL-3.0-or-later — matches the upstream Bioconductor package.
