# Pipelines

### GWAS Columns:

With each GWAS file, you can specify column names ex. `{"P":"your_gwas_pval_col", ...}`, if you do not specify header names it will assume your GWAS has the headers below.

If you do not specify the expected format of the GWAS, all pipelines will assume the following headers:
The standard column naming for GWASes are:

|           | CHR | BP  | BETA | SE  | P   | EA  | OA  | EAF | SNP | RSID |
|-----------|-----|-----|------|-----|-----|-----|-----|-----|-----|:-----|
| Mandatory | x   | x   | x    | x   | x   | x   | x   | x   |     |      |

### GWAS Standardisation

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as `CHR:POS_EA_OA`, where EA and OA are ordered alphabetically, and the subsequent BETA and EAF will be adjusted accordingly

## Collider Bias Pipeline

### input.json

* `incident`: Incident GWAS summary statistics file and optional column name map
* `subsequent`: Subsequent GWAS summary statistics file and optional column name map
* `ancestry`: ancestry of incident and subsequent GWAS, both must be same ancestry.  Options: "EUR", "EAS", "AFR", "AMR", "SAS"
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## Ancestry Comparison

### input.json
* List of input GWAS files:
  * `file`: GWAS summary statistics file
  * `columns`: optional column name map
  * `ancestry`: of incident and subsequent GWAS.  Options: "EUR", "EAS", "AFR", "AMR", "SAS"
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

