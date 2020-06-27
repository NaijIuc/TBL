# TargetByLasso
Predicting key therapeutic targets of a disease with Lasso regression, based on the drug test results matrix.
---
## Installation
```R
# install.packages("devtools")
 devtools::install_github("NaijIuc/TargetByLasso")
```
## Getting started

``` r
library(TargetByLasso)

TargetPre(single.effects, combo.effects, ref.matrix)
```
For the main function `TargetPre`, Three files are requested, the single drug effects that can be binomial or continuous value (single.effects), the combination drug effects that can be binomial only to illustrate either synergism or antagonism (combo.effects), and drug-target matrix (ref.matrix) that can be obtained from [STICH](http://stitch.embl.de/cgi/download.pl?UserId=UfyynCSx9VZy&sessionId=TpAyudNTkNKq). STITCH provides compound-protein relationships with the format of PubChem CIDs vs Ensembl peptide IDs, and users can convert to the desired type accordingly.
