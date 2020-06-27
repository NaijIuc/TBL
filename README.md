# TargetByLasso
---
Predicting key therapeutic targets of a disease with Lasso regression, based on the drug test results matrix.
---
## Installation
```R
# install.packages("devtools")
 devtools::install_github("tidyverse/tidyr")
```
## Getting started

``` r
library(TargetByLasso)
```
For the main function `TargetPre`, Three files are requested, the single drug effects that can be binomial or continuouse value, the combination drug effects that are binomial only, illustrating either synegism or antagism, and drug-target matrix that can be obtained from [STICH](http://stitch.embl.de/cgi/download.pl?UserId=UfyynCSx9VZy&sessionId=TpAyudNTkNKq). Among the three files, the drug-target matrix can be generated with funtion `Unknow`.
