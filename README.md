# Connectivity Mapping (cmap) Application

A web server-based application can be found [here](https://jessbwhite.shinyapps.io/cmap/).

## Call application in R

```
library(shiny)

runGitHub("cmap", "jessicaw9910")
```

## Input file

To use this workflow, enter a CSV file containing 3 columns of data from differential gene expression (DGE) analysis to create a rank ordered list of a given biological condition.  An example file is provided in the `test.csv` file.
1) **gene symbol (column 1)** - be sure these gene names map to the latest using `AnnotationHub()`
2) **logFC change (column 2)** - log<sub>2</sub> fold change relative to reference condition from DGE
3) **ranking metric (column 3)** - gene-level statistics, such as p-value, adjusted p-value, or absolute t-statistic, from a differential expression test.  Note that identical values will force arbitrary ranking so select a metric with minimum number of duplicates.  This is frequently an issue if adjusted p-values have been corrected using Benjamini-Hochberg.

## data

+ **LINCS_L1000_Chem_Pert.Rda** - contains cleaned [Enrichr](https://maayanlab.cloud/Enrichr/#stats) LINCS L1000 up-regulated and down-regulated gene lists in the form of matrices and lists
+ **new_moa.csv** - contains mapping of drugs to mechanisms of action from [clue.io](https://clue.io/data/REP#REP), [L1000 FWD](https://maayanlab.cloud/L1000FWD/), chemical supplier website (e.g., Tocris, Sigma Aldrich), and literature