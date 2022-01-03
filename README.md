
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEARS

<!-- badges: start -->

<!-- badges: end -->

The goal of GEARS is to reconstruct Gaussian Bayesian network and
compare two networks (identical or differential) with graph ordering
unknown.

## Installation

You can install the released version of GEARS from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GEARS")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("han16/GEARS")
```

In case people have installation problem, try this

    remove.packages(c("curl","httr"))
    install.packages(c("curl", "httr"))
    Sys.setenv(CURL_CA_BUNDLE="/usr/lib64/microsoft-r/3.4/lib64/R/lib/microsoft-r-cacert.pem")
    devtools::install_git("https://github.com/han16/GEARS")

## Example

``` r
library(GEARS)
## basic example code
head(exampleData1)
#>    V1  V2  V3  V4  V5 V6  V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19
#> 1 1.5 1.5 0.0 0.0 0.0  0 1.5  2  2 2.5 1.5   0 0.0   2   0 1.5 2.0 0.0   2
#> 2 1.5 1.5 2.0 1.5 2.0  2 1.5  2  2 2.5 0.0   0 0.0   0   0 1.5 0.0 2.0   2
#> 3 1.5 1.5 0.0 1.5 2.0  2 1.5  2  2 2.5 1.5   2 0.0   0   2 0.0 1.5 2.0   2
#> 4 1.5 1.5 0.0 0.0 1.5  2 0.0  0  0 1.5 0.0   0 1.5   2   2 0.0 1.5 0.0   0
#> 5 1.5 0.0 1.5 0.0 1.5  0 1.5  0  2 2.0 0.0   0 1.5   2   2 1.5 0.0 2.0   0
#> 6 1.5 1.5 2.0 1.5 2.0  2 1.5  2  2 2.5 0.0   0 1.5   0   0 0.0 0.0 1.5   0
#>   V20 V21 V22 V23 V24 V25 V26 V27 V28 V29 V30 V31 V32 V33 V34 V35 V36 V37
#> 1   0 0.0 0.0 1.5   2   2 0.0   0   0 1.5 0.0 0.0 2.0   0 2.0 0.0   0   0
#> 2   0 2.5 0.0 0.0   0   0 0.0   0   0 0.0 0.0 0.0 1.5   0 2.0 0.0   0   0
#> 3   0 2.5 0.0 1.5   0   0 0.0   2   0 0.0 0.0 2.0 0.0   0 0.0 2.5   0   0
#> 4   0 2.0 1.5 2.0   0   2 2.5   0   0 0.0 1.5 0.0 2.0   0 0.0 0.0   0   0
#> 5   2 2.5 1.5 0.0   0   2 0.0   0   2 0.0 0.0 1.5 2.0   0 2.0 2.5   0   0
#> 6   2 2.0 0.0 0.0   0   0 0.0   0   0 1.5 0.0 2.0 0.0   2 2.5 0.0   0   0
#>   V38 V39 V40 V41 V42 V43 V44 V45
#> 1   0 0.0 1.5   0 0.0   2   0 2.0
#> 2   0 0.0 1.5   2 2.0   0   0 2.5
#> 3   0 0.0 0.0   0 0.0   0   0 0.0
#> 4   0 1.5 2.0   2 2.5   0   0 0.0
#> 5   0 0.0 0.0   0 0.0   0   0 0.0
#> 6   0 0.0 1.5   0 2.0   0   0 0.0
head(exampleData4)
#>            V1         V2         V3         V4         V5         V6
#> 1  0.80426840  1.9006522  1.0777911 -0.1099156  6.7435888  1.0748047
#> 2  0.12403565  0.5072030 -0.6071371  0.9988323  1.4933037  1.3914729
#> 3 -0.09277538  1.6601093 -1.7514093  1.6036877  1.1269293  3.2696831
#> 4  0.24335230 -1.1334361 -0.2062619  0.6096576  0.8399068 -0.3126217
#> 5  0.02536413 -0.5756416 -1.3068243  0.1720098 -2.7572069  1.3700969
#> 6 -0.89784935 -1.7635564 -1.8967425 -0.8790398 -9.5008889 -2.9371123
#>          V7        V8         V9        V10
#> 1  4.144335  3.354748  2.8457608  13.357390
#> 2  1.487572  2.005783  6.7149070  18.386733
#> 3  6.053086  7.161020 10.9266216  37.144373
#> 4 -2.513313 -1.144628 -0.3666188  -5.634854
#> 5 -1.492838 -3.789733  2.9199732   3.416417
#> 6 -6.030453 -7.717611 -9.3254433 -29.471658
```

  - `exampleData1` is the regression coefficient in the network
  - `exampleData4` is one simulated data set with number of node as the
    number of the columns, 10.

## Reference

[Gaussian Bayesian network comparisons with graph ordering
unknown](https://www.sciencedirect.com/science/article/abs/pii/S0167947320302474)
