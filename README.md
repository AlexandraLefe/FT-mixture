# FT-mixture: Frequency trajectory mixture

FT-mixture is an R tool for clustering mutation frequency trajectories
from time series of sequencing data. More broadly, it is designed for
clustering frequency trajectories from time series of count data. It is
suited for pooled sequencing data characterizing, for example, viral
genetic material in wastewater samples although it can be applied to an
aggregation of clinical samples. FT-mixture takes, as input, time-series
of mutation count data and read depth at related positions. It is
therefore meant to be used after the process of variant calling. It
returns an estimated number of groups of mutations of similar frequency
trajectory along with their proportion, frequency at time origin and
selection coefficient estimates. We assume that one group, named the
neutral-group, is not under selection (constant frequency trajectory)
and *K* ≥ 0 groups, named non-neutral group(s), are under positive or
negative selection.

More detailes are provided in the working article \[1\] entitled
*Unsupervised detection and fitness estimation of emerging SARS-CoV-2
variants: Application to wastewater samples (ANRS0160)* (Lefebvre et
al., 2025). This page contains the R code along with brief explanations
over a limited subset of analyses presented in the article.

# Required packages

The package *VGAM* is required for running the main algorithm of
FT-mixture. The package *readr* is required for dataset preparation from
outputs of a variant caller, in text data format files, to count
matrices.

``` r
install.packages("VGAM")
install.packages("readr")
```

Packages *gtools* and *pROC* are needed in simulation studies,
respectively, for handling label switching and for computing Area Under
the ROC Curve.

``` r
install.packages("gtools") # function permutations is used for handling label switching in simulation studies
install.packages("pROC") # for computing Area Under the ROC Curve in simulated studies
```

# Example of use over a simulated dataset

In this section, we propose an illustration of the main functions on a
simulated dataset. The R code is provided in the folder named *R_code*.

The file *generate_simulated_data.R* can be used for simulating a
dataset either with the same model as FT-mixture (with function
*generate_simu*) or the hidden random walk model (with function
*simu_hidden_RW*) presented in \[1\].

## Dataset simulation

The following code gives an example of use of the function
*simu.dataset* for simulating a dataset composed of a neutral group and
*K* = 2 non-neutral groups.

``` r
set.seed(420)
source("R_code/generate_simulated_data.R")

n = 200 # number of mutations
time = c(0, 5, 12, 20) # vector of differences (usually in days) between sampling date and date of origin
date.origin = NULL # sampling date of origin
param = list( # parameter of the model as a list of:
  pi = c(0.6, 0.3, 0.1), # group proportion where the first slot is associated to the neutral group 
  mu = c(0.5,-3), # logit of mutation frequency at the date of origin for non-neutral groups
  s = c(-1,2) / max(time), # selection coefficient of non-neutral groups 
  alpha = 10, # alpha parameter of the beta-binomial distribution of the frequency of the neutral group
  beta = 50 # beta parameter of the beta-binomial distribution of the frequency of the neutral group
)
lambda = 40 # mean of the Poisson parameter for sampling read depths
simu.dataset = generate_simu(n = n, 
                             time = time,
                             date.origin = date.origin,
                             param = param,
                             lambda = lambda)
```

This function returns a list composed of

-   `data`: a dataset in the required format, as a list composed of

    -   a *n* × (*m*+1) matrix *x* of mutation count where *n* is the
        number of mutations and (*m*+1) is the number of samples (time
        points). Mutation names (respectively sampling dates) are given
        in the header column (respectively the header row) as character
        vectors. If `date.origin` is NULL, the default date of origin
        used is 2019-01-01.

    -   a *n* × (*m*+1) matrix *d* of read depth at related positions.

-   `z`: a vector of size *n*. *z* ∈ {0, …, *K*}<sup>*n*</sup> where
    value 0 stands for the neutral group.

``` r
head(simu.dataset$data$x)
```

    ##       2019-01-01 2019-01-06 2019-01-13 2019-01-21
    ## mut.1         35         25         17         21
    ## mut.2          3          9          6         16
    ## mut.3          6          5          3          6
    ## mut.4         11          9          8          6
    ## mut.5         20         23         26         17
    ## mut.6         34         26         21         12

``` r
head(simu.dataset$data$d)
```

    ##       2019-01-01 2019-01-06 2019-01-13 2019-01-21
    ## mut.1         47         39         43         60
    ## mut.2         29         49         41         44
    ## mut.3         34         45         40         58
    ## mut.4         54         53         36         33
    ## mut.5         36         45         46         30
    ## mut.6         52         43         40         33

``` r
simu.dataset$z[1:min(n,30)]
```

    ##  [1] 1 2 0 0 1 1 0 1 1 0 0 1 0 0 0 1 0 0 0 0 1 0 0 0 2 1 1 0 0 0

## Simulated data analysis

### The main function of FT-mixture

The file *EM_algorithm.R* contains the main function of FT-mixture
called *ft_clust*. This function combines the function *em_ft_clust*
which is the implementation of the EM algorithm over the model and
*em.initilization* for computing data-driven initial values for
*em_ft_clust* (see supporting information S2_Text in \[1\] for a
detailed explanation and a flowchart of the process).

``` r
source("R_code/EM_algorithm.R")
res <- ft_clust (data = simu.dataset$data, # dataset (mandatory)
                 K = 2, # choose a number of non-neutral groups (mandatory)
                 n.pre.init = c(5,5), 
                 n.pre.it = c(5,5),
                 n.init = c(10,10), 
                 n.it = c(5,5), 
                 check.initialization = TRUE, 
                 tol = 1e-3,
                 niter = 2000,
                 verbose = FALSE)
```

The function *ft_clust* performs a clustering of mutations composing the
header column of input `data` conditional on the dataset `data` and the
number of non-neutral groups given as input `K`. It takes as input:

-   `data` (mandatory): a dataset in the required format as a list of
    -   a *n* × (*m*+1) matrix *x* of mutation count whith mutation
        names as header column and sampling dates as header row.
    -   a *n* × (*m*+1) matrix *d* of read depth at related positions.
-   `K` (mandatory): a number of non-neutral groups (*K* ≥ 0).
-   `n.pre.init` and `n.pre.it`: respectively the number of random
    initial values and the number of iterations of *em_initialization*
    for computing data-driven initial values for *em_ft_clust* (see
    supporting information S2_Text in \[1\] for a detailed explanation
    and a flowchart of the process) (default: (5,5) and (5,5)
    respectively).
-   `n.init` and `n.it`: respectively the number of data-driven initial
    values and the number of iterations of *em_ft_clust* during the
    initialization step (see supporting information S2_Text in \[1\] for
    a detailed explanation and a flowchart of the process) (default:
    (10,10) and (5,5) respectively). It is recommended to rise `n.init`
    to (40,40) and/or `n.it` to (30,30) in case of unstable results for
    a tested number *K* ≥ 3 non-neutral groups.
-   `check.initialization`: if TRUE, the function returns the vector of
    log-likelihoods reached after the `n.it` iterations of *em_ft_clust*
    using the `n.init` inital values (default: TRUE).
-   `tol`: the tolerance on parameter differences for the algorithm to
    converge. The absolute value of the difference between parameters
    computed at the current and the previous iteration is used and
    divided by the absolute value of the parameter at current iteration
    for comparing quantities of no unit (default : 10<sup>−3</sup>).
-   `niter`: the maximal number of iterations of *em_ft_clust* in case
    of no convergence (default: 2000).  
-   `verbose`: if TRUE, the algorithm sequentially prints, the
    intermediate values of log-likelihood and parameter estimates
    obtained during the initialization step (default: FALSE).

The function *ft_clust* returns as output:

-   `loglik`: the log-likelihood of the model
-   `pi`, `mu`, `s`, `alpha` and `beta`: parameter estimates
-   `eta`: posterior group assignment as a *n* × (*K*+1) output matrix
    with mutation names in the header column and group number in the
    header row. The first column is associated to the neutral group.
-   `size.smallest.group`: the size of the smallest group in terms of
    number of mutations.

### A selection of results over a simulated dataset

The function *handle_label_switching* in the file
*generate_simulated_data.R* (solely usefull for simulation studies) uses
the mean square error between the parameter used for simulations and
estimated parameters `mu` and `s`, where `s` is multiplies by the
maximal time difference in order to sum quantities of similar range. It
takes as input:

-   `res` an output of ft_clust over a simulated dataset
-   `param` the parameter used for simulating the dataset in the same
    format as input `param` of function *generate_simu*.

``` r
sorted.res = handle_label_switching(res = res, param = param) 
```

We report in this paragraph, the table of computed parameter estimates
conditional on the simulated dataset *simu.dataset$data* (see first
chunck of this document) and conditional on K=2 non-neutral groups. We
also report the Area Under the ROC Curve of posterior group assignments
where output `z` of *generate_simu* is used as response and output `eta`
of *ft_clust* is used as predictor.

``` r
K = length(sorted.res$res$mu) # number of non-neutral groups # 2

# table of parameter estimates
knitr::kable(matrix(data = round(c(unlist(sorted.res$param), unlist(sorted.res$res[c("pi", "mu", "s", "alpha","beta")])), 3), nrow = 2, byrow = TRUE, 
                    dimnames = list(c("``True'' parameter", "Parameter estimates"), c(paste0("pi", 0:K), paste0("mu", 1:K), paste0("s", 1:K), "alpha", "beta"))),
             caption = "Table of the parameter used for simulations (True parameter) and parameter estimates.")
```

|                      |   pi0 |  pi1 |   pi2 |   mu1 |    mu2 |     s1 |    s2 |  alpha |   beta |
|:----------------|-----:|----:|-----:|-----:|------:|------:|-----:|------:|------:|
| \`\`True’’ parameter | 0.600 | 0.30 | 0.100 | 0.500 | -3.000 | -0.050 | 0.100 | 10.000 | 50.000 |
| Parameter estimates  | 0.602 | 0.33 | 0.068 | 0.493 | -2.836 | -0.049 | 0.095 | 11.034 | 54.564 |

Table of the parameter used for simulations (True parameter) and
parameter estimates.

``` r
# tables of AUC of posterior group assignment
AUC = rep(NA, K+1); names(AUC) = paste0("Group G", 0:K)
for(k in 0:K) 
  AUC[k+1] = pROC::roc(simu.dataset$z==k, sorted.res$res$eta[,k+1])$auc
knitr::kable(t(round(AUC,3)), caption = "Table of the Area Under the ROC Curve of posterior group assignment.")
```

| Group G0 | Group G1 | Group G2 |
|---------:|---------:|---------:|
|    0.997 |        1 |    0.989 |

Table of the Area Under the ROC Curve of posterior group assignment.

# Example over a wastewater treatment plant dataset

The folder *data* contains the two wastewater datasets analyzed in \[1\]
in the required format. Each file `data.Ifremer.WWTP1.RData` and
`data.Ifremer.WWTP2.Rdata` contains a dataset as a list composed of a
matrix *x* of mutation counts and a matrix *d* of read depths at related
positions with mutation names in the header column and sampling dates
(as character) in the header row.

``` r
load("data/Ifremer_WWTP1.RData")
head(data.Ifremer.WWTP1$x)
```

    ##      2020-10-20 2020-11-04 2020-11-17 2020-12-04 2020-12-25 2021-01-08
    ## A55T          0          0          1          0          0          0
    ## A55G          0          0          0          0          0          1
    ## G56A          0          0          0          0          0          1
    ## G56C          0          0          0          0          0          0
    ## A57T          0          0          0          0          0          0
    ## T58C          0          0          0          0          0          0
    ##      2021-02-02 2021-02-10 2021-02-23 2021-03-16 2021-03-24 2021-04-06
    ## A55T          0          0          0          0          0          0
    ## A55G          0          0          0          0          0          0
    ## G56A          0          0          0          0          0          0
    ## G56C          1          0          0          0          0          0
    ## A57T          1          0          0          0          0          0
    ## T58C          2          0          1          0          0          0

``` r
head(data.Ifremer.WWTP1$d)
```

    ##      2020-10-20 2020-11-04 2020-11-17 2020-12-04 2020-12-25 2021-01-08
    ## A55T        156        164        138        143         71        156
    ## A55G        156        164        138        143         71        156
    ## G56A        199        195        185        186        100        187
    ## G56C        199        195        185        186        100        187
    ## A57T        256        271        274        264        145        284
    ## T58C        307        318        317        321        180        332
    ##      2021-02-02 2021-02-10 2021-02-23 2021-03-16 2021-03-24 2021-04-06
    ## A55T        138        157        140         95         35        182
    ## A55G        138        157        140         95         35        182
    ## G56A        182        193        184        131         36        229
    ## G56C        182        193        184        131         36        229
    ## A57T        251        289        283        173         57        289
    ## T58C        323        333        318        209         71        321

## Dataset preparation

These two datasets were obtained by running the function
*from_csv_to_count_matrices*, implemented in the file *prepare_data.R*,
over raw data (.tsv output of the variant caller ivar and .txt output
samtools depth) stored in the folder *data_raw*. The package *readr* is
required to run the function *from_csv_to_count_matrices*.

``` r
source("R_code/prepare_data.R")
all.samples = read.table(file = "data_raw/sample_list.txt") # metadata of sample list
data.Ifremer.WWTP1 = from_csv_to_count_matrices(lsdt = all.samples[all.samples$Sampling.site == "WWTP1",], 
                                                path.to.csv = "data_raw/",
                                                length.seq = 29903)
data.Ifremer.WWTP2 = from_csv_to_count_matrices(lsdt = all.samples[all.samples$Sampling.site == "WWTP2",], 
                                                path.to.csv = "data_raw/",
                                                length.seq = 29903)
```

The function *reduce_data* in file *prepare_data.R* takes a dataset in
required format and reduce it according to requested time period and
thresholds.

``` r
source("R_code/prepare_data.R")
load("data/Ifremer_WWTP1.RData")
load("data/profile_matrix.RData")
data = reduce_data(data = data.Ifremer.WWTP1,
                   time.period = c("2020-11-04", "2020-11-17"), 
                   threshold.freq = 0.05, 
                   prop = NULL, 
                   threshold.depth = 10, 
                   reduce.profile = list(lineage = c("B.1", "B.1.1", "B.1.160"),
                                         threshold = 0.005,
                                         mutation.profile = profile.matrix)
                   )
```

This function takes as input:

-   `data` (mandatory): a dataset in the required format, as a list
    composed of a matrix *x* of mutation counts and a matrix *d* of read
    depths at related positions.
-   `time.period` (default: NULL): a vector composed of the starting and
    the ending date of the analysis. If NULL, the entire time period of
    the dataset entered as input `data` is used.
-   `threshold.freq` (default: 0.05): the minimal frequency of a
    mutation to be found in a proportion `prop` (default: NULL) of the
    samples. If `prop` is NULL, `threshold.freq` is required in at least
    one sample among all.
-   `threshold.depth` (default: 10): the threshold for read depth. Any
    read depth below that value in matrix `data$d` is set to zero along
    with associated mutation count in `data$x`.
-   `reduce.profile` (default: NULL): a list composed of a vector
    `lineage` of names of lineages, a single value `threshold` and a
    mutation profile matrix `mutation.profile`. Input `mutation.profile`
    is a mutations × lineages matrix containing the probability of
    mutations to belong to lineages. Mutations associated to a
    probability strictly below `reduce.profile$threshold` to belong to,
    at least, one lineage among `reduce.profile$lineage` are removed. If
    `reduce.profile` is NULL (default value), there is no reduction on
    mutation profile.

Note that for an unsupervised clustering, input `reduce.profile$lineage`
must only contain lineages known before the time period studied. The
file *profile_matrix.RData* in folder *data* contains a mutation profile
matrix, named *profile.matrix*, computed with Virpool package \[2\] (see
\[1\] for details).

## Selection of the number of groups

As mentioned in \[1\], one may compute the BIC or ICL of models composed
of various number of non-neutral groups *K* ≥ 0 using the output
`loglik` of the function *ft_clust* (log-likelihood of the model) and
use it as criteria to select the number of groups in the dataset.
However, this approach is inefficient on real datasets (wastewater
datasets) and an alternative strategy based on the size of the smallest
group is proposed and proved to be efficient on several wasterwater
datasets. One may therefore run the function *ft_clust* over the chosen
dataset conditional on various numbers of non-neutral groups *K* and
retrieve to output `size.smallest.group`.

``` r
size.of.the.smallest.group = rep(NA, 6);
names(size.of.the.smallest.group) = paste0("K=", 1:6)
for (K in 1:6) {
  tmp = ft_clust(data, K, n.init = c(40,40), n.it = c(30,30), verbose=FALSE)
  size.of.the.smallest.group[K] = tmp$size.smallest.group
}
knitr::kable(t(size.of.the.smallest.group), caption = "Size of the smallest group for increasing number of non-neutral groups.")
```

| K=1 | K=2 | K=3 | K=4 | K=5 | K=6 |
|----:|----:|----:|----:|----:|----:|
|   8 |   6 |   1 |   1 |   1 |   1 |

Size of the smallest group for increasing number of non-neutral groups.

For a computational boost, one may use parallel computation over `K`
and/or over `n.init`.

## Parameter estimates for WWTP1 dataset restricted to time period 2020-11-04 to 2020-11-17

``` r
res = ft_clust(data = data, K = 2, n.init = c(10,10), n.it = c(5,5), verbose = FALSE)

# sort groups in decreasing selection coefficients 
idx = sort.int(res$s, decreasing = TRUE, index.return = TRUE)$ix
pi = res$pi[c(1,idx+1)]; mu = res$mu[idx]; s = res$s[idx]; eta = res$eta[,c(1,idx+1)]
```

``` r
# print table of parameters
p = res$alpha / (res$alpha + res$beta) # frequency of the neutral group
hat.par = rbind(pi = pi, mu = c(log(p/(1-p)),mu), s = c(0,s))
dimnames(hat.par)[[2]] = paste0("Group ", 0:length(res$mu))
knitr::kable(round(hat.par, 3), caption = "Parameter estimates for WWTP1 dataset restricted to time period 2020-11-04 to 2020-11-17 conditional on K=2 non-neutral groups.")
```

|     | Group 0 | Group 1 | Group 2 |
|:----|--------:|--------:|--------:|
| pi  |   0.595 |   0.189 |   0.216 |
| mu  |  -0.009 |  -3.803 |   1.614 |
| s   |   0.000 |   0.258 |  -0.187 |

Parameter estimates for WWTP1 dataset restricted to time period
2020-11-04 to 2020-11-17 conditional on K=2 non-neutral groups.

``` r
# graphical representation of group frequency trajectories
# # vector of colors according to increasing and decreasing groups 
col = list(increasing = c("red", "darkorange", "magenta", "blue"), 
           decreasing = c("green", "cyan", "darkolivegreen", "pink", "yellow", "gray"))
col = c(col$increasing[1:sum(s>0)], col$decreasing[1:sum(s<0)])
# # retrieve sampling date for x-axis 
sampling.date = as.Date(colnames(data$x))

# plot 
K = length(res$mu)
par(mar=c(5.1+2, 5.1, 5.1, 2.1), xpd=TRUE) # for inset in legend
p = res$alpha / (res$alpha + res$beta) # frequency of the neutral group
plot(sampling.date, rep(p, length(sampling.date)), xaxt = "n", yaxt = "n",xlab="", ylab="mutation frequency",t="l", ylim=c(0,1), cex.lab = 1, col = "black", lwd = 1.4)
for (k in 1:K) {
  # compute frequency at different time point using the time differences between sampling date and date of origin
  p = mu[k] + s[k]*as.numeric(julian(sampling.date, origin = sampling.date[1])) # 
  lines(sampling.date, exp(p)/(1+exp(p)), col = col[k], lwd = 1.4)
}
axis(side = 2, las = 2, mgp = c(3, 0.75, 0))
text(x = sampling.date,
     y = exp(par("usr")[3])/(1+exp(par("usr")[3])) -0.58, #position on the y axis of the text
     adj=1,
     labels = as.character(sampling.date),
     xpd = NA,
     srt = 45, #rotate labels of ... degrees
     cex = 1)
axis(side = 1, at = sampling.date, labels = rep("",length(sampling.date)), tick=TRUE)
legend("topright", legend = paste0("G", (0:K)), lty=1, col = c("black", col), ncol = K+1, cex = 1)
```

<figure>
<img
src="README_files/figure-markdown_github/wwtp1%20restricted%20plot-1.png"
alt="Grapical representation of estimated group frequency trajectories in WWTP1 dataset restricted to time period 2020-11-04 to 2020-11-17 conditional on K=2 non-neutral groups." />
<figcaption aria-hidden="true">Grapical representation of estimated
group frequency trajectories in WWTP1 dataset restricted to time period
2020-11-04 to 2020-11-17 conditional on K=2 non-neutral
groups.</figcaption>
</figure>

## Parameter estimates for WWTP1 dataset over its entire time period

``` r
source("R_code/prepare_data.R")
source("R_code/EM_algorithm.R")
load("data/Ifremer_WWTP1.RData")
data = reduce_data(data = data.Ifremer.WWTP1,
                    time.period = NULL, 
                    threshold.freq = 0.05, 
                    prop = 1/4, 
                    threshold.depth = 10, 
                    reduce.profile = NULL
                    )
res = ft_clust(data = data, K = 6, n.init = c(40,40), n.it = c(5,5), verbose = FALSE)

# sort groups in decreasing selection coefficients 
idx = sort.int(res$s, decreasing = TRUE, index.return = TRUE)$ix
pi = res$pi[c(1,idx+1)]; mu = res$mu[idx]; s = res$s[idx]; eta = res$eta[,c(1,idx+1)]
```

|     | Group 0 | Group 1 | Group 2 | Group 3 | Group 4 | Group 5 | Group 6 |
|:----|--------:|--------:|--------:|--------:|--------:|--------:|--------:|
| pi  |   0.666 |   0.084 |   0.039 |   0.039 |   0.047 |   0.072 |   0.054 |
| mu  |  -1.518 |  -4.284 |  -5.967 |  -2.787 |   1.721 |  -0.861 |   0.933 |
| s   |   0.000 |   0.035 |   0.030 |   0.028 |  -0.017 |  -0.019 |  -0.021 |

Parameter estimates for WWTP1 dataset over its entire rime period
conditional on K=6 non-neutral groups.

<figure>
<img
src="README_files/figure-markdown_github/wwtp1%20entire%20plot-1.png"
alt="Estimated group frequency trajectories in WWTP1 dataset over its entire time period conditional on K=6 non-neutral groups." />
<figcaption aria-hidden="true">Estimated group frequency trajectories in
WWTP1 dataset over its entire time period conditional on K=6 non-neutral
groups.</figcaption>
</figure>

And the list of mutations composing each group (except the neutral one
that contains most mutations):

``` r
clust = vector("list", K+1); names(clust) = paste("Group", 0:K)
for (k in 0:K) {
  clust[[k+1]] = rownames(eta)[(apply(eta, 1, which.max)==k+1)]
}
print(clust[-1])
```

    ## $`Group 1`
    ##  [1] "C913T"   "C3267T"  "C5388A"  "C14676T" "C15279T" "A23063T" "C23604A"
    ##  [8] "G24914C" "C27972T" "G28048T" "A28281T" "G28882A" "G28883C"
    ## 
    ## $`Group 2`
    ## [1] "C2453T"  "G4136T"  "C6027T"  "T11296G" "C24418T" "T24506G"
    ## 
    ## $`Group 3`
    ## [1] "T16176C" "C23271A" "C23709T" "A28111G" "G28280C" "T28282A"
    ## 
    ## $`Group 4`
    ## [1] "C11497T" "G13993T" "G15766T" "A16889G" "G17019T" "C25710T" "C26735T"
    ## 
    ## $`Group 5`
    ##  [1] "A1087C"  "G12602T" "G13723T" "G15768T" "C22227T" "G23678T" "G25599T"
    ##  [8] "C26801G" "G26803T" "C27434T" "G29645T"
    ## 
    ## $`Group 6`
    ## [1] "C4543T"  "G5629T"  "G9526T"  "C18877T" "G22992A" "G25563T" "T26876C"
    ## [8] "G28975C" "G29399A"

# References

\[1\] *Unsupervised detection and fitness estimation of emerging
SARS-CoV-2 variants: Application to wastewater samples (ANRS0160)*
(Lefebvre, A. et al., 2025)
〈[arxiv-2501.06548](https://arxiv.org/abs/2501.06548)〉

\[2\] *VirPool: model-based estimation of SARS-CoV-2 variant proportions
in wastewater samples* (Gafurov, A. et al., 2022)
〈[https://doi.org/10.1186/s12859-022-05100-3](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05100-3)〉
