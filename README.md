# [An R package `PIC2`](https://github.com/hfang-bristol/PIC2)

## @ Overview

> The `PIC2` is an R package enabling a priority index solution to critical Covid-19 .


## @ Installation

### 1. Install R

Please install R (version 4.3.2 or above); see https://cran.r-project.org

If installed on `Ubuntu` (assuming you have a `ROOT (sudo)` privilege), please do so below

```ruby
sudo su
# here enter your password

wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.3.2.tar.gz
tar xvfz R-4.3.2.tar.gz
cd ~/R-4.3.2
./configure
make
make check
make install
R # start R
```

### 2. Install R packages

```ruby
R # start R

# if the package 'BiocManager' not installed, please do so
if(!("BiocManager" %in% rownames(installed.packages()))) install.packages("BiocManager")

# first, install basic packages: remotes, tidyverse
BiocManager::install(c('remotes','tidyverse'), dependencies=T)

# then, install the package 'PIC2' (now hosted at github)
BiocManager::install("hfang-bristol/PIC2", dependencies=T, force=T)

# check the package 'PIC2' successfully installed
library(PIC2)
```


## @ Showcase

> A showcase, describing a step-by-step protocol on how to generate the target index and the drug index for critical Covid-19, is made reproducible [`here`](http://www.genetictargets.com/PIC2/showcase).


## @ Contact

> Please drop [email](mailto:fh12355@rjh.com.cn) for bug reports or enquiries.


