# CopulaModel

Repository for the `CopulaModel` R package by Harry Joe and Pavel Krupskii.

## Source of content

The content in this repository is originally located at https://copula.stat.ubc.ca/, under the heading "CopulaModel software". It's made available on GitHub in an effort to improve accessibility.

## Version

The files you see here integrate the 2015.09.03 "patches and additions" with the original `CopulaModel` package.

## Installing the R package

:warning: These instructions have not been vetted yet :warning:. 

### Prerequisites

You'll need a couple things installed before installing `CopulaModel` from this GitHub repo:

1. gcc
2. devtools

#### 1\. gcc

For Mac OS, you can use HomeBrew:

1. Open Terminal
2. Check that you have HomeBrew installed by typing `which brew`. 
  - If it returns a path, such as `/usr/local/bin/brew`, then it's installed.
  - If it says `brew not found`, you'll have to [install HomeBrew](https://brew.sh/).
3. Install gcc by typing `brew install gcc`
  - If you get an error, hopefully HomeBrew will tell you how to fix it.
4. You can check that gcc is successfully isntalled by typing `which gcc`.

Sorry, I don't have Windows instructions yet. I think Linux should work the same way as Mac OS here.

#### 2\. devtools

Run the following R code to install the `devtools` package:

```
install.packages("devtools")
```

### Installing `CopulaModel`

Run the following R code to install the `CopulaModel` package:

```
devtools::install_github("vincenzocoia/CopulaModel")
```
