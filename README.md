# CopulaModel

Repository for the `CopulaModel` R package by Harry Joe and Pavel Krupskii.

## Source of content

The content in this repository is originally located at https://copula.stat.ubc.ca/, under the heading "CopulaModel software". It's made available on GitHub in an effort to improve accessibility.

## Version

The files you see here integrate the 2015.09.03 "patches and additions" with the original `CopulaModel` package.

## Installing the R package

:warning: These instructions have not been vetted yet :warning:. 

If you have the prerequisite software installed (see below), executing the following R code should install the `CopulaModel` package:

```
devtools::install_github("vincenzocoia/CopulaModel")
```

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

## Installation Troubleshooting

### "Library not loaded"

When tested on a Mac, the installation tries to look for "libraries" in the wrong place. You can tell if this is the case if you find the message "Library not loaded" somewhere in the error message. I figured out how to solve the problem, thanks to [this Stack Overflow post](https://stackoverflow.com/a/57225398). Here's what I ended up having to do.

The libraries that were missing on my mac were:

- `libgfortran.3.dylib`, although was present on my mac as `libgfortran.5.dylib` (a newer version?) 
    - R was looking for it at `/usr/local/gfortran/lib/libgfortran.3.dylib`
    - I found the file at `/usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libgfortran.5.dylib`
- `libquadmath.0.dylib`
    - R was looking for it at `/usr/local/gfortran/lib/libquadmath.0.dylib`
    - I found the file at `/usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libquadmath.0.dylib`
- `libR.dylib`
    - R was looking for it at `/Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libR.dylib`
    - I found the file at `/Library/Frameworks/R.framework/Versions/Current/Resources/lib/libR.dylib`

After finding the actual location of these files on my mac, I had to link them to where R was looking. I found the first two thanks to the stack overflow post above (although their suggested use of `locate` didn't work for me), and the last one I found by updating the R version in the path. Your paths may differ depending on your gcc version.

Linking the files involved using the terminal command `ln`, as in `ln actual/path/to/file path/where/R/is/looking`. You might need to precede all this with `sudo`. And for the last one (`libR.dylib`), I had to make the directory `3.5/Resources/lib/` in `/Library/Frameworks/R.framework/Versions/` using `mkdir`. So, all in all:

```
sudo ln /usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libgfortran.5.dylib /usr/local/gfortran/lib/libgfortran.3.dylib
ln /usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libquadmath.0.dylib /usr/local/gfortran/lib/libquadmath.0.dylib
mkdir -p /Library/Frameworks/R.framework/Versions/3.5/Resources/lib
ln /Library/Frameworks/R.framework/Versions/Current/Resources/lib/libR.dylib /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libR.dylib
```

### "Rank Mismatch"

Sometimes when I tried to install the package, I wouldn't get the "Library not loaded" error, but would get a "Rank mismatch" error in a Fortran script. You can see the error message below. 

To elicit the more desirable "library not loaded" error, I re-downloaded this very `CopulaModel` repository as a zip file, unzipped it, and called `install.packages("~/Downloads/CopulaModel-master", type = "source", repos = NULL)`

```
> devtools::install_github("vincenzocoia/CopulaModel")
Downloading GitHub repo vincenzocoia/CopulaModel@HEAD
✓  checking for file ‘/private/var/folders/dx/zr_sf92j4t95y6z16w_j6b780000gn/T/RtmpGaNFDT/remotes3fda601e7e12/vincenzocoia-CopulaModel-c0c2150/DESCRIPTION’ ...
─  preparing ‘CopulaModel’:
✓  checking DESCRIPTION meta-information ...
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  looking to see if a ‘data/datalist’ file should be added
─  building ‘CopulaModel_0.6.tar.gz’
   
* installing *source* package ‘CopulaModel’ ...
** using staged installation
** libs
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c R_exchmvn.c -o R_exchmvn.o
gfortran -mmacosx-version-min=10.13 -fno-optimize-sibling-calls  -fPIC  -Wall -g -O2  -c  Rgauss-trvine-nonuniq.f90 -o Rgauss-trvine-nonuniq.o
Rgauss-trvine-nonuniq.f90:107:43:

  103 |                 call intpr("within eps for ",15,j,1)
      |                                                2
......
  107 |                   call intpr("next row",8, A2(i,:),d)
      |                                           1
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)
Rgauss-trvine-nonuniq.f90:124:34:

  103 |                 call intpr("within eps for ",15,j,1)
      |                                                2
......
  124 |               call intpr("perm",4,aperm,d)
      |                                  1
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)
Rgauss-trvine-nonuniq.f90:154:33:

  103 |                 call intpr("within eps for ",15,j,1)
      |                                                2
......
  154 |         call intpr("next row",8, A1(i,:),d)
      |                                 1
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)
Rgauss-trvine-nonuniq.f90:164:51:

  103 |                 call intpr("within eps for ",15,j,1)
      |                                                2   
......
  164 |   call intpr("approx nonunique for j=1,...,d-2",32,nonuniq,d-2)
      |                                                   1
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)
make: *** [Rgauss-trvine-nonuniq.o] Error 1
ERROR: compilation failed for package ‘CopulaModel’
* removing ‘/Library/Frameworks/R.framework/Versions/4.0/Resources/library/CopulaModel’
* restoring previous ‘/Library/Frameworks/R.framework/Versions/4.0/Resources/library/CopulaModel’
Error: Failed to install 'CopulaModel' from GitHub:
  (converted from warning) installation of package ‘/var/folders/dx/zr_sf92j4t95y6z16w_j6b780000gn/T//RtmpGaNFDT/file3fda74f571f1/CopulaModel_0.6.tar.gz’ had non-zero exit status
```
