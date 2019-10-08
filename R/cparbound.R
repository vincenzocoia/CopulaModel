# Function to show domain of copula parameters with input of
# abbreviated copula name used in the pcop, dcop, pcondcop functions.

# e.g., input is gum or bb1 or bvncop etc

#copparbounds=read.table("../data/copparbounds.tab",header=T)

# Bounds on copula parameters for bivariate copula families included
# in this library
# copname = string for abbreviated copula name used in pcop, dcop, pcondcop
# Output: line number in copparbounds.tab with bounds of the copula parameter(s)
#   the lower and upper bounds are printed out
cparbound=function(copname)
{ data(copparbounds)
  icop=try(which(copname==copparbounds$copname),silent=T)
  if(length(icop)==0)
  { cat("copname ", copname, " not included\n"); return(NA) }
  else
  { #print(icop) 
    cat(copname,"\n")
    cat("parameter 1: ")
    cat("lower bound is ", copparbounds$lb1[icop])
    cat("; upper bound is ", copparbounds$ub1[icop],"\n")
    if(!is.na(copparbounds$lb2[icop]))
    { cat("parameter 2: ")
      cat("lower bound is ", copparbounds$lb2[icop])
      cat("; upper bound is ", copparbounds$ub2[icop],"\n")
    }
    if(!is.na(copparbounds$lb3[icop]))
    { cat("parameter 3: ")
      cat("lower bound is ", copparbounds$lb3[icop])
      cat("; upper bound is ", copparbounds$ub3[icop],"\n")
    }
    invisible(icop)
  }
}

#cparbound("bvncop")
#cparbound("pla")
#cparbound("gum")
#cparbound("bb1")
#cparbound("bb1rpow")
#cparbound("none")
