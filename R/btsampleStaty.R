# Functions for stationary bootstrap of Politis and Romano 1994. JASA
# Code written by Pavel Krupskii

##########################################################################

# n = length of series
# p = geometric rate such as n^(-1/3)
# Output: single random block
btblock= function(n,p)
{ ell = floor(n*runif(1))+1;
  b = rgeom(1,p)+1;
  if(ell+b-1 <= n) block0 = ell:(ell+b-1);
  if(ell+b-1 >  n) 
  { block1 = ell:n; block2 = 1:(ell+b-1-n); block0 = c(block1,block2); }
  block0
}

# sample blocks repeated until desired length is reached
# size = length of time series
# p = geometric rate such as n^(-1/3)
# Output: vector of length 'size', each element is an integer between 1 and 'size',
#  with indices of the bootstrapped sample 
btsampleStaty=function(size,p)
{ ib = 0;
  isample = NULL;
  while(ib < size)
  { btgrp = btblock(size,p);
    isample = c(isample, btgrp);
    ib = length(isample); 
    if(ib >= size) isample = isample[1:size]; 
  }
  isample
}

