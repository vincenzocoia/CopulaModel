# checks for conversions from vine array to binary representation and
# vice versa

library(CopulaModel)

d=5
d=6
d=7
d=8
bnum=2^((d-2)*(d-3)/2)-1  # largest bnum is for C-vine
cat("d=", d, "  #vine arrays in natural order =", bnum+1,"\n")
b10=floor(bnum/10)

# enumerate through vine arrays based on binary code
# only do this for up to d=8
cat("check that varray2bin and vbin2array are inverses\n")
for(b in 0:bnum)
{ A=vnum2array(d,b)
  bin=varray2bin(A) # if A is valid, this is dxd array with NA in lower triangle
  AA=vbin2array(d,bin) # should be same as A
  tem=max(abs(A-AA))
  if(b%%b10==0) cat("checking ", b,"\n")
  if(tem!=0) cat(b,tem,"\n") 
}
cat("done\n")

