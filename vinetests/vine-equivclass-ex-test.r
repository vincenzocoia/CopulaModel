# This code shows how sometimes there are two natural order representations
# of a vine equivalence class;
# reference is 
# Joe H, Cooke RM and Kurowicka D (2011). Regular vines: generation algorithm 
# and number of equivalence classes.
# In Dependence Modeling: Vine Copula Handbook, pp 219--231.
# World Scientific, Singapore.

library(CopulaModel)

d=7
perm=c(6,4,3,1,5,7,2) 
A=vnum2array(d,320)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    2    1    2    5
#[2,]    0    2    2    1    3    4    1
#[3,]    0    0    3    3    2    1    3
#[4,]    0    0    0    4    4    3    2
#[5,]    0    0    0    0    5    5    4
#[6,]    0    0    0    0    0    6    6
#[7,]    0    0    0    0    0    0    7

Aperm=varrayperm(A,perm)
print(Aperm)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    6    6    6    4    6    4    5
#[2,]    0    4    4    6    3    1    6
#[3,]    0    0    3    3    4    6    3
#[4,]    0    0    0    1    1    3    4
#[5,]    0    0    0    0    5    5    1
#[6,]    0    0    0    0    0    7    7
#[7,]    0    0    0    0    0    0    2

ANO1=varray2NO(Aperm)
ANO2=varray2NO(Aperm,irev=T) # the second form is equiv class has 2 NO arrays
print(ANO1) # $NO is same as above
#$NOa
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    6    6    6    4    6    4    5
#[2,]    0    4    4    6    3    1    6
#[3,]    0    0    3    3    4    6    3
#[4,]    0    0    0    1    1    3    4
#[5,]    0    0    0    0    5    5    1
#[6,]    0    0    0    0    0    7    7
#[7,]    0    0    0    0    0    0    2
#$NO
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    2    1    2    5
#[2,]    0    2    2    1    3    4    1
#[3,]    0    0    3    3    2    1    3
#[4,]    0    0    0    4    4    3    2
#[5,]    0    0    0    0    5    5    4
#[6,]    0    0    0    0    0    6    6
#[7,]    0    0    0    0    0    0    7
#$perm
#[1] 4 7 3 2 5 1 6
#$diag
#[1] 6 4 3 1 5 7 2

print(ANO2)
#$NOa
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    6    6    6    6    4    5    4
#[2,]    0    3    3    3    6    6    1
#[3,]    0    0    4    4    3    3    6
#[4,]    0    0    0    5    5    4    3
#[5,]    0    0    0    0    1    1    5
#[6,]    0    0    0    0    0    2    2
#[7,]    0    0    0    0    0    0    7
#$NO
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    1    3    4    3
#[2,]    0    2    2    2    1    1    5
#[3,]    0    0    3    3    2    2    1
#[4,]    0    0    0    4    4    3    2
#[5,]    0    0    0    0    5    5    4
#[6,]    0    0    0    0    0    6    6
#[7,]    0    0    0    0    0    0    7
#$perm
#[1] 5 6 2 3 4 1 7
#$diag
#[1] 6 3 4 5 1 2 7

# binary representations are different
bin1=varray2bin(ANO1$NO)
bin2=varray2bin(ANO2$NO)
print(bin1)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    1    1    1    1
#[2,]   NA    1    1    0    1    1    0
#[3,]   NA   NA    1    1    0    0    0
#[4,]   NA   NA   NA    1    1    0    0
#[5,]   NA   NA   NA   NA    1    1    0
#[6,]   NA   NA   NA   NA   NA    1    1
#[7,]   NA   NA   NA   NA   NA   NA    1
print(bin2)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    1    1    1    1
#[2,]   NA    1    1    1    0    0    1
#[3,]   NA   NA    1    1    0    0    0
#[4,]   NA   NA   NA    1    1    0    0
#[5,]   NA   NA   NA   NA    1    1    0
#[6,]   NA   NA   NA   NA   NA    1    1
#[7,]   NA   NA   NA   NA   NA   NA    1

# get the bnum decimal representation
bvec1=NULL
for (i in 4:d) 
{ bvec1=c(bvec1,bin1[2:(i - 2), i]) }
# 0 1 0 1 0 0 0 0 0 0
# 256+256/4=320
dec1=v2d(bvec1,2)
print(dec1)

bvec2=NULL
for (i in 4:d) 
{ bvec2=c(bvec2,bin2[2:(i - 2), i]) }
#  1 0 0 0 0 0 1 0 0 0
# 512+8=520
dec2=v2d(bvec2,2)
print(dec2)

vnum2array(7,520)
# same as ANO2$NO)
