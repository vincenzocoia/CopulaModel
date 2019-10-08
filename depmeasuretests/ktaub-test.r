# check of taub for 2 ordinal variables 

library(CopulaModel)

#n=200
#n=20
n=100
set.seed(123)
x1=floor(runif(n,0,3))
x2=floor(runif(n,0,3))

otab=table(x1,x2) # for used with taub (Kendall taub for 2-way table)
print(otab)
#   x2
#x1   0  1  2
#  0  7 14 12
#  1 14 14  6
#  2 10 12 11


tau1=cor(x1,x2,method="kendall")
tau2=taucor(x1,x2)  
tau3=taub(otab)    
cat("cor, taucor and taub\n")
cat(tau1,tau2,tau3,"\n")
#cor, taucor and taub
#-0.05820342 -0.05820342 -0.05820342 
