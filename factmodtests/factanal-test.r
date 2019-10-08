# check factanal for number of factors for best fit

library(CopulaModel)

data(euro07gf)
r07=cor(euro07gf$zscore)

#          OSEAX      FTSE       AEX      FCHI      SSMI     GDAXI       ATX
#OSEAX 1.0000000 0.6872951 0.7367698 0.7151750 0.6606776 0.7155215 0.6691907
#FTSE  0.6872951 1.0000000 0.8913743 0.9073861 0.8332023 0.8902659 0.7618218
#AEX   0.7367698 0.8913743 1.0000000 0.9186479 0.8454261 0.8898639 0.7675997
#FCHI  0.7151750 0.9073861 0.9186479 1.0000000 0.8789436 0.9249947 0.7706429
#SSMI  0.6606776 0.8332023 0.8454261 0.8789436 1.0000000 0.8461359 0.7458799
#GDAXI 0.7155215 0.8902659 0.8898639 0.9249947 0.8461359 1.0000000 0.7934639
#ATX   0.6691907 0.7618218 0.7675997 0.7706429 0.7458799 0.7934639 1.0000000

# factanal

# 7-par fit to 21 saturated
fa1=factanal(covmat=r07,factors=1)
tem1=as.matrix(fa1$loadings[,1])
rfact1=tem1%*%t(tem1)+diag(fa1$uniq)
cat("max abs difference1 =", max(abs(r07-rfact1)),"\n")
cat("corDis1=", corDis(rfact1,r07),"\n")
#max abs difference1 = 0.05753447 
#corDis1= 0.1145649 

# 13-par fit to 21 saturated
fa2=factanal(covmat=r07,factors=2)
tem2=as.matrix(fa2$loadings[,1:2])
rfact2=tem2%*%t(tem2)+diag(fa2$uniq)
cat("max abs difference2 =", max(abs(r07-rfact2)),"\n")
cat("corDis2=", corDis(rfact2,r07),"\n")
#max abs difference2 = 0.02526525 
#corDis2= 0.05096038 

# 18-par fit to 21 saturated
fa3=factanal(covmat=r07,factors=3)
tem3=as.matrix(fa3$loadings[,1:3])
rfact3=tem3%*%t(tem3)+diag(fa3$uniq)
cat("max abs difference3 =", max(abs(r07-rfact3)),"\n")
cat("corDis3=", corDis(rfact3,r07),"\n")
#max abs difference3 = 0.008528348 
#corDis1= 0.0121464 


