# check if given A arrays are valid vine arrays

library(CopulaModel)

# simplest 4-dimemsional non-vine array with 1:d on diagonal
A=matrix(0,4,4)
diag(A)=1:4
A[1,2:4]=c(1,1,2)
A[2,3:4]=c(2,3)
A[3,4]=1
print(varraycheck(A)) # -1

A1=matrix(0,5,5)
diag(A1)=1:5
A1[1,2:4]=1
A1[2,3:4]=2
A1[3,4]=3
A1[4,5]=4
A2=A1
A2[1:2,4]=c(2,1)

# check 2 improper arrays for first 3 entries of column 5 of A1 and A2
A1[1:3,5]=c(2,3,1)
print(varraycheck(A1)) # -1
A1[1:3,5]=c(3,2,1)
print(varraycheck(A1)) # -1

A2[1:3,5]=c(2,3,1)
print(varraycheck(A2)) # -1
A2[1:3,5]=c(3,2,1)
print(varraycheck(A2)) # -1

A2=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,3,2,4,0,0, 4,3,2,1,5,0,
        5,3,4,2,1,6),6,6)
b2=varraycheck(A2)
print(b2)  # -1

