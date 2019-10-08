# Sequential MST based on partial correlations for best truncated Gaussian vine 
# This depends on the R package igraph0 (e.g., version 0.5.6)
# for the minimum spanning tree algorithm and operations on graphs.

# This is a modification of code written by Eike Brechmann.

# The following functions are not exported
# initFirstG = function(rmat,iprint=F)
# fitFirstT = function(mst,iprint=F)
# buildNextG = function(oldVineGraph,rmat,iprint=F)
# fitT = function(mst,oldVineGraph,progress = FALSE)
# asRVA = function(RVine)

# rmat = correlation matrix
# ntrunc = upper bound in truncation level to consider
# iprint = print flag for intermediate results
# Output: 
#  $RVM with 
#  $RVM$VineA  vine array 
#  $RVM$pc  partial correlations by tree
#  $RVM$Matrix vine array in VineCopula format [d:1,d:1]
#  $RVM$Cor   partial correlations in VineCopula format [d:1,d:1]
#  $mst with
#  $mst[[1]]  (first tree)
#  ...
#  $mst[[d-1]]  (last tree where d=dimension)
#  $treeweight  vector of length d-1 with 
#                 sum_edge log(1-rho[edge]^2) for trees 1,...d-1
#  $trunclevel  same as inputted ntrunc
#  $truncval = sum_{1:trunclevel} treeweight/ sum_{1:(d-1)} treeweight
gausstrvine.mst=function(rmat,ntrunc,iprint=F)
{ rmat=as.matrix(rmat)
  d=dim(rmat)[1]
  colnames(rmat)= rownames(rmat) = paste("V",1:d,sep="")
  RVine=list(Tree=NULL, Graph=NULL)
  mst=list()
  treeweight=rep(NA,d-1)
  logdet=log(det(rmat))
  truncated=FALSE
  
  # tree 1
  g= initFirstG(rmat)
  mst[[1]] = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
  treeweight[1] = sum(log(1-E(mst[[1]])$weight^2))
  fitval= treeweight[1]/logdet
  if(ntrunc==1)
  { trunclevel=1
    truncval=fitval
    truncated=TRUE
  }
  
  VineTree= fitFirstT(mst[[1]])
  RVine$Tree[[1]] = VineTree
  RVine$Graph[[1]] = g
  
  # trees 2:(d-1)
  for(i in 2:(d-1))
  { g= buildNextG(VineTree,rmat,iprint=iprint)
    mst[[i]] = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
    treeweight[i] = sum(log(1-E(mst[[i]])$weight^2))
    fitval= sum(treeweight[1:i])/logdet 
    if(!truncated)
    { if(ntrunc==i)
      { trunclevel=i
        truncval=fitval
        truncated=TRUE
      }
    }
    VineTree= fitT(mst[[i]],VineTree,progress=F)
    dd= ecount(VineTree)
    RVine$Tree[[i]] = VineTree
    RVine$Graph[[i]] = g
  }
  
  RVA=asRVA(RVine)
  return(list(RVM=RVA,mst=mst,treeweight=treeweight,trunclevel=trunclevel,truncval=truncval))
}

# initialize the first graph
initFirstG= function(rmat,iprint=F)
{ g= graph.adjacency(rmat, mode="lower",weighted=TRUE,diag=FALSE)
  E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep=",")
  for(i in 1:ecount(g)) E(g)$conditionedSet[[i]] = get.edges(g,i)
  # for(i in 1:ecount(g)) E(g)$conditionedSet[[i]] = get.edges(g,i-1)
  if(iprint) print(g)
  return(g)
}

# Fit the first tree of the vine.
# This function adds attribute $CondName1, $CondName2, $Copula.CondName
# but otherwise doesn't change the graph
fitFirstT= function(mst,iprint=F)
{ d=ecount(mst)
  for(i in 1:d)
  { a= get.edges(mst,i)
    # a= get.edges(mst,i-1)+1
    # starts as  NULL so simplify, $name comes from initFirstG 
    E(mst)[i]$CondName1 = V(mst)[a[1]]$name
    E(mst)[i]$CondName2 = V(mst)[a[2]]$name
    # E(mst)[i-1]$CondName1 = V(mst)[a[1]-1]$name
    # E(mst)[i-1]$CondName2 = V(mst)[a[2]-1]$name
  }
  if(iprint) 
  { print(E(mst)[1]$CondName1); print(E(mst)[2]$CondName1)
    print(E(mst)[1]$CondName2); print(E(mst)[2]$CondName2)
    # print(E(mst)[0]$CondName1); print(E(mst)[1]$CondName1)
    # print(E(mst)[0]$CondName2); print(E(mst)[1]$CondName2)
  }
  return(mst)
}

# Build next tree with vertices that are edges of previous tree
buildNextG= function(oldVineGraph,rmat,iprint=F)
{ EL= get.edgelist(oldVineGraph)
  dd= ecount(oldVineGraph) # number of edges in previous tree
  if(iprint) cat("\n***ecount=", dd,"\n")
  g= graph.full(dd) # all posible edges
  V(g)$name = E(oldVineGraph)$name # previous graph had dd edges
  V(g)$conditionedSet = E(oldVineGraph)$conditionedSet
  if(!is.null(E(oldVineGraph)$conditioningSet)) 
  { V(g)$conditioningSet = E(oldVineGraph)$conditioningSet }
  
  # for(i in 0:(ecount(g)-1))
  for(i in 1:ecount(g)) # g is full graph so it goes thru all pairs
  { con= get.edge(g,i) # internal indices for edge
    temp= get.edges(oldVineGraph,con) # edges of prev tree 
    # this is a ways of pairing up nodes (which were prev edges)
    temp2=unique(c(temp))
    ok=(length(temp2)==3) 
    if(ok)
    { same=intersect(temp[1,],temp[2,])
      other1=setdiff(temp[1,],same)
      other2=setdiff(temp[2,],same)
      E(g)[i]$nodes = paste(as.character(c(other1,other2,same)),collapse=",")
      # E(g)[i]$nodes =paste(as.character(c(other1,other2,same)+1),collapse=",")
      name.node1= strsplit( V(g)[con[1]]$name,split=" *[,|] *")[[1]]
      name.node2= strsplit( V(g)[con[2]]$name,split=" *[,|] *")[[1]]
      schnitt=intersect(name.node1,name.node2)
      symdiff=setdiff(union(name.node1,name.node2),schnitt)
      # new edge
      E(g)[i]$name = paste(paste(symdiff, collapse= ","),paste(schnitt, collapse= ","),sep= " | ")
      l1= c(unlist(V(g)[con[1]]$conditionedSet),unlist(V(g)[con[1]]$conditioningSet))
      l2= c(unlist(V(g)[con[2]]$conditionedSet),unlist(V(g)[con[2]]$conditioningSet))
      # l1= c(V(g)[con[1]]$conditionedSet,V(g)[con[1]]$conditioningSet)
      # l2= c(V(g)[con[2]]$conditionedSet,V(g)[con[2]]$conditioningSet)
      schnitt=intersect(l1,l2)
      symdiff=setdiff(union(l1,l2),schnitt)
      suppressWarnings({E(g)$conditionedSet[i] = list(symdiff)})
      suppressWarnings({E(g)$conditioningSet[i] = list(schnitt)})
      # suppressWarnings({E(g)$conditionedSet[i+1] = list(symdiff)})
      # suppressWarnings({E(g)$conditioningSet[i+1] = list(schnitt)})
      conditioned=symdiff  # size 2 
      conditioning=schnitt # given
      # conditioned=symdiff+1  # size 2
      # conditioning=schnitt+1 # given
      needed= sort(c(conditioned,conditioning))
      ind= which(needed %in% conditioned)
      pp=partcor(rmat,conditioning,conditioned[1],conditioned[2]) 
      E(g)[i]$weight = pp
    }
    E(g)[i]$todel = !ok
  }
  
  g= delete.edges(g, E(g)[E(g)$todel])
  if(iprint) print(g)
  return(g)
}

# Fit tree
fitT= function(mst,oldVineGraph,progress=FALSE)
{ dd=ecount(mst)
  # for(i in 0:(dd-1))
  for(i in 1:dd)
  { con= get.edge(mst,i)
    temp= get.edges(oldVineGraph,con)
    same=intersect(temp[1,],temp[2,])
    # n1, n2 are the variable names in the symmetric differnce
    if(temp[1,1]==same)
    { n1= E(oldVineGraph)[con[1]]$CondName2 }
    else
    { n1= E(oldVineGraph)[con[1]]$CondName1 }
    if(temp[2,1]==same)
    { n2= E(oldVineGraph)[con[2]]$CondName2 }
    else
    { n2= E(oldVineGraph)[con[2]]$CondName1 }
    if(progress == TRUE) message(n1," + ",n2," --> ", E(mst)[i]$name)
    E(mst)[i]$CondName2 = n1
    E(mst)[i]$CondName1 = n2
  }
  return(mst)
}

# Convert RVine object to R-vine array
asRVA= function(RVine)
{ n=length(RVine$Tree)+1
  con= list()
  names= V(RVine$Tree[[1]])$name
  conditionedSets=NULL
  correspondingCors=list()
  conditionedSets[[n-1]][[1]] = (E(RVine$Tree[[n-1]])$conditionedSet)
  for(k in 1:(n-2))
  { conditionedSets[[k]] = E(RVine$Tree[[k]])$conditionedSet
    correspondingCors[[k]] = as.list(E(RVine$Tree[[k]])$weight)
  }
  correspondingCors[[n-1]] = list()
  correspondingCors[[n-1]][[1]] = 0
  A=matrix(0,n,n); pc=matrix(NA,n,n)
  # A=matrix(-1,n,n); pc=matrix(NA,n,n)
  for(j in n:2)
  { k=n+1-j
    w=unlist(conditionedSets[[j-1]][[1]])[1]
    # w=conditionedSets[[j-1]][[1]][1]
    A[j,j]=w
    A[(j-1),j]=unlist(conditionedSets[[j-1]][[1]])[2]
    # A[(j-1),j]=conditionedSets[[j-1]][[1]][2]
    pc[(j-1),j]=correspondingCors[[j-1]][[1]][1]
    if(j==2)
    { A[(j-1),(j-1)]=unlist(conditionedSets[[j-1]][[1]])[2] }
    # { A[(j-1),(j-1)]=conditionedSets[[j-1]][[1]][2] }
    else
    { for(ii in (j-2):1)
      { i=n+1-ii
        for(jj in 1:length(conditionedSets[[ii]]))
        { cs=unlist(conditionedSets[[n-i+1]][[jj]])
          # cs=conditionedSets[[n-i+1]][[jj]]
          if(cs[1]==w)
          { A[ii,j]=cs[2]; break }
          else if(cs[2]==w)
          { A[ii,j]=cs[1]; break }
        }
        pc[ii,j]=correspondingCors[[ii]][[jj]][1]
        conditionedSets[[ii]][[jj]] = NULL
        correspondingCors[[ii]][[jj]] = NULL
      }
    }
  }
  # A=A+1
  # M is vinearray in format of VineCopula, Cor with M=A[n:1,n:1] etc
  M=A[n:1,n:1]; Cor=pc[n:1,n:1]
  return(list(VineA=A,pc=pc,Matrix=M,Cor=Cor))
}

# The following is in partialcorr.R and is only needed if
# this file is for standalone use.
# S = covariance or correlation matrix
# given = vector indices for the given or conditioning variables
# j,k = indices for the conditioned variables
# Output: partial correlation of variables j,k given indices in 'given'
#partcor= function(S,given,j,k)
#{ S11=S[given,given]
#  jk=c(j,k)
#  S12=S[given,jk]
#  S21=S[jk,given]
#  S22=S[jk,jk]
#  if(length(given)>1) { tem=solve(S11,S12); Om212=S21%*%tem }
#  else { tem=S12/S11; Om212=outer(S21,tem) }
#  om11=1-Om212[1,1]
#  om22=1-Om212[2,2]
#  om12=S[j,k]-Om212[1,2]
#  om12/sqrt(om11*om22)
#}

