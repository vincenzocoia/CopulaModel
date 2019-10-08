# updated 2015.08.19 for CopulaModel
# function name change , 
# common functions with gausstrvineMST.R have been commented out

# Much of the code and documentation in this file is written by E Brechmann.

# This function is based on the algorithm in
#  Brechmann and Joe (2015), Truncation of vine copulas using fit indices
#     J Multivariate Analysis, v 138, 19-33
# A genetic algorithm is combined with minimum spanning trees.
# Inputs: 
#   cormat = correlation matrix
#   target = target fit value (near but less than 1; 1-alpha)
#   nneigh = number of nearest neighbors of minimum spanning tree
#   nbmst = number of best minimum spanning trees
#   selectcrit = selection criterion ("prop" for proportional or "rank")
#   selectfact = selection factor 
#    (total fraction of solutions retained after elitism and random selection)
#   elite = percentage for elitism (default 0.05)
#   method = fit index ("nfi", "cfi" or "ifi")
#       NFI=normed fit index, CFI=comparative fit index,
#       IFI=incremental fit index
#   n = sample size of data which the correlation matrix is computed from
#   iprint = print progress (Boolean)
# Output: 
#  $RVM with 
#  $RVM$VineA  vine array in CopulaModel format, variables [1:d,1:d]
#  $RVM$pc  partial correlations by tree [1:d,1:d]
#  $RVM$Matrix Rvine Matrix in VineCopula format, variables [d:1,d:1]
#  $RVM$Cor   partial correlations in VineCopula format [d:1,d:1]
#  $msts with
#  $msts[[1]]  (first tree)
#  ...
#  $msts[[d-1]]  (last tree where d=dimension)
#  $treeweight  vector of length d-1 with 
#                 sum_edge log(1-rho[edge]^2) for trees 1,...d-1
#  $trunclevel  same as inputted ntrunc
#  $truncval = sum_{1:trunclevel} treeweight/ sum_{1:(d-1)} treeweight

gausstrvine.galg=function(cormat,target=0.99,nneigh=10,nbmst=10,
  selectcrit="prop",selectfact=0.1,elite=0.05,method="nfi",n=0,iprint=T)
{ 
  cormat = as.matrix(cormat)
  # determine log determinant of the empirical correlation matrix
  true = log(det(cormat))
  # determine dimension and standardize column names
  d = dim(cormat)[1]
  colnames(cormat) = rownames(cormat) = paste("V",1:d,sep="")
  # calculate lambda_0 for CFI and IFI
  if(method %in% c("cfi","ifi"))
  { lam0 = -(n*true+d*(d-1))
    if(lam0 < 0) stop("Sample size is too small.")
  }
  # number of candidate models per tree
  ntrees0 = nneigh+nbmst
  # round to the upper bound on the number of neighbors
  nneigh = min((d-2)*(d-1)/2,nneigh)
  # actual number of candidate models per tree
  ntrees = nneigh+nbmst
  # initialize variables
  truncated = FALSE
  RVineTree = list()
  msts = list()
  
  # tree 1
  if(iprint) cat("Tree",1)
  # determine complete graph with edge weights according to the correlation matrix
  g = initFirstG(cormat)
  # determine candidate trees
  msts[[1]] = newT(g,nneigh,nbmst)
  msts[[d]] = list()
  # determine the weight of the candidate trees (1-truncated vines)
  msts[[d]]$cumweight = msts[[1]]$weight
  # determine the fitness of the candidate trees as: fit_index/(1-alpha)
  if(method == "nfi")
  { msts[[d]]$cumfitness = (msts[[d]]$cumweight/true)/target }
  else if(method == "cfi")
  { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1)+2-2*d)
    cfi = 1-pmax(lam,0)/pmax(lam,lam0,0)
    msts[[d]]$cumfitness = cfi/target
  }
  else if(method == "ifi")
  { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1))
    ifi = 1-lam/lam0
    msts[[d]]$cumfitness = ifi/target
  }
  msts[[d]]$cumlabel = msts[[1]]$label
  
  # number of candidate models to be selected
  nselect = ntrees
  if(iprint) cat(":",nselect,"candidate(s)\n")
  
  # check if truncation rule is satisfied
  if(any(msts[[d]]$cumfitness>=1))
  { # set truncation level to 1
    trunclevel = 1
    truncval = msts[[1]]$weight[1]/true
    truncated = TRUE
    # store best fit: the MST
    mstsOld = msts
    msts[[1]]$mst = list()
    msts[[1]]$mst[[1]] = mstsOld[[1]]$mst[[1]]
    msts[[1]]$weight = mstsOld[[1]]$weight[1]
    msts[[1]]$label = mstsOld[[1]]$label[1]
    
    msts[[d]] = list()
    msts[[d]]$cumweight = msts[[1]]$weight
    if(method == "nfi")
    { msts[[d]]$cumfitness = (msts[[d]]$cumweight/true)/target }
    else if(method == "cfi")
    { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1)+2-2*d)
      cfi = 1-pmax(lam,0)/pmax(lam,lam0,0)
      msts[[d]]$cumfitness = cfi/target
    }
    else if(method == "ifi")
    { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1))
      ifi = 1-lam/lam0
      msts[[d]]$cumfitness = ifi/target
    }
    msts[[d]]$cumlabel = msts[[1]]$label
    
    # set models to be selected to 1 (no other candidates needed after truncation)
    nselect = 1
    nbmst = 1
    nneigh = 0
  }
  
  # fit partial correlation vine tree 1 for each candidate model
  VineTree = list()
  for(j in 1:nselect) VineTree[[j]] = fitFirstT(msts[[1]]$mst[[j]])
  RVineTree[[1]] = list()
  for(j in 1:nselect) RVineTree[[1]][[j]] = VineTree[[j]]
  
  # tree 2:(d-1)
  for(i in 2:(d-1))
  { if(iprint) cat("Tree",i)
    g = list()
    mstlist = list()
    weightmat = matrix(NA,ntrees,nselect)
    labelmat = matrix(NA,ntrees,nselect)
    for(j in 1:nselect)
    { # determine line graph with partial correlations as edge weights
      g[[j]] = buildNextG(VineTree[[j]],cormat)
      # determine candidate trees
      mstlist[[j]] = newT(g[[j]],nneigh,nbmst)
      # determine weight of the i-truncated partial correlation vines
      weightmat[1:length(mstlist[[j]]$label),j] = msts[[d]]$cumweight[j]+mstlist[[j]]$weight
      labelmat[1:length(mstlist[[j]]$label),j] = paste(msts[[d]]$cumlabel[j],mstlist[[j]]$label,sep="_")
    }
    
    # determine the fitness of the candidate models as: fit_index/(1-alpha)
    if(method == "nfi")
    { fitness = (weightmat/true)/target }
    else if(method == "cfi")
    { lam = -(n*(true-weightmat)+d*(d-1)+i*(i+1-2*d))
      cfi = matrix(NA,ntrees,nselect)
      for(jj in 1:nselect) cfi[,jj] = 1-pmax(lam[,jj],0)/pmax(lam[,jj],lam0,0)
      fitness = cfi/target
    }
    else if(method == "ifi")
    { lam = -(n*(true-weightmat)+d*(d-1))
      ifi = matrix(NA,ntrees,nselect)
      for(jj in 1:nselect) ifi[,jj] = 1-lam[,jj]/lam0
      fitness = ifi/target
    }
    
    # check if truncation rule is satisfied (if not yet truncated)
    if(!truncated & any(fitness >= 1, na.rm=TRUE))
    { # set truncation level to i
      trunclevel = i
      truncval = max(fitness,na.rm=TRUE)*target
      truncated = TRUE
      # determine best fit
      selected = which.max(fitness)
      # set models to be selected to 1 (no other candidates needed after truncation)
      nselect = 1
      nbmst = 1
      nneigh = 0
    }
    else if(nselect == 1)
    { # if already truncated: select MST1 to finish tree sequence (doesn't matter which trees are selected)
      selected = 1
    }
    else
    { # (random) selection of candidate solutions
      if(selectcrit == "prop")
      { sfit = sum(fitness,na.rm=TRUE)
        fit_prob = c(fitness)/sfit
        fit_prob[is.na(fit_prob)] = 0    # probabilities very similar, since fitness similar
      }
      else if(selectcrit == "rank")
      { nnna = sum(!is.na(weightmat))
        fit_prob = 2*(rank(fitness,na.last="keep")-1)/(nnna-1)/nnna
        fit_prob[is.na(fit_prob)] = 0     # probabilities more distinct, since wider range of ranks
      }
      nselectOld = nselect
      # retain x% of candidates but at least ntrees0 many
      nselect = min(sum(fit_prob>0),max(ntrees0,ceiling(sum(fit_prob>0)*selectfact)))
    
      #elitism
      preselected = tail(order(fit_prob),ceiling(sum(fit_prob>0)*elite))
      fit_prob[preselected] = 0
      selected = sample(1:(nselectOld*ntrees), nselect-length(preselected), prob=fit_prob)
      selected = sort(c(preselected,selected)) #elitism
    }
  
    if(iprint) cat(":",nselect,"candidate(s)\n")
    mstsOld = msts
    VineTreeOld = VineTree
    RVineOldTree = RVineTree
  
    # store results
    msts[[i]] = list()
    msts[[i]]$mst = list()
    for(k in 1:i) msts[[k]]$mst = list()
    for(k in 1:i) msts[[k]]$weight = rep(NA,nselect)
    for(k in 1:i) msts[[k]]$label = rep(NA,nselect)
    for(k in 1:i) RVineTree[[k]] = list()
    msts[[d]]$cumweight = rep(NA,nselect)
    msts[[d]]$cumlabel = rep(NA,nselect)
  
    for(j in 1:nselect)
    { # column and row of the selected models
      COL = ceiling(selected[j]/ntrees)
      ROW = selected[j]-(COL-1)*ntrees
      for(k in 1:(i-1))
      { msts[[k]]$mst[[j]] = mstsOld[[k]]$mst[[COL]]
        msts[[k]]$weight[j] = mstsOld[[k]]$weight[[COL]]
        msts[[k]]$label[j] = mstsOld[[k]]$label[[COL]]
        RVineTree[[k]][[j]] = RVineOldTree[[k]][[COL]]
      }
      msts[[i]]$mst[[j]] = mstlist[[COL]]$mst[[ROW]]
      msts[[i]]$weight[j] = mstlist[[COL]]$weight[ROW]
      msts[[i]]$label[j] = mstlist[[COL]]$label[ROW]
      msts[[d]]$cumweight[j] = weightmat[selected[j]]
      msts[[d]]$cumlabel[j] = labelmat[selected[j]]
      # determine partial correlation vine tree
      VineTree[[j]] = fitT(msts[[i]]$mst[[j]],VineTreeOld[[COL]])
    }
  
    # determine the fitness of the candidate models as: fit_index/(1-alpha)
    if(method == "nfi")
    { msts[[d]]$cumfitness = (msts[[d]]$cumweight/true)/target }
    else if(method == "cfi")
    { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1)+i*(i+1-2*d))
      cfi = 1-pmax(lam,0)/pmax(lam,lam0,0)
      msts[[d]]$cumfitness = cfi/target
    }
    else if(method == "ifi")
    { lam = -(n*(true-msts[[d]]$cumweight)+d*(d-1))
      ifi = 1-lam/lam0
      msts[[d]]$cumfitness = ifi/target
    }
  
    # store partial correlation vine tree
    RVineTree[[i]] = list()
    for(j in 1:nselect) RVineTree[[i]][[j]] = VineTree[[j]]
  }

  # prepare output
  RVine = list(Tree = NULL)
  for(i in 1:(d-1)) RVine$Tree[[i]] = RVineTree[[i]][[1]]
  # transform list of trees to R-vine matrix
  #RVM = asRVM(RVine)
  RVA = asRVA(RVine)
  return(list(RVM=RVA,msts=msts,trunclevel=trunclevel,truncval=truncval,
    truncvine=msts[[d]]$cumlabel))
}

# initialize the first graph
# build complete graph with correlation coefficients as edge weights and with appropriate vine labels
# Input: cormat = correlation matrix
# Output: graph g
#initFirstG= function(cormat,iprint=F)
#{ g = graph.adjacency(cormat, mode="lower",weighted=TRUE,diag=FALSE)
#  E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep=",")
#  for(i in 1:ecount(g)) E(g)$conditionedSet[[i]] = get.edges(g,i)
#  if(iprint) print(g)
#  return(g)
#}

# determine best spanning trees for given graph
# Inputs:
#   g = graph
#   nneigh = number of nearest neighbors of minimum spanning tree
#   nbmst = number of best minimum spanning trees
# Output:
#   mstlist
newT= function(g,nneigh,nbmst)
{ mst_bmst = kMST(g,nbmst)
  nbmst1 = length(mst_bmst$label)
  # determine neighbors of the MST
  if(nneigh > 0)
  { mst_neigh = treeNeighbors(g,mst_bmst$mst[[1]],nneigh)
    nneigh1 = length(mst_neigh$label)
  }
  # prepare output as list
  mstlist = list()
  mstlist$mst = list()
  for(i in 1:nbmst1) mstlist$mst[[i]] = mst_bmst$mst[[i]]
  if(nneigh == 0)
  { mstlist$weight = mst_bmst$weight
    mstlist$label = mst_bmst$label
  }
  else
  { for(i in 1:nneigh1) mstlist$mst[[nbmst1+i]] = mst_neigh$mst[[i]]
    mstlist$weight = c(mst_bmst$weight,mst_neigh$weight)
    mstlist$label = c(mst_bmst$label,mst_neigh$label)
  }
  return(mstlist)
}

# Determine edges included and not included in the tree
# Inputs:
#   g = graph
#   mst = minimum spanning tree
#   nneigh = number of nearest neighbors of minimum spanning tree
# Output:
#  mst = mst_neigh = neighboring MST?
#  weight = weight,
#  label = label
treeNeighbors= function(g,mst,nneigh)
{ 
  incl = which((E(g)$name %in% E(mst)$name))
  not.incl = which(!(E(g)$name %in% E(mst)$name))
  # weights of the edges
  ord = order(E(g)$weight)
  # order edges that are not included by weight
  which.max.not.incl = rev(ord[which(!(ord %in% incl))])
  # round to the max. possible number of neighbors
  ntrees = min(length(which.max.not.incl),nneigh)
  
  if(ntrees > 0)
  { mst_neigh = list()
    weight = rep(NA,ntrees)
    for(i in 1:ntrees)
    { # determine neighbors by adding and removing edges
      mst_neigh[[i]] = addNode(g,mst,which.max.not.incl[i])
      # calculate new tree weight
      weight[i] = sum(log(1-E(mst_neigh[[i]])$weight^2))
    }
    label = paste("N",1:ntrees,sep="")
  }
  else
  { mst_neigh = NA
    weight = NA
    label = NA
  }
  return(list(mst=mst_neigh,weight=weight,label=label))
}

# Add edge to tree
# Inputs:
#   g = graph
#   mst = minimum spanning tree
#   edge.incl = edge to include
# Output:
#   a new MST with a replacement edge (of a deleted edge)
addNode= function(g,mst,edge.incl)
{ 
  incl = which((E(g)$name %in% E(mst)$name))
  pos = rep(0,ecount(g))
  pos[c(incl,edge.incl)] = 1
  mst1 = delete.edges(g, E(g)[!pos])
  
  # remove edge to form new tree (remove an edge of the circle formed by adding the new edge)
  circ = girth(mst1)$circle
  circ.graph = induced.subgraph(mst1, circ, impl="auto") #subgraph(mst1, circ)
  del.names = E(circ.graph)$name
  which.del = rep(1,ecount(mst1))
  which.del[which(E(mst1)$name == min(del.names[del.names!=E(g)$name[edge.incl]]))] = 0
  mst = delete.edges(mst1, E(mst1)[!which.del])
  return(mst)
}

# determine names of the vine tree edges (no actual fitting here) (by J. Dissmann)
# This function adds attribute $CondName1, $CondName2, $Copula.CondName
# but otherwise doesn't change the graph
#fitFirstT= function(mst,iprint=F)
#{ d = ecount(mst)
#  for(i in 1:d)
#  { a = get.edges(mst,i)
#    E(mst)[i]$CondName1 = V(mst)[a[1]]$name
#    E(mst)[i]$CondName2 = V(mst)[a[2]]$name
#  }
#  if(iprint) 
#  { print(E(mst)[0]$CondName1); print(E(mst)[1]$CondName1)
#    print(E(mst)[0]$CondName2); print(E(mst)[1]$CondName2)
#  }
#  return(mst)
#}

# determine line graph (original version by J. Dissmann)
# some condensation by H Joe with intersect() and setdiff()
# Inputs:
#   oldVineGraph = previous graph
#   rmat = correlation matrix
#   iprint = print flag
#buildNextG= function(oldVineGraph,rmat,iprint=F)
#{ EL= get.edgelist(oldVineGraph)
#  dd= ecount(oldVineGraph) # number of edges in previous tree
#  if(iprint) cat("\n***ecount=", dd,"\n")
#  g= graph.full(dd) # all posible edges
#  V(g)$name = E(oldVineGraph)$name # previous graph had dd edges
#  V(g)$conditionedSet = E(oldVineGraph)$conditionedSet
#  if(!is.null(E(oldVineGraph)$conditioningSet)) 
#  { V(g)$conditioningSet = E(oldVineGraph)$conditioningSet }
#  
#  for(i in 1:(ecount(g))) # g is full graph so it goes thru all pairs
#  { con= get.edge(g,i) # internal indices for edge
#    temp= get.edges(oldVineGraph,con) # edges of prev tree
#    # this is a way of pairing up nodes (which were prev edges)
#    temp2=unique(c(temp))
#    ok=(length(temp2)==3)
#    if(ok)
#    { same=intersect(temp[1,],temp[2,])
#      other1=setdiff(temp[1,],same)
#      other2=setdiff(temp[2,],same)
#      E(g)[i]$nodes = paste(as.character(c(other1,other2,same)),collapse=",")
#      name.node1= strsplit( V(g)[con[1]]$name,split=" *[,|] *")[[1]]
#      name.node2= strsplit( V(g)[con[2]]$name,split=" *[,|] *")[[1]]
#      schnitt=intersect(name.node1,name.node2)
#      symdiff=setdiff(union(name.node1,name.node2),schnitt)
#      # new edge
#      E(g)[i]$name = paste(paste(symdiff, collapse= ","),paste(schnitt, collapse= ","),sep= " | ")
#      
#      l1= c(unlist(V(g)[con[1]]$conditionedSet),unlist(V(g)[con[1]]$conditioningSet))
#      l2= c(unlist(V(g)[con[2]]$conditionedSet),unlist(V(g)[con[2]]$conditioningSet))
#      schnitt=intersect(l1,l2)
#      symdiff=setdiff(union(l1,l2),schnitt)
#      suppressWarnings({E(g)$conditionedSet[i] = list(symdiff)})
#      suppressWarnings({E(g)$conditioningSet[i] = list(schnitt)})
#      conditioned=symdiff  # size 2
#      conditioning=schnitt # given
#      needed= sort(c(conditioned,conditioning))
#      ind= which(needed %in% conditioned)
#      #pp=partcor(rmat,conditioning,conditioned[1],conditioned[2]) 
#      pcornew is version if this file is source'd
#      pp=pcornew(rmat,conditioning,conditioned[1],conditioned[2]) 
#      E(g)[i]$weight = pp
#    }
#    E(g)[i]$todel = !ok
#  }
#  
#  g= delete.edges(g, E(g)[E(g)$todel])
#  if(iprint) print(g)
#  return(g)
#}

# S = covariance or correlation matrix
# different interface where unused indices of S is OK.
# given = vector indices for the given or conditioning variables
# j,k = indices for the conditioned variables
# Output: partial correlation of variables j,k given indices in 'given'
#pcornew= function(S,given,j,k)
#{ S11=S[given,given]
#  jk=c(j,k)
#  S12=S[given,jk]
#  S21=S[jk,given]
#  S22=S[jk,jk]
#  if(length(given)>1) { tem=solve(S11,S12); Om212=S21%*%tem }
#  else { tem=S12/S11; Om212=outer(S21,tem) }
#  om11=S[j,j]-Om212[1,1]
#  om22=S[k,k]-Om212[2,2]
#  om12=S[j,k]-Om212[1,2]
#  om12/sqrt(om11*om22)
#}

# Fit tree
# Inputs:
#   mst = minimum spanning tree
#   oldVineGraph = previous graph
#   progress = flag
#fitT= function(mst,oldVineGraph,progress=FALSE)
#{ dd=ecount(mst)
#  for(i in 1:dd)
#  { con= get.edge(mst,i)
#    temp= get.edges(oldVineGraph,con)
#    same=intersect(temp[1,],temp[2,])
#    # n1, n2 are the variable names in the symmetric differnce
#    if(temp[1,1]==same)
#    { n1= E(oldVineGraph)[con[1]]$CondName2 }
#    else
#    { n1= E(oldVineGraph)[con[1]]$CondName1 }
#    if(temp[2,1]==same)
#    { n2= E(oldVineGraph)[con[2]]$CondName2 }
#    else
#    { n2= E(oldVineGraph)[con[2]]$CondName1 }
#    if(progress == TRUE) message(n1," + ",n2," --> ", E(mst)[i]$name)
#    E(mst)[i]$CondName2 = n1
#    E(mst)[i]$CondName1 = n2
#  }
#  return(mst)
#}


# transform list of R-vine trees to R-vine matrix (by J. Dissmann)
# also return is the vine array A with reversed column/row order
#  and the matrix of partial correlations in the indexing of A
#  compared with the RvineMatrix.
#asRVA = function(RVine)
#{ n = length(RVine$Tree)+1
#  con = list()
#  names = V(RVine$Tree[[1]])$name
#  conditionedSets = NULL
#  corresppondingCors = list()
#  conditionedSets[[n-1]][[1]] = (E(RVine$Tree[[n-1]])$conditionedSet)
#  for(k in 1:(n-2))
#  { conditionedSets[[k]] = E(RVine$Tree[[k]])$conditionedSet
#    corresppondingCors[[k]] = as.list(E(RVine$Tree[[k]])$weight)
#  }
#  corresppondingCors[[n-1]] = list()
#  corresppondingCors[[n-1]][[1]] = 0
#  
#  Cor = matrix(NA,n,n)
#  M = matrix(NA,n,n)
#  for(k in 1:(n-1))
#  { w = unlist(conditionedSets[[n-k]][[1]])[1]
#    M[k,k] = w
#    M[(k+1),k] = unlist(conditionedSets[[n-k]][[1]])[2]
#    Cor[(k+1),k] = corresppondingCors[[n-k]][[1]][1]
#    if(k == (n-1))
#    { M[(k+1),(k+1)] = unlist(conditionedSets[[n-k]][[1]])[2] }
#    else
#    { for(i in (k+2):n)
#      { for(j in 1:length(conditionedSets[[n-i+1]]))
#        { cs = unlist(conditionedSets[[n-i+1]][[j]])
#          if(cs[1] == w)
#          { M[i,k] = cs[2]
#            break
#          }
#          else if(cs[2] == w)
#          { M[i,k] = cs[1]
#            break
#          }
#        }
#        Cor[i,k] = corresppondingCors[[n-i+1]][[j]][1]
#        conditionedSets[[n-i+1]][[j]] = NULL
#        corresppondingCors[[n-i+1]][[j]] = NULL
#      }
#    }
#  }
#  
#  M = M
#  M[is.na(M)] = 0
#  # HJ addition for compatibility with CopulaModel package
#  A=M[n:1,n:1]; pc=Cor[n:1,n:1]
#  return(list(Matrix=M,Cor=Cor,VineA=A,pc=pc))
#}

# algorithm by N. KATOH, T. IBARAKI AND H. MINES: "AN ALGORITHM FOR FINDING K MINIMUM SPANNING TREES"
# (similar version in H. N. GABOW: "TWO ALGORITHMS FOR GENERATING WEIGHTED SPANNING TREES IN ORDER")
# Inputs:
#   g = graph
#   ntrees = number of best spanning trees
# Output:
#    mst = list of msts
#    weight= list of weights of the MSTs
#    label= labels for the MSTs
kMST= function(g,ntrees)
{ weight_vec = log(1-E(g)$weight^2)
  P = list()
  Q = list()
  msts = list()
  msts_weight = rep(NA,ntrees)
  mst = minimum.spanning.tree(g, weights=weight_vec)
  w_mst = sum(log(1-E(mst)$weight^2))
  msts[[1]] = mst
  msts_weight[1] = w_mst
  
  if(ntrees > 1)
  { mst1 = kMST_EX(g,mst,IN=NULL,OUT=NULL)
    mst1.min = which.min(mst1$weight)
    if(length(mst1.min) == 0)
    { ntrees = 1 }
    else
    { P[[1]] = list(oldweight=w_mst, newweight=w_mst+mst1$weight[mst1.min], del=mst1$del[mst1.min], incl=mst1$incl[mst1.min], oldmst=mst1$oldmst, newmst=mst1$newmst[[mst1.min]], IN=NULL, OUT=NULL)
      Q[[1]] = mst1
      for(j in 2:ntrees)
      { Pold = P
        Qold = Q
        P = list()
        Q = list()
        P_weights = rep(NA,j-1)
        for(i in 1:(j-1)) P_weights[i] = Pold[[i]]$newweight
        if(all(P_weights == 1e300))
        { ntrees = j-1
          break
        }
        whichmin = which.min(P_weights)
        P_min = Pold[[whichmin]]
        msts[[j]] = P_min$newmst
        msts_weight[j] = P_min$newweight
        Qupdate = Qold[[whichmin]]
        update_min = which.min(Qupdate$weight)
        #Qupdate$newmst[[update_min]] = NULL #PRODUCES ERRORS!!!
        Qupdate$weight[update_min] = NA
        #Qupdate$incl[update_min] = NA
        #Qupdate$del[update_min] = NA
        Q[[whichmin]] = Qupdate
        
        #Q[[j]][[j]] = kMST_EX(g,P_min$newmst,IN=c(P_min$IN,P_min$incl),OUT=P_min$OUT)
        Q[[j]] = kMST_EX(g,P_min$newmst,IN=P_min$IN, OUT=c(P_min$OUT,P_min$del)) #GABOW
        for(i in (1:(j-1))[-whichmin]) Q[[i]] = Qold[[i]]
        if(all(is.na(Q[[whichmin]]$weight)))
        { P[[whichmin]] = list(newweight = 1e300) }
        else
        { new_min = which.min(Q[[whichmin]]$weight)
          P[[whichmin]] = list(oldweight=P_min$oldweight, newweight=P_min$oldweight+Q[[whichmin]]$weight[new_min],
          del=Q[[whichmin]]$del[new_min], incl=Q[[whichmin]]$incl[new_min],
          #oldmst=Q[[whichmin]]$oldmst, newmst=Q[[whichmin]]$newmst[[new_min]], IN=P_min$IN, OUT=c(P_min$OUT,P_min$incl))
          oldmst=Q[[whichmin]]$oldmst, newmst=Q[[whichmin]]$newmst[[new_min]], IN=c(P_min$IN,P_min$del), OUT=P_min$OUT) #GABOW
        }
        
        if(all(is.na(Q[[j]]$weight)))
        { P[[j]] = list(newweight = 1e300) }
        else
        { j_min = which.min(Q[[j]]$weight)
          P[[j]] = list(oldweight=P_min$newweight, newweight=P_min$newweight+Q[[j]]$weight[j_min],
          del=Q[[j]]$del[j_min], incl=Q[[j]]$incl[j_min],
          #oldmst=Q[[j]]$oldmst, newmst=Q[[j]]$newmst[[j_min]], IN=c(P_min$IN,P_min$incl), OUT=P_min$OUT)
          oldmst=Q[[j]]$oldmst, newmst=Q[[j]]$newmst[[j_min]], IN=P_min$IN, OUT=c(P_min$OUT,P_min$del)) #GABOW
        }
        for(i in (1:(j-1))[-whichmin]) P[[i]] = Pold[[i]]
      }
    }
  }
  
  label = paste("MST",1:ntrees,sep="")
  msts_weight = msts_weight[1:ntrees]
  return(list(mst=msts, weight=msts_weight, label=label))
}

# Inputs:
#   g = graph
#   mst = minimum spanning tree
#   IN = edges to remain in the new trees
#   OUT = edges to remain out of the new trees
# Output:
#   oldmst = old MST,
#   newmst = list of new MSTs
#   weight = list of weights of the MSTs
#   incl = included edges (w.r.t. g)
#   del = deleted edges (w.r.t. g)
kMST_EX= function(g,mst,IN,OUT)
{ ec = ecount(mst)
  ec_g = ecount(g)
  incl = which((E(g)$name %in% E(mst)$name))
  elig = rep(1,ec)
  elig[incl %in% IN] = 0
  L = order(log(1-E(g)$weight^2))
  L.not.incl = L[!(L %in% c(incl,OUT))]
  
  tex.weight = rep(NA,ec)
  tex.incl = rep(NA,ec)
  tex.del = rep(NA,ec)
  tex.mst = list()
  
  for(i in 1:length(L.not.incl))
  { for(j in 1:ec)
    { if(elig[j] == 1)
      { which.del = rep(1,ec_g)
        which.del[incl[-j]] = 0
        which.del[L.not.incl[i]] = 0
        mst1 = delete.edges(g, E(g)[which(which.del==1)])
        if(is.connected(mst1))
        { tex.weight[j] = log(1-E(g)$weight^2)[L.not.incl[i]]-log(1-E(g)$weight^2)[incl[j]]
          tex.incl[j] = L.not.incl[i]
          tex.del[j] = incl[j]
          tex.mst[[j]] = mst1
          elig[j] = 0
        }
      }
    }
    if(all(elig == 0)) break
  }
  
  #tex.min = which.min(tex.weight)
  #return(list(oldmst=mst, newmst=tex.mst[[tex.min]], weight=tex.weight[tex.min], incl=tex.incl[tex.min], del=tex.del[tex.min]))
  return(list(oldmst=mst, newmst=tex.mst, weight=tex.weight, incl=tex.incl, del=tex.del))
}
