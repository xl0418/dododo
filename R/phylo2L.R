library(ape)

phylo2L = function(emdata){
  # compute the relative branching times 
  brt=branching.times(emdata)
  if(min(brt)<0){
    brt = brt+abs(min(brt))
  }
  # number of species including extinct species.
  num.species=emdata$Nnode+1
  brt_preL = c(brt[emdata$edge[,1]-length(emdata$tip.label)])
  # check if the relative branching times are equal to the real branching times.
  # if not correct it to the real branching times.
  if(min(brt_preL) == 0){
  correction = max(emdata$edge.length[which(brt_preL==0)])
  brt_preL = brt_preL+correction
  }
  # preliminary L table
  pre.Ltable = cbind(brt_preL,emdata$edge,emdata$edge.length,brt_preL-emdata$edge.length)
  # identify the extant species and the extinct species
  extantspecies.index = pre.Ltable[which(pre.Ltable[,5]<=1e-10),3]
  tipsindex = c(1:num.species)
  extinct.index3 = subset(tipsindex,!(tipsindex %in% extantspecies.index))
  # assigen the extinct species with extinct times; the extant species with -1 and
  # the internal nodes with 0.
  eeindicator = matrix(0,length(emdata$edge.length),1)
  eeindicator[match(extantspecies.index,pre.Ltable[,3])]=-1
  ext.pos = match(extinct.index3,pre.Ltable[,3])
  eeindicator[ext.pos]= pre.Ltable[ext.pos,5]
  pre.Ltable=cbind(pre.Ltable,eeindicator)
  
  sort.L = pre.Ltable[order(pre.Ltable[,1],decreasing = TRUE),]
  nodesindex = unique(emdata$edge[,1])
  L = sort.L
  realL=NULL
  do=0
  while(do == 0){
  j = which.min(L[,3])
  daughter = L[j,3]
  parent = L[j,2]
  if(parent %in% nodesindex){
  L[which(L[,2]==parent),2] = daughter
  if(length(which(L[,3]==parent))==0){
    realL = rbind(realL,L[j,],row.names = NULL)
    L = L[-j,,drop=FALSE]
  }else{
  L[which(L[,3]==parent),6] = L[j,6]
  L[which(L[,3]==parent),3] = daughter
  L = L[-j,,drop=FALSE]
  }
  }else{
    realL = rbind(realL,L[j,],row.names = NULL)
    L = L[-j,,drop=FALSE]
  }
  
  if(nrow(L)==0){
    do = 1
  }
  }
  realL = realL[order(realL[,1],decreasing = T),]
  L = realL[,c(1,2,3,6)]
  
  daughter.index = L[,3]
  daughter.realindex = c(1:nrow(L))
  parent.index = L[,2]
  parent.realindex = match(parent.index, daughter.index)
  
  L[,2]=parent.realindex
  L[,3]=daughter.realindex
  L[1,2] = 0
  L[1,3] = -1
  L[2,2] = -1
  for(i in c(2:nrow(L))){
    if(L[i-1,3]<0){
      mrows = which(L[,2]==abs(L[i-1,3]))
      L[mrows,2] = L[i-1,3]
      L[mrows,3] = -1* L[mrows,3]
    }
  }
  dimnames(L) = NULL
  return(L)
}

