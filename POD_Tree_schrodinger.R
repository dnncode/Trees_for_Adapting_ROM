# R code for:
# Regression Trees on Grassmann Manifold for Adapting Reduced-Order Models

# Code Written by: Xiao Liu (liuxiao314923@gmail.com)
# Updated on: Oct 2022

# Reference paper:
# Paper Authors: Xiao Liu (corresponding author) and Xinchao Liu, 
# Department of Industrial Engineering, University of Arkansas, Fayetteville, AR, 72701 
# arXiv: https://arxiv.org/abs/2206.11324

# -----------------------
# Description:
# 1) This R code is used to create the regression trees on Grassmann manifold for adapting ROM
# 2) Sample data are provided to run the code (i.e., Example II in the paper)
# 3) The code can be recycled and modified for other examples
# -----------------------

# ------------------------------------
# ------------------------------------
# Code:
# ------------------------------------
# ------------------------------------

rm(list=ls())
library(corpcor)
library(RSpectra)
library(MASS)
library(R.matlab)
wd = " " # user-defined work directory
setwd(wd)

# sub-routines are needed
source("rSVD.R") # randomized SVD
source("node.split.R") # code for splitting a tree node

# ------------------------------------
# ------------------------------------
# Load the Sample Data (Example II: the Schrodinger equation)
# ------------------------------------
# ------------------------------------
# number of conditions:
p = seq(0.05,0.5,0.01) # this is the parameter set used to simulate the data
n = length(p) # number of parameters

# get the data (this is the simulated training data)
# Note that: the data can be generated from Data_Generation.R. 
# Hence, the users need to run Data_Generation.R first
load("u.list.schrodinger.RData")
# Here:
# u.list is a list. The length of the list equals "n"
# each component of the u.list contains the snapshot data generated under one parameter setting


if (FALSE){ # This part obtains the local POD basis for each parameter setting
  K = 10
  POD.list = list()
  for (i in 1:n){
    svd.output = fast.svd(u.list[[i]])
    U = svd.output$u
    Sigma = diag(svd.output$d)
    V = svd.output$v
    POD.list[[i]] = U[,1:K] 
    print(i)
  }
}
load("POD20.list.schrodinger.RData") # The local POD bases are provided; Alternatively, the users may run the code above to generated the bases

# ----------
# Plot of the distance between two parameters and the distance on Grassmann manifold
# ----------
if (FALSE){
  ii = 1
  dist.p = dist.G = array()
  for (i in 1:(length(p)-1)){
    for (j in (i+1):length(p)){
      dist.p[ii] = sqrt( sum( (p[i]-p[j])^2 ) )
      POD.j = POD.list[[j]]
      POD.i = POD.list[[i]]
      part = t(POD.i) %*% POD.j
      tmp = t( part ) %*% part
      ev = eigen(tmp)
      Sigma = ev$values
      dist.G[ii] = sqrt( sum(acos(Sigma)^2, na.rm=TRUE) )
      ii = ii + 1
    }
  }
  
  par(mar=c(6,4,1,1))
  plot(dist.p, dist.G, col="blue",
       xlab="Euclidean distance between parameters",
       ylab="Riemannian distance between two subspaces")
  cor(dist.p, dist.G)
}

K = ncol(POD.list[[1]]) # number of POD mode
data = POD.list # a list of POD bases on the node;
u = u.list # a list of snapshot data
n.p = length(p) # parameters
n.data = n.p

# standardize p (this step depends on the problem)
# In our code, all parameters are standardized on [0,1]
p = (p-0.05)/(0.5-0.05)
p = matrix(p, ncol=1)

# prepare the training/testing data
set.seed(10)
case.train = c( seq(1, n, 4), n)
data.train = data[case.train]
data.test = data[-case.train]
u.train = u[case.train] 
u.test = u[-case.train]
p.train = matrix( p[case.train,], ncol=1)
p.test = matrix( p[-case.train,], ncol=1)
n.data.train = nrow(p.train)
n.data.test = nrow(p.test)

# get the global POD for training data only
if (FALSE){
  K = 10
  tmp = crossprod(do.call(cbind,u.train))
  ev = eigs_sym(tmp, K, which = "LM")
  ev.values = ev$values
  tmp.sv = sqrt(ev.values[1:K])
  tmp.v <- ev$vectors[,1:K]
  u.global = do.call(cbind,u.train) %*% tmp.v %*% diag(1/tmp.sv)
}
load("POD20.global.train.schrodinger.RData")
# The global POD bases are provided; Alternatively, the users may run the code above to generated the bases


# ------------------------------------
# ------------------------------------
# Growing the tree
# ------------------------------------
# ------------------------------------
iteration.tree.list=list()
repeatition = 1 # because rSVD is random, one may choose to run the algorithm for a few times
for (i.time in 1:repeatition){ # by default, we only run the algorithm once. 
  print(i.time)
  
  leaf.list = list() # a list of leaf nodes
  n.leaf = 1 # number of leaf nodes
  non.leaf.list = list()
  threshold = 5 # below which a tree node is no longer split
  
  # create the list for the root node
  node.list = list()
  node.list[[1]] = array(0,dim=c(1,ncol(p.train))) # lower bounds of search range for each parameter
  node.list[[2]] = array(1,dim=c(1,ncol(p.train))) # upper bounds of search range for each parameter
  node.list[[3]] = c(1:n.data.train) # data id on this node
  node.list[[4]] = 0 # label of this node. 0: non-terminal node; 1: terminal node
  node.list[[5]] = u.global # global POD on the root node (this is the global POD from training data only)
  
  non.leaf.list[[1]] = node.list
  
  flag = TRUE
  i.while = 0
  while (flag == TRUE){
    i.while = i.while + 1
    print(rep(i.while, 20))
    
    n.non.leaf = 0
    non.leaf.list.tmp = list()
    # tree node splitting process start
    for (ii in 1:length(non.leaf.list)){
      
      split.data = non.leaf.list[[ii]]
      output = node.split(node.data = split.data, para=p.train, 
                          snapshot=u.train, POD=data.train, order=K) # a user-defined node splitting function
      
      if (output[[3]]==0){
        left.daughter = output[[1]]
        right.daughter = output[[2]]
        
        n.left = length( left.daughter[[3]] )
        if (n.left <= threshold ){
          left.daughter[[4]] =1
          leaf.list[[n.leaf]] = left.daughter
          n.leaf = n.leaf + 1
        }else{
          n.non.leaf = n.non.leaf + 1
          non.leaf.list.tmp[[n.non.leaf]] = left.daughter
        }
        
        n.right = length( right.daughter[[3]] )
        if (n.right <= threshold ){
          right.daughter[[4]] =1
          leaf.list[[n.leaf]] = right.daughter
          n.leaf = n.leaf + 1
        }else{
          n.non.leaf = n.non.leaf + 1
          non.leaf.list.tmp[[n.non.leaf]] = right.daughter
        }
      }else{
        
        split.data[[4]] = 1
        leaf.list[[n.leaf]] = split.data
        n.leaf = n.leaf + 1
        
      }
      
      
    }
    
    non.leaf.list = non.leaf.list.tmp
    flag = (length(non.leaf.list)>0)
    
  }
  
  iteration.tree.list[[i.time]] = leaf.list
}

#load("schrodinger/tree5.rSVD120.POD20.schrodinger.RData")




# ------------------------------------
# ------------------------------------
# The following code is used for comparing different methods.
# ------------------------------------
# ------------------------------------

f.error = array(0/0, dim=c(nrow(p.test),length(iteration.tree.list))) 
f.error.inf = f.error.rel = f.error #different errors
# f.error: Frobenious error
# f.error.inf: errors based on L_infinity
# f.error.rel: relative error

# errors for the proposed method
for (i.time in 1:length(iteration.tree.list)){
  leaf.list = iteration.tree.list[[i.time]]
  
  n.leaf = length(leaf.list)
  for (ii in 1:nrow(p.test)){
    print(ii)
    # find the leaf
    for (i in 1:n.leaf){
      p.range.low = leaf.list[[i]][[1]]
      p.range.high = leaf.list[[i]][[2]]
      tmp1 = ((p.test[ii,1]-p.range.low[1]) * (p.test[ii,1]-p.range.high[1]) <=0)
      if (tmp1 == 1){i.select = i}
    }
    
    # get the accuracy
    POD.pred = leaf.list[[i.select]][[5]]
    u.low = POD.pred %*% ( ginv(POD.pred) %*% u.test[[ii]] )
    f.error[ii, i.time] = norm(u.test[[ii]]-u.low, type = "F")
    
    Sum = 0
    for (jj in 1:ncol(u.low)){
      Sum = Sum + sum( (u.test[[ii]][,jj]-u.low[,jj])^2 )
    }
    f.error.L2[ii, i.time] = sqrt(Sum)
    
    Sum = -Inf
    for (jj in 1:ncol(u.low)){
      tmp = max(  abs(u.test[[ii]][,jj]-u.low[,jj])  )
      if (tmp>Sum){Sum=tmp}
    }
    f.error.inf[ii, i.time] = Sum
    
    Sum = 0
    for (jj in 1:ncol(u.low)){
      tmp1 = sqrt(sum( (u.test[[ii]][,jj]-u.low[,jj])^2))
      tmp2 = sqrt(sum( (u.test[[ii]][,jj])^2))
      Sum = Sum + tmp1/tmp2
    }
    f.error.rel[ii, i.time] = norm(u.test[[ii]]-u.low, type = "F")/norm(u.test[[ii]], type = "F")
    
  }
}

# errors when the global POD basis is used
f.global = array(0/0,dim=c(n.data.test,1))
f.global.inf = f.global.rel = f.global #different errors
for (i in 1:n.data.test){
  u.low = u.global %*% ( ginv(u.global) %*% u.test[[i]] )
  f.global[i] = norm(u.test[[i]]-u.low, type = "F")
  
  Sum = -Inf
  for (jj in 1:ncol(u.low)){
    tmp = max(  abs(u.test[[i]][,jj]-u.low[,jj])  )
    if (tmp>Sum){Sum=tmp}
  }
  f.global.inf[i] = Sum
  
  Sum = 0
  for (jj in 1:ncol(u.low)){
    tmp1 = sqrt(sum( (u.test[[i]][,jj]-u.low[,jj])^2))
    tmp2 = sqrt(sum( (u.test[[i]][,jj])^2))
    Sum = Sum + tmp1/tmp2
  }
  f.global.rel[i] = norm(u.test[[i]]-u.low, type = "F")/norm(u.test[[i]], type = "F")
  
  print(i)
}

# ---------------------
# Plot the errors
# ---------------------
# Plot the Frobenious error
par(mar=c(6,4,1,1))
p.test.original = p.test * (0.5-0.05) + 0.05
plot(p.test.original, f.error,type="h",
     ylim=c(0,5),col="darkgreen",
     xlab=expression(alpha),ylab="error (Frobenius norm)",
     sub="(minimum samples in a leaf:5, rank of POD basis:20)")
lines(p.test.original, f.global,col="red",ylim=c(0,5),
      lwd=2,lty=2)
legend("topright", legend=c("Predicted POD basis by trees","Global POD"),
       col=c("darkgreen","red"),lty=c(1,2), bty = "n")

# Plot the error based on the L_infinity norm
par(mar=c(6,4,1,1))
p.test.original = p.test * (0.1-0.001) + 0.001
#plot(f.error/f.global)
plot(p.test.original, f.error.inf,type="h",
     ylim=c(0,0.15),col="darkgreen",
     xlab=expression(gamma),ylab="error (L-infinity norm)",
     sub="(minimum samples in a leaf:5, rank of POD basis:20)")
lines(p.test.original, f.global.inf,col="red",ylim=c(0,0.15),
      lwd=2,lty=2)
legend("topright", legend=c("Predicted POD basis by trees","Global POD"),
       col=c("darkgreen","red"),lty=c(1,2), bty = "n")

# Plot the relative error 
par(mar=c(6,4,1,1))
p.test.original = p.test * (0.1-0.001) + 0.001
#plot(f.error/f.global)
plot(p.test.original, f.error.rel*100,type="h",
     ylim=c(0,8),col="darkgreen",
     xlab=expression(gamma),ylab="error (in percentage)",
     sub="(minimum samples in a leaf:5, rank of POD basis:10)")
lines(p.test.original, f.global.rel*100,col="red",ylim=c(0,8),
      lwd=2,lty=2)
legend("topright", legend=c("Predicted POD basis by trees","Global POD"),
       col=c("darkgreen","red"),lty=c(1,2), bty = "n")




# ---------------
# ---------------
# Comparison between the proposed method, interpolation method and the use of global POD basis
# ---------------
# ---------------
# ---- sub-routine
# polynomial bases for polynomial interpolation
lbase = function(x,ii,seq){
  tmp = (x-seq)/(seq[ii] - seq)
  case = which( (tmp==Inf)+(tmp==-Inf) > 0 )
  tmp = tmp[-case]
  output =  prod(tmp)
  return(output)
}

# compute the Frobenious error for the interpolation method
f.error.int = array(0/0, dim=c(n.data.train, n.data.test))
stable.condition = array(0/0, dim=c(n.data.train, n.data.test))
for (i.ref in 1:n.data.train){
  
  POD.ref = as.matrix( data.train[[i.ref]] )
  
  # compute all Z.list matrices
  Z.list = list()
  for (i in 1:n.data.train){
    POD.i = as.matrix( data.train[[i]] )
    part1 = solve( t(POD.ref)%*%POD.i )
    part2 = POD.i %*% part1 - POD.ref
    
    tmp = fast.svd(part2)
    Sigma = tmp$d
    u = tmp$u
    v = tmp$v
    
    Z.list[[i]] = u %*% atan(diag(Sigma)) %*% t(v)
  }
  
  # interpolation and error calculation
  for (i.target in 1:n.data.test){
    Z.inter = array(0,dim=dim(Z.list[[1]]))
    
    if (TRUE){ #lagrangian
      for (i in 1:n.data.train){
        w = lbase(x=p.test[i.target], ii=i,seq=p.train)
        Z.inter = Z.inter + w * Z.list[[i]]
        #print(w)
        # rowSums(p.test[i.target,]-p.train)^2
      }
    }
    
    # interpolated POD
    svd.Z.inter = fast.svd(Z.inter)
    U.Z = svd.Z.inter$u
    Sigma.Z = diag(svd.Z.inter$d)
    V.Z = svd.Z.inter$v
    POD.inter = POD.ref %*% V.Z %*% diag(cos( diag(Sigma.Z) )) + U.Z %*% diag(sin( diag(Sigma.Z) ))
    if (svd.Z.inter$d[1]>pi/2){
      stable.condition[i.ref, i.target] = 1
    }else{
      stable.condition[i.ref, i.target] = 0
    }
    
    # error:
    tmp = POD.inter %*% ginv(POD.inter)
    u.low = tmp %*% u.test[[i.target]]
    f.error.int[i.ref, i.target] = norm(u.test[[i.target]]-u.low, type = "F")
    #print(c(i.ref,i.target))
  }
}


# plot the unstable locations
p.test.original = p.test * (0.5-0.05) + 0.05
p.train.original = p.train * (0.5-0.05) + 0.05
grid.x = rep( p.train.original, each=length(p.test.original))
grid.y = rep( p.test.original, length(p.train.original)) 
grid.unstable = array(0/0,dim=c(1,2))
for (i in 1:nrow(stable.condition)){
  for (j in 1:ncol(stable.condition)){
    if (stable.condition[i,j]==1){
      grid.unstable = rbind(grid.unstable,
                            matrix(c(p.train.original[i],p.test.original[j]),nrow=1))
    }
  }
}
grid.unstable = grid.unstable[-1,]

par(mar=c(4,4,1,1))
plot(grid.x, grid.y,pch=20,col="darkgrey",cex=0.3,
     xlab="alpha for the reference point", ylab="target alpha")
points(grid.unstable,col="red",pch=4)

# remove unstable interpolations and plot the accuracy
case = which(stable.condition==1)
f.error.int[case] = 0/0

f.error.int.max = apply(f.error.int, 2, max, na.rm=TRUE)
f.error.int.min = apply(f.error.int, 2, min, na.rm=TRUE)
case = which(f.error.int.max=="-Inf")
f.error.int.max[case] = 0/0
case = which(f.error.int.min=="Inf")
f.error.int.min[case] = 0/0

y.max = max( log(c(f.error.int.max,f.error, f.global)), na.rm=TRUE )*1.05
y.min = min( log(c(f.error.int.min,f.error, f.global)), na.rm=TRUE )*0.95
f.error.int.mean = colMeans(f.error.int,na.rm=TRUE)

p.test.original = p.test * (0.5-0.05) + 0.05
par(mar=c(4,4,1,1))
plot(p.test.original,log(f.error.int.max),ylim=c(y.min,y.max),
     col="darkgrey",pch="-",xlab=expression(alpha),ylab="error (log scale)")
points(p.test.original,log(f.error.int.min),ylim=c(y.min,y.max),pch="-",col="darkgrey")
lines(p.test.original,log(f.error.int.mean),ylim=c(y.min,y.max),
      col="blue",lty=2,lwd=2)
lines(p.test.original,log(f.global),ylim=c(y.min,y.max),col="darkgreen",
      lty=3,lwd=2)
lines(p.test.original,log(f.error),ylim=c(y.min,y.max),col="red",
      lwd=2,lty=1)
legend("topleft",legend=c("proposed method", "global POD", "interpolation (mean)","interpolation (worse/best cases)"),
       lty=c(1,3,2,0),pch=c(NA,NA,NA,"-"),col=c("red","darkgreen","blue","darkgrey"),
       lwd=c(2,2,2,NA),bty="n")

# ------------------------------------
# ------------------------------------
# Additional visualizations
# ------------------------------------
# ------------------------------------
# plot the solutions of the Schrodinger equation and the reconstructed solutions from the predicted POD bases
library(colorRamps)
par(mar=c(4,4,2,1))
par(mfcol=c(2,3))
for (ii in c(5,20,33)){
  image(t(u.test[[ii]]),col=blue2green2red(400),
        ylab = "distance", xlab="time",
        main = paste(expression(alpha),"=", p.test[ii,1]*(0.5-0.05) + 0.05,sep="")  )
  n.leaf = length(leaf.list)
  
  # find the leaf
  for (i in 1:n.leaf){
    p.range.low = leaf.list[[i]][[1]]
    p.range.high = leaf.list[[i]][[2]]
    tmp1 = ((p.test[ii,1]-p.range.low[1]) * (p.test[ii,1]-p.range.high[1]) <=0)
    if (tmp1 == 1){i.select = i}
  }
  POD.pred = leaf.list[[i.select]][[5]]
  u.low.tree = POD.pred %*% ( ginv(POD.pred) %*% u.test[[ii]] )
  image(t(u.low.tree),col=blue2green2red(400),
        ylab = "distance", xlab="time",
        main = paste(expression(alpha),"=", p.test[ii,1]*(0.5-0.05) + 0.05,sep="") )
}






