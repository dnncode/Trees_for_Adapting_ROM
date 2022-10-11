# the function to split a tree node:

node.split = function(node.data, para, snapshot, POD, order){
  
  # input data structure
  # node.list[[1]] # lower bounds of search range for each parameter
  # node.list[[2]] # upper bounds of search range for each parameter
  # node.list[[3]] # data id on this node
  # node.list[[4]] # label 
  
  opt.value = Inf
  POD.L = NA
  POD.R = NA
  
  
  for (i.p in 1:ncol(para)){
    r = node.data[[2]][i.p]-node.data[[1]][i.p]
    #for (i.cut in seq(node.data[[1]][i.p],node.data[[2]][i.p],r/10)){
      for (i.cut in seq(node.data[[1]][i.p]+r/10,node.data[[2]][i.p]-r/10,r/10)){
      #print(c(i.p, i.cut))
      
      case = which( para[node.data[[3]],i.p]<= i.cut   )
      case.r = which( para[node.data[[3]],i.p] > i.cut   )
      
      if (length(case)>0){ 
        
        # obtain the data on the left nodes
        u.L = snapshot[node.data[[3]]][case]
        # local global POD
        u.combine = do.call(cbind,u.L)
        tmp.u = rSVD(u.combine, r=120)
        
        #POD.L = U[,1:K]
        POD.L = tmp.u[,1:order]
        
        # compute the loss on the left node
        d = array(0/0,dim=c(1,length(case)))
        j = 1
        for (i in node.data[[3]][case]){
          POD.j = POD[[i]]
          part = t(POD.L) %*% POD.j
          # svd.output.2 = fast.svd(part)
          # Sigma = svd.output.2$d
          
          tmp = t( part ) %*% part
          ev = eigen(tmp)
          Sigma = ev$values
          Sigma[Sigma>1] = 1
          
          if (length(case)==1){
            d[j] = 0
          }else{
            d[j] = sqrt( sum(acos(Sigma)^2, na.rm=TRUE) )
          }
          j = j+1
          #print(j)
        }
        loss.l = sum(d)
      }else{
        loss.l = 0
      }
      
      
      
      if (length(case.r)>0){
        # obtain the data on the right node
        u.R = snapshot[node.data[[3]]][case.r]
        # local global POD
        u.combine = do.call(cbind,u.R)
        tmp.u = rSVD(u.combine, r=120)
        
        #POD.R = U[,1:K]
        POD.R = tmp.u[,1:order]
        
        # compute the loss on the right node
        d = array(0/0,dim=c(1,length(case.r)))
        j = 1
        for (i in node.data[[3]][case.r]){
          POD.j = POD[[i]]
          part = t(POD.R) %*% POD.j
          #svd.output.2 = fast.svd(part)
          #Sigma = svd.output.2$d
          tmp = t( part ) %*% part
          ev = eigen(tmp)
          Sigma = ev$values
          Sigma[Sigma>1] = 1
          
          if (length(case.r)==1){
            d[j] = 0
          }else{
            d[j] = sqrt( sum(acos(Sigma)^2, na.rm=TRUE) )
          }
          j = j+1
          #print(j)
        }
        loss.r = sum(d)
      }else{
        loss.r = 0
      }
      
      
      # total loss
      loss = loss.l + loss.r
      
      # update
      if (loss < opt.value){
        opt.value = loss
        opt.p = i.p
        opt.cut = i.cut
        opt.POD.l = POD.L
        opt.POD.r = POD.R
        opt.loss.r = loss.r
        opt.loss.l = loss.l
        opt.case.l = case
        opt.case.r = case.r
      }
      
    }
  }
  
  
  output = list()
  node.list = list()
  
  # obtain left.daughter node information
  node.list[[1]] = node.data[[1]] # lower bounds of search range for each parameter
  tmp = node.data[[2]]
  tmp[opt.p] = opt.cut
  node.list[[2]] = tmp # upper bounds of search range for each parameter
  node.list[[3]] = opt.case.l # data id on this node
  node.list[[4]] = 0# label 
  node.list[[5]] = opt.POD.l
  output[[1]] = node.list
  
  # obtain right.daughter node information
  tmp = node.data[[1]]
  tmp[opt.p] = opt.cut
  node.list[[1]] = tmp # lower bounds of search range for each parameter
  node.list[[2]] = node.data[[2]] # upper bounds of search range for each parameter
  node.list[[3]] = opt.case.r # data id on this node
  node.list[[4]] = 0# label 
  node.list[[5]] = opt.POD.r
  output[[2]] = node.list
  
  if ( (length(opt.case.l)>0)+(length(opt.case.r)>0)<2){
    output[[3]]=1
  }else{
    output[[3]]=0
  }
  return(output)
}