


# for laplace & gumbel distr
library(VGAM)
library(extraDistr)

############ generate data for a pair (x, y), where y = x + error
# n: number of data points

pairdatagen = function(n = NULL, x = NULL, error = T, mx = NULL, sdx = NULL, noninv = T, ninv.method = "gp", inv.method = "random",
                             se = 1, lo = -1, hi = 1, tau.p = NULL, noninv_split=T, error.distr = NULL){
  scale01 = function(x){
    (x-min(x))/diff(range(x))
  }
  
  if(is.null(mx)) mx = 0 # mean of x
  if(is.null(sdx)) sdx = 1 # sd of x
  if(is.null(x)) x = rnorm(n, mx, sdx) 

  ## random error
  if(error){
    se_normal <- runif(1, 0.5, 0.75)
    if(is.null(error.distr)){
      e = rnorm(n, 0, se_normal)
    }else if(error.distr == "normal"){
      e = rnorm(n, 0, se_normal)
    }else if(error.distr == "laplace"){
      e = rlaplace(n=n, location = 0, scale = 1)
    }else if(error.distr == "gumbel"){
      e = rgumbel(n=n, location = 0, scale = 1)
    }else if(error.distr == "t"){
      e = rt(n=n, df=5)
    }
  }else{
    e <- rep(0, n)
  }
  
  if(noninv){ #if we want non-invertible functions ONLY for nonlinear relations
    if(noninv_split){ # this will always be true when noninv_split==TRUE
      # gaussian process
      if(ninv.method == "gp"){
        l = runif(1, 5, 5.25)
        d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
        Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
        # generate 1 sample; returns 500*1 dim vector, since x is 500*1
        y = as.vector(mvtnorm::rmvnorm(n=1, sigma=Sigma_SE) + e)
      }
      
      if(ninv.method == "random"){
        ninv.method = sample(c("piecewise", "quadratic", "cubic", "x4"), 1)
      }
      
      if(ninv.method == "piecewise") {
        sign = sample(c(-1,1), 1)
        
        # tau.t = if(is.null(tau.p)){
        #   mean(sample(x, 10))
        # }else{
        #   quantile(x, probs = tau.p)
        # } 
        # tau.t = unname(quantile(x,probs=c(0.1, 0.125, 0.15, 0.2, 0.8, 0.85, 0.875, 0.9)))
        tau.t = unname(quantile(x,probs=c(0.15, 0.2, 0.25, 0.3, 0.7, 0.75, 0.8, 0.85)))
        x1 = as.matrix(x[which(x < tau.t)])
        x2 = as.matrix(x[which(x >= tau.t)])
        
        b12 = sample(tail(abs(runif(15, lo, hi))), 2)
        b1 = b12[1] * sign
        b2 = b12[2] * sign * -1/4
        a1 = runif(1, lo, hi)
        a2 = a1 + tau.t * (b1 - b2)
        
        y = e
        y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
        y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
        
      }
      
      
      if(ninv.method == "nl1") y = sample(c(-1,1), 1)*max(runif(2, lo, hi))*x^2 + sample(c(-1,1), 1)*min(runif(2, lo, hi))*x + e
      
      if(ninv.method == "nl2") y = sample(c(-1,1), 1)*cos(max(2, 2*runif(1, lo, hi))*x) + e
      
      if(ninv.method == "nl3") y = sample(c(-1,1), 1)*runif(1, lo, hi)*x^3 + 2*runif(1, lo, hi)*x^2 + e
      
      if(ninv.method == "nl4") y = sample(c(-1,1), 1)*tanh(x) + sample(c(-1,1), 1)*cos(2.5*x) + sample(c(-1,1), 1)*x^2 + e
      
      if(ninv.method == "cubic"){
        y = sample(c(-1,1), 1)*(x-1)^2*(runif(1, 0.75, 1.25)*x+3) + e
      }
      if(ninv.method == "quadratic") y= sample(c(-1,1), 1)*(x+runif(1, 0.1, 0.65))^2 + e
      
      if(ninv.method == "x4") y = sample(c(-1,1), 1)*(0.75*x^4-2*x^2-x^3) + e
      
      
    }else{
      sign = sample(c(-1,1), 1)
      if(inv.method == "random"){
        inv.method = sample(c("piecewise", "exp", "cubic", "sinh"), 1)
        # inv.method = sample(c("piecewise", "exp", "cubic", "log", "sinh", "tanh", "sigmoid", "mt1", "mt2", "mt3"), 1)
      }
      if(inv.method == "piecewise") {
        tau.t = mean(sample(x, 10))
        x1 = as.matrix(x[which(x < tau.t)])
        x2 = as.matrix(x[which(x >= tau.t)])
        
        b12 = sample(tail(abs(runif(15, lo, hi))), 2)
        b1 = b12[1] * sign
        b2 = b12[2] * sign
        a1 = runif(1, lo, hi)
        a2 = a1 + tau.t * (b1 - b2)
        
        y = e
        y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
        y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
      }else if(inv.method == "exp"){
        y = runif(1, 0.05, 0.35)*exp(x) + e
      }else if(inv.method == "cubic"){
        y = sign * x^3 * runif(1, 0.1, 0.5) + e
      }else if(inv.method == "sinh"){
        y = sign * sinh(x) + e
      }else if(inv.method == "sigmoid"){
        # y = runif(1, 0, 5) * runif(1, 0, 1)*x/(1+abs(runif(1, 0, 1))*x)
        # y = 1/(1+exp(x-runif(1,0,0.5)))^runif(1, 0.01, 5)
        # y = 1/(1+exp(x-runif(1,0,0.3))-runif(1,-1,1)) + e
        y = 1/(1+exp(x)) + e
      }else if (inv.method == "tanh"){
        y = sign * tanh(x) + e
      }else if(inv.method == "log"){
        x = scale01(x)+0.01
        y = log(x) + e
      }else if(inv.method == "cam_sigmoid"){
        a <- rexp(1, rate=4)+1
        b <- runif(1, 0.5, 2) * sample(c(-1,1), 1)
        c <- runif(1, -2, 2)
        y <- a * (b*(x+c))/(1+abs(b*(x+c))) + e
      }
    }
  }else{ # just linear vs nonlinear split in edges
    if(noninv){ # nonlinear edges
      sign = sample(c(-1,1), 1)
      if(ninv.method == "random"){
        # ninv.method = sample(c("piecewise", "nl1", "nl2", "nl3", "nl4", "piecewise", "exp", "cubic", "log", "sinh", "tanh", "sigmoid"), 1) 
        ninv.method = sample(c("piecewise", "exp", "cubic", "log", "sinh", "tanh", "sigmoid"), 1) 
      }
      if(ninv.method == "nl1") y = sample(c(-1,1), 1)*max(runif(2, lo, hi))*x^2 + sample(c(-1,1), 1)*min(runif(2, lo, hi))*x + e
      
      if(ninv.method == "nl2") y = sample(c(-1,1), 1)*cos(max(2, 2*runif(1, lo, hi))*x) + e
      
      if(ninv.method == "nl3") y = sample(c(-1,1), 1)*runif(1, lo, hi)*x^3 + sample(c(-1,1), 1)*3.5*runif(1, lo, hi)*x^2 + e
      
      if(ninv.method == "nl4") y = sample(c(-1,1), 1)*tanh(x) + sample(c(-1,1), 1)*cos(2.5*x) + sample(c(-1,1), 1)*x^2 + e
      
    
      if(ninv.method == "piecewise") {
        # tau.t = mean(sample(x, 10))
        tau.t <- sample(quantile(x, c(runif(1, 0.15, 0.35), runif(1, 0.7, 0.9))), 1)
        x1 = as.matrix(x[which(x < tau.t)])
        x2 = as.matrix(x[which(x >= tau.t)])
        
        b12 = sample(tail(abs(runif(15, lo, hi))), 2)
        b1 =  b12[1] * sign
        b2 = b12[2] * sign
        a1 = runif(1, lo, hi)
        a2 = a1 + tau.t * (b1 - b2)
        
        y = e
        y[which(x < tau.t)] = y[which(x < tau.t)] + a1 + b1 * x1
        y[which(x >= tau.t)] = y[which(x >= tau.t)] + a2 + b2 * x2
      }else if(ninv.method == "sl"){
        a = runif(1, lo, hi)
        b = max(runif(10, lo, hi))
        y = a + b*x + e
      }else if(ninv.method == "exp"){
        y = exp(x) + e
      }else if(ninv.method == "cubic"){
        y = x^3 + e
      }else if(ninv.method == "log"){
        x = scale01(x)+0.01
        y = log(x) + e
      }else if(ninv.method == "sinh"){
        y = sinh(x) + e
      }else if (ninv.method == "tanh"){
        y = tanh(x) + e
      }else if(ninv.method == "sigmoid"){
        y = 1/(1+exp(-x))
      }else if(ninv.method == "mt1"){
        y = ifelse(x < 0, -x^2, x^2) + e
      }
    }else{
      # simple linear regression
      a = runif(1, lo, hi)
      b = max(runif(10, lo, hi))
      y = a + b*x + e
    }
    
    
  }
  
  
  dat = NULL
  dat$x = scale(x) #scale(x)
  dat$y = scale(y) #scale(y)
  return(dat)
}


# plot(pairdatagen(n = 500, se = 0.5, noninv = T), cex = 0.5)
# plot(pairdatagen(n = 500, se = 0.5, noninv = F), cex = 0.5)






################ generate data for node Y=f(Xinv)+f(Xninv)


nodedatagen = function(n = NULL, Xinv = matrix(numeric(0),0,0), Xninv = matrix(numeric(0),0,0), 
                       nodeN.inv = NULL, nodeN.ninv = NULL, se = 1, ninv.method = "gp", 
                       inv.method = "random", noninv_split=noninv_split, error.distr = NULL, edge_fx_list = NULL){
  
  Xinv = as.matrix(Xinv); Xninv = as.matrix(Xninv)
  nodeN.inv = ncol(Xinv)
  nodeN.ninv = ncol(Xninv)
  # cat("\n invertible parent: ", colnames(Xinv))
  # cat("\n invertible parent fx: ", edge_fx_list[colnames(Xinv)])
  # print(c(nodeN.inv, nodeN.ninv))
  # cat("\n non-invertible parent: ", colnames(Xninv))
  # cat("\n non-invertible parent fx: ", edge_fx_list[edge_fx_list!="0"])
  # meant for linear edges
  if(nodeN.inv > 0) {
    Yinv = sapply(1:nodeN.inv, function(i) pairdatagen(n, x = Xinv[, i], error = F, noninv = F, ninv.method = ninv.method, 
                                                       inv.method = inv.method, se = se, noninv_split=noninv_split, error.distr = error.distr)$y)
  }else {
    Yinv = matrix(0, n, 1)
  }
  # meant for non-linear edges
  if(length(edge_fx_list)==0){ # no pre-specified edge fx
    edge_fx_list_ninv <- rep(ninv.method, nodeN.ninv)
    edge_fx_list_inv <- rep(inv.method, nodeN.ninv)
  }else{ # use pre-specified functions
    edge_fx_list_ninv <- edge_fx_list
    edge_fx_list_inv <- edge_fx_list
  }
  if(nodeN.ninv > 0) {
    Yninv = sapply(1:nodeN.ninv, function(i) pairdatagen(n, x = Xninv[, i], error = F, noninv = T, ninv.method = edge_fx_list_ninv[i], 
                                                         inv.method = edge_fx_list_inv[i], se = se, noninv_split=noninv_split, error.distr = error.distr)$y)
  }else{
    Yninv = matrix(0, n, 1)
  }
  
  # gaussian_se <- runif(1, min=0.75, max=1.25)
  gaussian_se <- runif(1, min=0.2, max=0.3)
  non_gaussian_se <- runif(1, min=0.5, max=0.75)
  if(is.null(error.distr)){
    e = rnorm(n, 0, gaussian_se)
  }else if(error.distr == "normal"){
    e = rnorm(n, 0, gaussian_se)
  }else if(error.distr == "laplace"){
    e = VGAM::rlaplace(n=n, location = 0, scale = non_gaussian_se)
  }else if(error.distr == "gumbel"){
    e = VGAM::rgumbel(n=n, location = 0, scale = non_gaussian_se)
  }else if(error.distr == "t"){
    e = extraDistr::rlst(n=n, df=5, mu=0, sigma = non_gaussian_se)
  }

  
  y = rowSums(Yinv) + rowSums(Yninv) + e
  
  dat = NULL
  dat$Xinv = Xinv
  dat$Xninv = Xninv
  dat$y = scale(y)
  return(dat)
}



DAGdatagen = function(n = NULL, dag.bn = NULL, nodeN = NULL, labels = NULL, dagSparsityProb = NULL, 
                      ninv.perct = 1, nonlinear.perct = 1, nonlinear.arcs = NULL, se = 1, 
                      ninv.method = "gp", inv.method = "random", error.distr=NULL, edge_fx_amat=NULL) {
  
  # n - data size
  library(igraph)
  library(bnlearn)
  
  if(is.null(labels)) labels = paste("X", 1:nodeN, sep = "")
  
  if(is.null(dag.bn)){
    if(is.null(dagSparsityProb)) dagSparsityProb = 2/(nodeN-1)
    dag.bn = bnlearn::random.graph(as.character(1:nodeN), prob = dagSparsityProb)
    nodes(dag.bn) = labels
  }else{
    nodeN = length(dag.bn$nodes)
    labels = names(dag.bn$nodes)
    # labels = paste("X", 1:nodeN, sep = "")
    # dag.bn <- rename.nodes(dag.bn, labels)
    # nodes(dag.bn) = as.character(1:nodeN)
  }


  arcs.ind = dag.bn$arcs #apply(dag.bn$arcs, 2, as.numeric)
  
  n.arcs = nrow(arcs.ind)
  # not used if non-invertible relations are pre-specified
  n.arcs.nonlinear = round(n.arcs *  nonlinear.perct) # non-linear edge percentage, default 100%
  if(ninv.perct > 0){
    noninv_split <- TRUE
  }else{
    noninv_split <- FALSE
  }

  # declaring which are non-linear edges
  if(is.null(nonlinear.arcs)){
    nonlinear.arcs = matrix(arcs.ind[sample(1:n.arcs, n.arcs.nonlinear),], ncol= 2, byrow = F)
  }else{
    nonlinear.arcs = matrix(nonlinear.arcs, ncol = 2, byrow = F)
  }
  
  dat = scale(sapply(1:nodeN, function(i)  rnorm(n, 0, 1)))
  colnames(dat) = nodes(dag.bn) = labels
  
  nodes.order = names(topological.sort(as.igraph(dag.bn)))#as.numeric(topo_sort(as.igraph(dag.bn)))
  
  
  for(i in nodes.order) {
    if(i %in% arcs.ind[, 2]) {
      if(i %in% nonlinear.arcs[, 2]) {
        pa.ninv = labels %in% nonlinear.arcs[nonlinear.arcs[, 2] %in% i, 1]
        pa.inv = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1][!arcs.ind[arcs.ind[, 2] %in% i, 1]%in%labels[pa.ninv]]
        if(!is.null(edge_fx_amat)){
          edge_fx_list <- edge_fx_amat[pa.ninv,i]
        }else{
          edge_fx_list <- NULL
        }
        dat[, labels %in% i] = nodedatagen(n = n, Xinv = dat[, pa.inv], Xninv = dat[, pa.ninv], se = se, 
                               ninv.method = ninv.method, inv.method = inv.method, noninv_split=noninv_split, error.distr=error.distr, edge_fx_list = edge_fx_list)$y
      }else{
        pa.inv = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1]
        dat[, labels %in% i] = nodedatagen(n = n, Xinv = dat[, pa.inv], se = se, 
                               ninv.method = ninv.method, inv.method = inv.method, noninv_split=noninv_split, error.distr=error.distr)$y
      }
    }
  }
  
  # ground truth. cpdag with fixed non-invertible edge 
  
  
  # wl = matrix(labels[ninv.arcs], ncol = 2)
  # 
  # dag.bn$learning$whitelist = wl
  # cpdag.bn = cpdag(dag.bn, wlbl = TRUE)
  # cpdagAmat = amat(cpdag.bn)
  # 
  # cpdagAmat.ninv = matrix(0, nodeN, nodeN)
  # colnames(cpdagAmat.ninv) = rownames(cpdagAmat.ninv) = labels
  # cpdagAmat.ninv[ninv.arcs] = 1 
  
  
  # return(list(DAGdata = scale(dat), dag.bn = dag.bn, cpdag.bn = cpdag.bn, 
              # cpdagAmat = cpdagAmat, cpdagAmat.ninv = cpdagAmat.ninv))
  return(list(DAGdata = scale(dat), nonlinear.arcs= nonlinear.arcs, dag.amat = amat(dag.bn)))
}



### MODIFIED DATA GENERATING FUNCTION
# topological sort of dag
# take nodes [1...k] to have linear parents, [k+1...N] to have nonlinear parents
# purpose: generate sorted data to test our alg vs CAM
# want to see if ordering of nodes affects estimated structure
DAGdatagen_sorted = function(n = NULL, dag.bn = NULL, nodeN = NULL, labels = NULL, dagSparsityProb = NULL, ninv.perct = 1,
                      nonlinear.perct = 1, nonlinear.arcs=NULL, nonlinear.nodes=NULL, se = 1, ninv.method = "random", inv.method = "random") {
  
  # n - data size
  library(igraph)
  library(bnlearn)
  
  if(is.null(dag.bn)){
    if(is.null(dagSparsityProb)) dagSparsityProb = 2/(nodeN-1)
    dag.bn = bnlearn::random.graph(as.character(1:nodeN), prob = dagSparsityProb)
    labels = names(dag.bn$nodes)
  }else{
    nodeN = length(dag.bn$nodes)
  }
  
  if(is.null(labels)) labels = paste("X", 1:nodeN, sep = "")
  
  
  arcs.ind = dag.bn$arcs #apply(dag.bn$arcs, 2, as.numeric)
  
  n.arcs = nrow(arcs.ind)
  # not used if non-invertible relations are pre-specified
  n.nodes.nonlinear = round(nodeN *  nonlinear.perct) # non-linear edge percentage, default 100%
  if(ninv.perct > 0){
    noninv_split <- TRUE
  }else{
    noninv_split <- FALSE
  }
  
  nodes.order = names(topo_sort(as.igraph(dag.bn)))#as.numeric(topo_sort(as.igraph(dag.bn)))
  nonlinear.nodes <- tail(nodes.order, n.nodes.nonlinear)
  # declaring which are non-invertible edges
  if(is.null(nonlinear.arcs)){
    # if no node specified to have nonlinear parents, then randomly sample arcs
    nonlinear.arcs = matrix(arcs.ind[arcs.ind[,2] %in% nonlinear.nodes,], ncol= 2, byrow = F)
  }else{
    # otherwise set the parents as nonlinear edges
    nonlinear.arcs = matrix(nonlinear.arcs, ncol = 2, byrow = F)
  }
  
  dat = scale(sapply(1:nodeN, function(i)  rnorm(n, 0, 1)))
  
  
  for(i in nodes.order) {
    if(i %in% arcs.ind[, 2]) {
      if(i %in% nonlinear.arcs[, 2]) {
        pa.ninv = labels %in% nonlinear.arcs[nonlinear.arcs[, 2] %in% i, 1]
        pa.inv = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1][!arcs.ind[arcs.ind[, 2] %in% i, 1]%in%labels[pa.ninv]]
        dat[, labels %in% i] = nodedatagen(n = n, Xinv = dat[, pa.inv], Xninv = dat[, pa.ninv], se = se, 
                                           ninv.method = ninv.method, inv.method = inv.method, noninv_split=noninv_split)$y
      }else{
        pa.inv = labels %in% arcs.ind[arcs.ind[, 2] %in% i, 1]
        dat[, labels %in% i] = nodedatagen(n = n, Xinv = dat[, pa.inv], se = se, 
                                           ninv.method = ninv.method, inv.method = inv.method, noninv_split=noninv_split)$y
      }
    }
  }
  
  # ground truth. cpdag with fixed non-invertible edge 
  
  colnames(dat) = nodes(dag.bn) = labels
  
  return(list(DAGdata = scale(dat), nonlinear.nodes=nonlinear.nodes, nonlinear.arcs= nonlinear.arcs, dag.amat = amat(dag.bn)))
  # return(list(DAGdata = dat, ninv.arcs= ninv.arcs, dag.amat = amat(dag.bn)))
}

