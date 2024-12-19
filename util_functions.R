


# loads rda file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# get JI
# percentage of correct edges among the union of edges in the estimated & true graphs
get_ji <- function(true_cpdag, temp_cpdag){
  true_amat <- amat(true_cpdag)
  temp_amat <- amat(temp_cpdag)
  # get true pos and remove duplicated edges (undirected edges are counted twice)
  tp <- bnlearn::compare(true_cpdag, temp_cpdag, arcs=TRUE)$tp
  tp <- nrow(matrix(tp[!duplicated(t(apply(tp, 1, sort))),], ncol=2))
  
  
  # get edge set size of true restricted cpdag & learned cpdag
  true_total_edge <- nrow(directed.arcs(true_cpdag)) + nrow(undirected.arcs(true_cpdag))/2
  temp_total_edge <- nrow(directed.arcs(temp_cpdag)) + nrow(undirected.arcs(temp_cpdag))/2
  ji <- tp/(true_total_edge+temp_total_edge-tp)
  
  return(ji)
}

# get F1 score
get_f1 <- function(true_cpdag, temp_cpdag_list){
  res <- c()
  for(i in c(1:length(temp_cpdag_list))){
    temp_cpdag <- temp_cpdag_list[[i]]
    # extract adjacency matrices
    if("matrix" %in% class(true_cpdag)){
      true_amat <- true_cpdag
    }else{
      true_amat <- amat(true_cpdag)
    }
    if("matrix" %in% class(temp_cpdag)){
      temp_amat <- temp_cpdag
    }else{
      temp_amat <- amat(temp_cpdag)
    }
    
    # get shd score
    # shd_score <- bnlearn::shd(true=true_cpdag, learned=temp_cpdag, wlbl = FALSE)
    
    # get true pos and remove duplicated edges (undirected edges are counted twice; should be one edge only)
    tp <- bnlearn::compare(true_cpdag, temp_cpdag, arcs=TRUE)$tp
    tp <- sum(!duplicated(t(apply(tp, 1, sort))))
    
    
    # get edge set size of true restricted cpdag & learned cpdag
    true_total_edge <- nrow(directed.arcs(true_cpdag)) + nrow(undirected.arcs(true_cpdag))/2
    # number of undirected edges in graph
    n_undir <- nrow(undirected.arcs(temp_cpdag))/2
    temp_total_edge <- nrow(directed.arcs(temp_cpdag)) + n_undir
    
    precision <- tp/temp_total_edge
    recall <- tp/true_total_edge
    f1 <- 2*precision*recall/(precision+recall)
    
    # get incorrectly oriented edges; check w/ skeleton
    tp_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag), 
                                    bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$tp
    tp_skeleton <- matrix(tp_skeleton[!duplicated(t(apply(tp_skeleton, 1, sort))),], ncol=2)
    
    # number of incorrectly oriented edges
    r <- 0
    # number of incorrectly oriented edges that are undirected (supposedly should be directed)
    r_undir <- 0
    if(nrow(tp_skeleton) > 0){
      for(i in 1:nrow(tp_skeleton)){
        n1 <- tp_skeleton[i,1]
        n2 <- tp_skeleton[i,2]
        # need edge direction to be same
        if(true_amat[n1,n2]!=temp_amat[n1,n2] || true_amat[n2,n1]!=temp_amat[n2,n1]){
          r <- r+1
          # check if wrong dir is b/c undirected edge
          if(temp_amat[n1,n2] + temp_amat[n2,n1]==2){
            r_undir <- r_undir + 1
          }
        }
      }
    }
    fp <- temp_total_edge-tp-r
    fn <- true_total_edge-tp-r
    
    # another way to calculate f1
    f1_2 <- 2*tp/(2*tp+fp+fn+2*r)
    
    # calculate SHD manually; fp+fn+wrong_dir
    shd_score <- fp+fn+r
    
    res <- rbind(res, c(shd_score, f1, tp, fp, fn, r, r_undir, n_undir))
    colnames(res) <- c("shd","f1", "tp", "fp", "fn", "r", "r_undir", "n_undir")
  }
  return(res)
}


get_wrong_dir <- function(true_cpdag, temp_cpdag){
  true_amat <- amat(true_cpdag)
  temp_amat <- amat(temp_cpdag)
  # get incorrectly oriented edges
  tp_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag), 
                                  bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$tp
  tp_skeleton <- matrix(tp_skeleton[!duplicated(t(apply(tp_skeleton, 1, sort))),], ncol=2)
  
  r_edges <- c()
  r <- 0
  num_temp_undir <- 0
  num_temp_undir2 <- 0
  for(i in 1:nrow(tp_skeleton)){
    if(nrow(tp_skeleton)==0){
      break
    }
    n1 <- tp_skeleton[i,1]
    n2 <- tp_skeleton[i,2]
    same_orientation <- (true_amat[n1,n2]==temp_amat[n1,n2]) + (true_amat[n2,n1]==temp_amat[n2,n1])
    if(same_orientation != 2){
      r <- r+1
      if(temp_amat[n1,n2]==temp_amat[n2,n1]){ # check if edge in learned DAG is undirected
        num_temp_undir2 <- num_temp_undir2+1
      }
    }
    if(true_amat[n1,n2]!=true_amat[n2,n1]){
      num_temp_undir <- num_temp_undir+1
    }
  }
  
  # get wrong dir, and directed edges
  return(c(r, r-num_temp_undir2, num_temp_undir2))
}

  # returns: percentage of undirected edges captured, shd, ji
get_undir_pct <- function(true_amat, cpdag_amat, temp_amat){
  # number of undirected edges undetected
  N.det.edge = sum((abs(temp_amat - cpdag_amat) + abs(true_amat - cpdag_amat))==2)
  # number of undirected edges total
  N.undir.edge = nrow(matrix(which((cpdag_amat == t(cpdag_amat) & (cpdag_amat != 0)), arr.ind = TRUE), ncol = 2))/2
  # percentage of undirected edges from CPDAG captured by GAM
  return(N.det.edge/N.undir.edge)
}


# extracts results from the GAM edge orientation step
analyze_cpdag <- function(temp.cpdag, true.cpdag){
  # temp.cpdag <- as.bn(graph_from_adjacency_matrix(temp$amat))
  # gam_cpdag <- temp$gam_cpdag
  # if(!is.null(gam_cpdag)){
  #   temp.cpdag$learning$whitelist <- gam_cpdag$learning$whitelist
  # }
  
  # if(!is.null(temp$sigpair)){
  #   # change_name <- function(x){
  #   #   return(paste("X",x,sep=""))
  #   # }
  #   # temp.cpdag$learning$whitelist <- matrix(temp$sigpair, ncol = 2, dimnames = list(NULL, c("from", "to")))
  #   # exclude_edges <- apply(temp$sigpair[,c(2,1)], 2, change_name)
  #   temp.cpdag$learning$whitelist <- rbind(temp.cpdag$learning$whitelist, apply(matrix(temp$sigpair[,c(2,1)], ncol = 2, dimnames = list(NULL, c("from", "to"))), 2, change_name))
  # }
  # edgeorient_time <- temp$total_time
  shd_cpdag <- bnlearn::shd(true=true.cpdag, learned=temp.cpdag, wlbl=TRUE)
  f1_cpdag <- get_f1(true.cpdag, temp.cpdag)
  ji_cpdag <- get_ji(true.cpdag, temp.cpdag)
  noninv_fn <- count_fn(true.cpdag, temp.cpdag)
  # return(c(shd_cpdag, f1_cpdag, ji_cpdag, noninv_fn, temp$gam_edge_time, temp$edgeorient_time))
  return(c(shd_cpdag, f1_cpdag, ji_cpdag, noninv_fn))
}


# get number of missing edges that are non-invertible/non-linear
count_fn <- function(true_cpdag, temp_cpdag){
  true_amat <- amat(true_cpdag)
  temp_amat <- amat(temp_cpdag)
  
  fn_skeleton <- bnlearn::compare(bnlearn:::skeleton(true_cpdag), 
                                  bnlearn:::skeleton(temp_cpdag), arcs=TRUE)$fn
  if(!is.null(true_cpdag$learning$whitelist)){
    non_inv_fn <- nrow(inner_join(as.data.frame(true_cpdag$learning$whitelist), as.data.frame(fn_skeleton), by = c("from", "to")))
  }else{
    non_inv_fn <- 0
  }
  return(non_inv_fn)
}


# get true dag of network
# input: network name
# output: bnlearn obj DAG
get_true_dag <- function(network='child', fixed_edges=TRUE){
  if(network=="1"){
    # EXAMPLE 1
    true_amat <- rbind(c(0,1,1,0), c(0,0,0,0), c(0,0,0,1), c(0,0,0,0))
    colnames(true_amat) <- c("x1","x2","x3","x4")
    rownames(true_amat) <- c("x1","x2","x3","x4")
    ninv.arcs <- rbind(c("x1", "x2"), c("x1", "x3"))
  }else if(network=="2"){
    # EXAMPLE 2
    # same network as network 3, but with node 6 removed
    true_amat <- rbind(c(0,1,0,0,0), c(0,0,1,1,0),c(0,0,0,0,0), 
                       c(0,0,0,0,0), c(0,0,0,1,0))
    colnames(true_amat) <- c("x1","x2","x3","x4","x5")
    rownames(true_amat) <- c("x1","x2","x3","x4","x5")
    # ninv.arcs <- rbind(c("x2", "x3"))
    ninv.arcs <- c("x2", "x4")
  }else if(network=="3"){
    # EXAMPLE 3
    # true dag: x1->x2->x3, x2->x4<-x5, x4->x6
    true_amat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0), 
                       c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
    colnames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
    rownames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
    ninv.arcs <- rbind(c("x2", "x4"))
  }else if(network=="asia"){
    asia <- loadRData("asia.rda")
    true_amat <- amat(asia)
    ninv.arcs <- rbind(c("tub", "either"),c("smoke","lung"))
  }else if(network=="sachs"){
    sachs <- loadRData("/Users/StellaHuang/Documents/UCLA/Research/GAM/sachs.rda")
    true_amat <- amat(sachs)
    ninv.arcs <- rbind(c("PKC", "Mek"),c("PIP3","PIP2"), c("PKA","P38"), c("PKA","Erk"))
  }else if(network=="child"){
    child <- loadRData("child.rda")
    true_amat <- amat(child)
    ninv.arcs <- rbind(c("BirthAsphyxia", "Disease"),c("LungParench","ChestXray"), 
                       c("CardiacMixing","HypoxiaInO2"), c("Disease","Sick"))
    # c("Disease","LungParench") c("BirthAsphyxia", "Disease")
  }else if(network=="sang"){
    insurance <- loadRData("sangiovese.rda")
    true_amat <- amat(insurance)
  }else if(network=="mehra"){
    mehra <- loadRData("mehra.rda")
    true_amat <- amat(mehra)
  }else if(network=="mildew"){
    mildew <- loadRData("mildew.rda")
    true_amat <- amat(mildew)
    # c("Disease","LungParench") c("BirthAsphyxia", "Disease")
  }else if(network=="alarm"){
    insurance <- loadRData("alarm.rda")
    true_amat <- amat(insurance)
  }else if(network=="water"){
    water <- loadRData("water.rda")
    true_amat <- amat(water)
  }else if(network=="magic2"){
    insurance <- loadRData("magic2.rda")
    true_amat <- amat(insurance)
  }else if(network=="ecoli"){
    ecoli <- loadRData("ecoli.rda")
    true_amat <- amat(ecoli)
  }else if(network=="hf"){
    hf <- loadRData("hailfinder.rda")
    true_amat <- amat(hf)
  }else if(network=="magic"){
    magic <- loadRData("magic.rda")
    true_amat <- amat(magic)
  }else if(network=="hepar2"){
    hepar2<- loadRData("hepar2.rda")
    true_amat <- amat(hepar2)
  }else if(network=="pathfinder"){
    data<- loadRData("pathfinder.rda")
    true_amat <- amat(data)
  }else if(network=="win"){
    data<- loadRData("win.rda")
    true_amat <- amat(data)
  }else if(network=="water_rand_order"){
    data<- loadRData("water.rda")
    true_amat <- amat(data)
  }else if(network=="magic2_rand_order"){
    data<- loadRData("magic2.rda")
    true_amat <- amat(data)
  }else if(network=="mehra_rand_order"){
    data<- loadRData("mehra.rda")
    true_amat <- amat(data)
  }else if(network=="magic_rand_order"){
    data<- loadRData("magic.rda")
    true_amat <- amat(data)
  }else if(network=="mildew_rand_order"){
    data<- loadRData("mildew.rda")
    true_amat <- amat(data)
  }else if(network=="alarm_rand_order"){
    data<- loadRData("alarm.rda")
    true_amat <- amat(data)
  }else{
    print("else network 3 chosen; creating true amat")
    true_amat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0), 
                       c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
    colnames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
    rownames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
  }
  
  
  # create bnlearn obj of true dag
  true.dag = empty.graph(colnames(true_amat))
  amat(true.dag) = true_amat
  nodeN = length(true.dag$nodes)
  labels = paste("X", 1:nodeN, sep = "")
  if(network %in% c("water_rand_order", "magic2_rand_order", "mehra_rand_order", "magic_rand_order")){
    labels = rev(labels)
  }
  true.dag <- rename.nodes(true.dag, labels)
  
  return(true.dag)
}

get_metrics <- function(true.cpdag, init=NULL, edgedel_obj=NULL, edgedel_cpdag=NULL, edgedir_obj=NULL, edgedir_cpdag, edgerec_obj=NULL, edgerec_cpdag=NULL){
  # 1. init vs true
  init_shd <- bnlearn::shd(true=true.cpdag, learned=cpdag(init))
  init_f1 <- get_f1(true.cpdag, cpdag(init))[1]
  
  # 2. edgedel vs true
  shd_edgedel <- bnlearn::shd(true=true.cpdag, learned=edgedel_cpdag, wlbl=TRUE)
  f1_edgedel <- get_f1(true.cpdag, edgedel_cpdag)[1]
  skeleton_time <-attr(edgedel_obj[[1]],"gam_time")
  # learning_time <-attr(edgedel_obj[[1]],"skeleton_time")
  
  # 3. edgedir vs true
  edgedir_results <- analyze_cpdag(edgedir_cpdag, true.cpdag)
  edgedir_time <- c(edgedir_obj$gam_edge_time, edgedir_obj$edgeorient_time)
  
  
  # 4. edgerec vs true
  edgerec_results <- analyze_cpdag(edgerec_cpdag, true.cpdag)
  edgerec_time <- edgerec_obj$runtime
  
  all_results <- c(init_shd, init_f1, shd_edgedel, f1_edgedel, skeleton_time,
                   edgedir_results, edgedir_time, edgerec_results, edgerec_time)
  
  # a <- c('shd_pc', 'f1_pc', 'shd_skeleton', 'f1_skeleton', 'edgedel_time', 'shd_edgedir', 'f1_edgedir',
  #   'tp_edgedir', 'fp_edgedir', 'fn_edgedir', 'r_edgedir', 'ji_edgedir', 'noninv_fn_edgedir', 'edgedir_time', 'edgedir_test_time',
  #   'shd_final', 'f1_final', 'tp_final', 'fp_final', 'fn_final', 'r_final', 'ji_final', 'noninv_fn_final', 'oos_time', 'alg_time')
  return(all_results)
}


classify_edges <- function(init_cpdag, final_cpdag, nonlinear_edges){
  init_edges <- cbind(init_cpdag$arcs, "nonlinear"= rep(0,nrow(init_cpdag$arcs)))
  final_edges <- final_cpdag$arcs
  # get linear edges (i.e. from initial learning alg and preserved till the end)
  linear_edges <- bnlearn::compare(final_cpdag, init_cpdag,  arcs=TRUE)$tp
  linear_edges <- cbind(linear_edges, "nonlinear"= rep(0,nrow(linear_edges)))
  # get nonlinear edges (i.e. found via oriented, added procedures)
  # if(!is.null(edgerec_res)){
  #   edgerec_res <- edgerec_res[,c(1,2)]
  # }
  # colnames(nonlinear_edges) <- c("from", "to")
  if(!is.null(nonlinear_edges)){
    nonlinear_edges <- cbind(nonlinear_edges, "nonlinear"= rep(1,nrow(nonlinear_edges)))
    
    # left join to get linear/nonlinear status
    final_edges <- merge(x=final_edges,y=rbind(linear_edges, nonlinear_edges), 
                         by=c("from","to"), all.x=TRUE)
    final_edges <- final_edges[!duplicated(final_edges[,c(1,2)], fromLast=TRUE),]
    # if edge status is missing, then it's linear (not added from GAM procedures)
    final_edges[is.na(final_edges)] <- 0
  }else{
    final_edges <- cbind(final_edges, "nonlinear"= rep(0,nrow(final_edges)))
  }

  return(final_edges)
}


# for calculating the likelihood under GAM
est_normsd <- function(mod){
  est_sd <- sqrt(sum((mod$y-mod$fitted.values)^2)/mod$df.residual)
  return(est_sd)
}


cpdag2dag <- function(cpdag){
  # convert CPDAG to random DAG in equivalence class
  counter <- 0
  # need to supply ordering, so we pass in random order
  # permutate order until cpdag is converted into dag
  while(counter < 100){
    dag <- try(bnlearn::pdag2dag(cpdag, names(cpdag$nodes)[sample(1:length(names(cpdag$nodes)))]), silent=TRUE)
    if(class(dag) != "try-error"){
      return(dag)
    }
    counter <- counter+1
  }
  return(cpdag)
}

get_gam_loglik <- function(train_data, test_data, y, x, final_edges){
  # no parents
  # just use test data parameters to estimate log likelihood
  if(length(x)==0){
    # get loglike of train data
    y_train <- train_data
    ll_train <- dnorm(y_train, mean = mean(y_train), sd = sd(y_train), log = TRUE)
    # get loglik of test data
    y_data <- test_data
    ll_test <- dnorm(y_data, mean = mean(y_data), sd = sd(y_data), log = TRUE)
    return(c(sum(ll_train), sum(ll_test)))
  }
  
  # final_edges <- final_edges[final_edges["to"]==y,]
  # has parents, build GAM
  
  # number of base functions for smoothing function
  p <- dim(as.matrix(train_data))
  if(p[1]/p[2] < 3*10){
    k_basis <- ceiling(p[1]/(3*p[2]))
  }else{
    k_basis <- min(4, ceiling(p[2]/3))
  }
  
  covar <- "~"
  for(i in 1:length(x)){
    nonlinear <- as.numeric(final_edges[final_edges[,"from"]==x[i], "nonlinear"])
    if(nonlinear==0 || nonlinear=="0"){
      covar <- paste(covar, x[i], "+", sep="")
    }else{
      covar <- paste(covar, "s(", x[i], ", k=", k_basis, ", bs='cr') +", sep="")
    }
  }
  covar <- substr(covar, 1, nchar(covar)-1)
  fm <- as.formula(paste(y, covar, sep=""))
  gam.control.param = list(nthreads=4, epsilon=1e-07, maxit=150, mgcv.tol=1e-07)
  gam_mod = gam(formula = as.formula(fm), data = train_data, optimizer = c("outer","bfgs"),
                  select = TRUE, gam.control=gam.control.param, debug=TRUE)
  
  # predict on test data
  y_pred <- predict.gam(gam_mod, newdata=test_data)
  y_data <- test_data[,y]
  # print(mean(abs(train_data[,y]-gam_mod$fitted.values)))
  # print(mean(abs(y_pred-y_data)))

  # y_data <- gam_mod$y
  # estsd_mod <- est_normsd(gam_mod)
  # nobs_prior <- length(y_data)
  # prior_res <- gam_mod$residuals
  
  # find likelihood of train & test data under estimated parameters (mean & sd) under gam
  ll_gam_train <- logLik(gam_mod)
  
  est_sd <- sqrt(sum((y_data-mean(y_data))^2)/gam_mod$df.residual)
  est_sd <- sqrt(sum((y_data-mean(y_pred))^2)/gam_mod$df.residual)
  est_sd <- sd(y_pred)

  y_train_pred <- predict.gam(gam_mod, newdata=train_data)
  ll_gam_train <- dnorm(train_data[,y], mean = y_train_pred, sd = sd(y_train_pred), log = TRUE)
  ll_gam_test <- dnorm(y_data, mean = y_pred, sd = est_sd, log = TRUE)
  return(c(sum(ll_gam_train), sum(ll_gam_test)))
  
}


compute_dag_ll <- function(train_data, test_data, dag, final_edges){
  # get list of nodes
  joint_ll_train <- 0
  joint_ll_test <- 0
  
  train_data <- as.data.frame(train_data)
  test_data <- as.data.frame(test_data)

  for(curr_node in nodes(dag)){
    # print(curr_node)
    curr_node_pa <- dag$nodes[[curr_node]]$parents
    edge_set <- final_edges[final_edges[,"to"]==curr_node,,drop = FALSE]
    loglik_res <- get_gam_loglik(train_data[,c(curr_node, curr_node_pa)], test_data[,c(curr_node, curr_node_pa)], curr_node, curr_node_pa, edge_set)
    joint_ll_train <- joint_ll_train + loglik_res[1]
    joint_ll_test <- joint_ll_test + loglik_res[2]
    print(c(joint_ll_train/nrow(train_data), joint_ll_test/nrow(test_data)))
    # print(paste('ll_train: ', joint_ll_train))
    # print(paste('ll_test: ', joint_ll_test))
    # joint_ll <- sum(get_gam_loglik(train_data[,c(curr_node, curr_node_pa)], test_data[,c(curr_node, curr_node_pa)], curr_node, curr_node_pa, edge_set))
  }
  
  return(c(joint_ll_train, joint_ll_test))
}
# 
# classify_edges(ges_dag, oos_ges_res$gam_cpdag, ll_ges_cpdag$sigpair, oos_ges_res$added_edges)
# classify_edges(pc_cpdag, oos_pc_res$gam_cpdag, ll_pc_cpdag$sigpair, oos_pc_res$added_edges)
# classify_edges(ccdr_dag, osos_ccdr_res$gam_cpdag, ll_ccdr_cpdag$sigpair, oos_ccdr_res$added_edges)
# 
# 
# status <- try(bnlearn::pdag2dag(oos_ccdr_res$gam_cpdag, names(oos_ccdr_res$gam_cpdag$nodes)[sample(1:length(names(oos_ccdr_res$gam_cpdag$nodes)))]), silent=TRUE)
# # 
# compute_dag_ll(example_data, example_data, oos_ccdr_res$gam_cpdag, classify_edges(ccdr_dag, oos_ccdr_res$gam_cpdag, ll_ccdr_cpdag$sigpair, oos_ccdr_res$added_edges))
# compute_dag_ll(example_data, example_data, ges_final_dag, classify_edges(ges_dag, ges_final_dag, ll_ges_cpdag$sigpair, oos_ges_res$added_edges))
# 
# ges_final_dag <-cpdag2dag(oos_ges_res$gam_cpdag)
# classify_edges(ges_dag, cpdag2dag(oos_ges_res$gam_cpdag), ll_ges_cpdag$sigpair, oos_ges_res$added_edges)
# 


# get_chipseq_res(train_data, test_data, pc_cpdag, gam_pc_obj$dsep_set, ll_pc_cpdag$sigpair, oos_pc_res$added_edges, oos_pc_res$gam_cpdag)

get_chipseq_res <- function(train_data, test_data, init_cpdag, edgedel_res, edgedir_res, edgeadd_res, final_cpdag){
  # details about graph structure
  init_undir_edge <- ifelse(is.null(nrow(undirected.arcs(init_cpdag))), 0, nrow(undirected.arcs(init_cpdag)))
  n_undir_edge <- ifelse(is.null(nrow(undirected.arcs(final_cpdag))), 0, nrow(undirected.arcs(final_cpdag)))
  n_nonlinear_del <- ifelse(is.null(nrow(edgedel_res)), 0, nrow(edgedel_res))
  n_nonlinear_dir <- ifelse(is.null(nrow(edgedir_res)), 0, nrow(edgedir_res))
  n_nonlinear_add <- ifelse(is.null(nrow(edgeadd_res)), 0, nrow(edgeadd_res))
  
  # get edge classification for final cpdag
  if(nrow(edgeadd_res) > 0){
    nonlinear_edges <- rbind(edgedir_res, edgeadd_res[,c("from", "to")])
  }else{
    nonlinear_edges <- edgedir_res
  }
  
  # extend cpdag to random dag in equivalence class
  if(init_undir_edge > 0){
    init_dag <- cpdag2dag(init_cpdag)
  }else{
    init_dag <- init_cpdag
  }
  final_dag <- cpdag2dag(final_cpdag)
  
  # classify linear/nonlinear edges in final dag
  classified_edges <- classify_edges(init_cpdag, final_dag, nonlinear_edges)
  n_nonlinear_edges <- nrow(nonlinear_edges)
  
  # get joint likelihood
  test_joint_ll_final <- compute_dag_ll(train_data, test_data, final_dag, classified_edges)
  test_joint_ll_init <- compute_dag_ll(train_data, test_data, init_dag, classify_edges(init_dag, init_dag, NULL))
  
  return(c(nrow(train_data), nrow(test_data), nrow(init_cpdag$arcs)-init_undir_edge, init_undir_edge/2, 
           test_joint_ll_init, test_joint_ll_final, nrow(final_cpdag$arcs)-n_undir_edge, 
           n_undir_edge/2, n_nonlinear_del, n_nonlinear_dir, n_nonlinear_add))
  
  # c(test_alpha, pc_alpha, linear_alg, fold_num, nrow(test_data), 
  #   nrow(pc_cpdag$arcs)-init_undir_edge, init_undir_edge/2,  
  #   test_joint_ll_init, test_joint_ll, nrow(final_cpdag$arcs)-n_undir_edge, n_undir_edge/2,
  #   n_nonlinear_del, n_nonlinear_dir, n_nonlinear_add)
}


