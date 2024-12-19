library(bnlearn)
library(RBGL)
library(pcalg)
library(mgcv)
library(devtools)
library(igraph)
# load package to use pPC function
# install_github("jirehhuang/phsl")
# print("installed all packages; phsl now")
library(phsl)
suppressWarnings(suppressMessages(library(RBGL)))
suppressWarnings(suppressMessages(library(devtools)))
suppressWarnings(suppressMessages(library(igraph)))


# GAM RESOURCE PAGE: https://content.naic.org/sites/default/files/call_materials/NAIC%20GAM%20Presentation%20-%20Insurance%20Summit%202021%20-Final.pdf


# purpose: gets parents (directed/undirected edges) of node
# input: skeleton from pPC_skeleton()/list of nodes from bnlearn obj, node
# output: directed/undirected parents of node (with labels of being pa or not) or FALSE if no parents/node
get_curr_pa = function(x, node) {
  # check if node is in skeleton
  if(node %in% names(x)){
    pa_set <- x[[node]]$nbr[!x[[node]]$nbr %in% x[[node]]$children]
    is_confirmed_pa <- pa_set %in% x[[node]]$parents
  }else{
    return(NULL)
  }
  if(length(pa_set) > 0){
    return(cbind(pa_set, is_confirmed_pa))
  }
  else{
    return(NULL)
  }
}

# purpose: removes edge from skeleton & updates properties
# input: pPC skeleton, nodes to update, p-value of coef from GAM
# output: skeleton (modifies skeletons then returns it)
drop_arc = function(dag, node1, node2, p_value, dsep_set){
  # update neighbor & markov blanket for node1, node2
  dag[[node1]]$nbr <- dag[[node1]]$nbr[dag[[node1]]$nbr != node2]
  dag[[node2]]$nbr <- dag[[node2]]$nbr[dag[[node2]]$nbr != node1]
  dag[[node1]]$mb <- dag[[node1]]$mb[dag[[node1]]$mb != node2]
  dag[[node2]]$mb <- dag[[node2]]$mb[dag[[node2]]$mb != node1]
  
  # update p-value & dsep.set in attributes
  if(is.null(attr(dag, "dsep.set"))){
    return(dag)
  }
  
  p_value <- unname(p_value)
  node_pairs <- attr(dag, "dsep.set")
  pair_to_drop <- sort(c(node1, node2))
  for(i in 1:length(node_pairs)){
    if(all(pair_to_drop == sort(node_pairs[[i]]$arc))){
      attr(dag, "dsep.set")[[i]]$p.value <- p_value
      # add found conditional indep. relation to current dsep_set
      attr(dag, "dsep.set")[[i]]$dsep.set <- c(attr(dag, "dsep.set")[[i]]$dsep.set, dsep_set)
    }
  }
  return(dag)
}

# purpose: builds GAM formula
# input: output nodeY, predictors vector x
# output: formula for GAM

get_gam_formula <- function(y, x, k, smooth_fx="tp"){
  covar <- "~"
  for(i in 1:length(x)){
    if(length(x)==0){
      break
    }
    covar <- paste(covar, "s(", x[i], ", k=", k, ", bs='", smooth_fx, "') +", sep="")
  }
  covar <- substr(covar, 1, nchar(covar)-1)
  fm <- as.formula(paste(y, covar, sep=""))
  return(fm)
}



# purpose: uses GAM to detect non-linear relations & further delete edges from skeleton
# input: data
# output: list of
# 1. GAM reduced skeleton
# 2. deleted edges from GAM and conditioning set
# 3. PC estimated skeleton
# note: pPC object can be converted into bnlearn object
prune_dag <- function(dat, init_dag, alpha=0.05, standardized=FALSE, cluster=NULL, nodelabels=colnames(dat)){
  
  start_time = Sys.time()
  # add names to nodes
  nodelabels = if(is.null(nodelabels)){
    paste("X", 1:ncol(dat), sep = "")
    colnames(dat) = nodelabels 
  }else{
    nodelabels <- nodelabels
  }
  
  # standardize data if needed
  if(standardized==FALSE){
    dat <- apply(dat, 2, function(x) if(is.numeric(x)){scale(x, center=TRUE, scale=TRUE)} else x)
  }
  dat <- as.data.frame(dat)
  
  # save for comparison later
  curr_struct <- init_dag$nodes
  cond_set <- as.data.frame(matrix(ncol = 0, nrow = 0))
  full_dsep_vec <- c()
  
  
  # number of base functions for smoothing function
  p <- dim(as.matrix(dat))
  if(p[1]/p[2] < 3*10){
    k_basis <- ceiling(p[1]/(3*p[2]))
  }else{
    k_basis <- min(7, ceiling(p[2]/5))
  }
  
  gam.control.param = list(nthreads=4, epsilon=1e-04, maxit=100, mgcv.tol=1e-04)
  
  # try to extend the CPDAG to a DAG and estimate a topological order
  temp_dag <- tryCatch(bnlearn::cextend(init_dag), error=function(e) e, warning=function(w) w)
  if(is(temp_dag, 'error')){
    attempt <- 1
    while(TRUE && attempt < 6){
      attempt <- attempt + 1
      temp_dag <- tryCatch(bnlearn::pdag2dag(init_dag, ordering=nodelabels[sample(c(1:length(nodelabels)))]),
                           error=function(e) e, warning=function(w) w)
      if(!is(temp_dag, 'error')) break
    }
    # if no DAG extension is possible, then use random order
    node_order <- nodelabels[sample(c(1:length(nodelabels)))]
  }else{
    node_order <- names(topological.sort(as.igraph(temp_dag)))
  }
  
  for(curr_node in node_order){
    # this include parents and neighbors (undirected edges)
    curr_node_pa <- get_curr_pa(curr_struct, curr_node)
    if(is.null(curr_node_pa)){
      next # no parents/no node
    }
    is_confirmed_pa <- curr_node_pa[,2] # col of true/false
    curr_node_pa <- curr_node_pa[,1] # col of parent/neighbor nodes
    # run gam model
    fm <- get_gam_formula(curr_node, curr_node_pa, k_basis)
    m1 <- gam(fm, method = "REML", optimizer = c("outer","newton"), select = TRUE, data=as.data.frame(dat), gam.control=gam.control.param)
    # check if any variable is insignificant, remove edge if so
    m1_pval <- summary(m1)$s.table[,4]
    insig_edges <- which(m1_pval >= alpha)
    insig_edges_pval <- unname(m1_pval[insig_edges])

    
    # iteratively delete insignificant edges
    # delete from bnlearn object, then update "curr_struct" with new structure
    for(j in 1:length(insig_edges)){
      if(length(insig_edges) == 0){
        break
      }
      insig_node <- curr_node_pa[insig_edges[j]]
      # update skeleton by dropping arc & updating conditioning set, with info on if edge was directed or undirected
      dsep_set <- curr_node_pa[-insig_edges]
      cond_set <- rbind(cond_set, c(insig_node, curr_node, is_confirmed_pa[insig_edges[j]]))
      full_dsep_vec <- c(full_dsep_vec, list(dsep_set))
      # update the bnlearn obj, depending on relation (check if insig node is parent of neighbor)
      if(is_confirmed_pa[insig_edges[j]]){
        init_dag <- drop.arc(init_dag, from=insig_node, to=curr_node)
      }
    }
    curr_struct <- init_dag$nodes
  }
  # rename conditioning set
  if(nrow(cond_set)>0){
    colnames(cond_set) <- c("from", "to", "directed_edge")
    cond_set$dsep <- full_dsep_vec   
    # process insignificant undirected edges
    if(sum(cond_set[,3]=="FALSE") > 0){
      cond_set_nbr <- cond_set[cond_set[,3]=="FALSE", ]
      for(i in c(1:nrow(cond_set_nbr))){
        temp_edge <- c(cond_set_nbr[i, 2], cond_set_nbr[i, 1], cond_set_nbr[i, 3], list(cond_set_nbr[i, 4]))
        init_dag <- drop.arc(init_dag, from=cond_set_nbr[i, 1], to=cond_set_nbr[i, 2])
        # add opposite direction to cond set
        cond_set <- rbind(cond_set, temp_edge)
      }
    }
    # add everything to blacklist
    bl <- cond_set[,c(1,2)]
    bl <- matrix(bl, ncol = 2, byrow=TRUE, dimnames = list(NULL, c("from", "to")))
    if(is.null(init_dag$learning$blacklist)){
      init_dag$learning$blacklist <- bl
    }else{
      init_dag$learning$blacklist <- rbind(init_dag$learning$blacklist, bl)
    }
  }

  # record time for GAM procedure
  end_time = Sys.time()
  total_time <- as.numeric(end_time - start_time, unit = "secs")
  
  return(list(cpdag=init_dag, cond_set=cond_set, total_time=total_time))
}




# function to plot graphs
plot_network <- function(true_bn=NULL, pc_skeleton=NULL, gam_skeleton=NULL){
  par(mfrow = c(1,3))
  plot(true_bn, main="True BN")
  plot(pc_skeleton, main="Post-PC, CPDAG Estimate")
  plot(gam_skeleton, main="Post-GAM, CPDAG Estimate")
}



# function to load RDA file containing asia network
# https://www.bnlearn.com/bnrepository/discrete-small.html#asia
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



# function to compare true skeleton and estimated skeleton
# input: true network, estimated cpdag
# output: number of edges deleted correctly & incorrectly
skeleton_eval <- function(true_bn, est_cpdag){
  # unname(unlist(bnlearn::compare(bnlearn:::skeleton(true_bn), bnlearn:::skeleton(pc_cpdag))))
  dag_diff <- bnlearn::compare(bnlearn:::skeleton(true_bn), bnlearn:::skeleton(est_cpdag))
  return(c(dag_diff$tp, dag_diff$fp, dag_diff$fn))
  
}



