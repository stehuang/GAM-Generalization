library(bnlearn)
library(RBGL)
library(pcalg)
library(mgcv)
library(devtools)
# load package to use pPC function
# install_github("jirehhuang/phsl")
library(phsl)

# GAM RESOURCE PAGE: https://content.naic.org/sites/default/files/call_materials/NAIC%20GAM%20Presentation%20-%20Insurance%20Summit%202021%20-Final.pdf


# purpose: gets neighbors of node in pPC skeleton obj
# bnlearn implementation: https://github.com/cran/bnlearn/blob/master/R/frontend-nodes.R
# input: skeleton from pPC_skeleton(), node
# output: neighbors of node in skeleton
get_nbr = function(x, node) {
  # check if node is in skeleton
  if (node %in% names(x))
    return(x[[node]]$nbr)
  else
    print("node not found in skeleton")
  break
  
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
get_gam_formula <- function(y, x){
  covar <- "~"
  for(i in 1:length(x)){
    if(length(x)==0){
      break
    }
    covar <- paste(covar, "s(", x[i], ", k=-1, bs='cr') +", sep="")
  }
  covar <- substr(covar, 1, nchar(covar)-1)
  fm <- as.formula(paste(y, covar, sep=""))
  return(fm)
}



# purpose: uses GAM to detect non-linear relations & further delete edges from skeleton
# input: data
# output: list of
   # 1. GAM reduced reduced skeleton
   # 2. deleted edges from GAM and conditioning set
   # 3. PC estimated skeleton
# note: pPC object can be converted into bnlearn object
GAM_edgedel <- function(data, alpha=0.05, standardized=FALSE, cluster=NULL){
  # standardize data if needed
  if(standardized==FALSE){
    data <- lapply(data, function(x) if(is.numeric(x)){scale(x, center=TRUE, scale=TRUE)} else x)
  }
  data <- as.data.frame(data)
  
  # gam parameters
  gam.control.param = list(nthreads=4, epsilon=1e-07, maxit=150, mgcv.tol=1e-07)
  
  # run pPC algorithm to estimate skeleton
  pc_skeleton <- phsl:::ppc_skeleton(x = data, cluster = cluster, whitelist = NULL, blacklist = NULL,
                                test = "cor", alpha = alpha, B = NULL,
                                complete = bnlearn:::check.data(data)$complete.nodes,
                                sort_pval = FALSE, max_groups = 20, debug = TRUE, max.sx = ncol(data)-2)
  # save for comparison later
  gam_skeleton <- pc_skeleton
  cond_set <- as.data.frame(matrix(ncol = 0, nrow = 0))
  full_dsep_vec <- c()

  start_time = Sys.time()
  for(curr_node in names(gam_skeleton)){
    curr_nbr <- get_nbr(gam_skeleton, curr_node)
    # only consider nodes with more than 2 neighbors
    if(length(curr_nbr) <= 2) {next}
    
    # run gam model
    fm <- get_gam_formula(curr_node, curr_nbr)
    m1 <- gam(fm, method = "REML", optimizer = c("outer","newton"), 
              select = TRUE, data=as.data.frame(data), gam.control=gam.control.param)
  
    
    # check if any variable is insignificant, remove edge if so
    m1_pval <- summary(m1)$s.table[,4]
    insig_edges <- which(m1_pval > alpha)
    insig_edges_pval <- unname(m1_pval[insig_edges])

    # iteratively delete insignificant edges
    for(j in 1:length(insig_edges)){
      if(length(insig_edges) == 0){
        break
      }
      insig_node <- curr_nbr[insig_edges[j]]
      # update skeleton by dropping arc & updating conditioning set
      dsep_set <- curr_nbr[-insig_edges]
      gam_skeleton <- drop_arc(gam_skeleton, curr_node, insig_node, insig_edges_pval[j], dsep_set)
      cond_set <- rbind(cond_set, c(curr_node, insig_node))
      full_dsep_vec <- c(full_dsep_vec, list(dsep_set))
    }
  }
  end_time = Sys.time()
  # rename conditioning set
  if(nrow(cond_set)>0){
    colnames(cond_set) <- c("n1", "n2")
    cond_set$dsep <- full_dsep_vec   
  }else{
    cond_set = NULL
  }
  # record time for GAM procedure
  attr(gam_skeleton , "gam_time") <- as.numeric(end_time - start_time, unit = "secs")
  
  return(list(gam_skeleton, cond_set, pc_skeleton))
}


# purpose: performs v-structure detection & edge orientation to get PDAG
# input: data, skeleton
# output: adjacency matrix of learned pdag
get_cpdag <- function(data, skeleton, blacklist=NULL, 
                     whitelist=NULL, test='cor', alpha=0.05, debug=FALSE){
  # this function does v-structure detection first, then edge orientation
  cpdag <- bnlearn:::learn.arc.directions(x=data, local.structure = skeleton, 
                                          blacklist=blacklist, whitelist=whitelist, 
                                          test=test, alpha=alpha, max.sx = NULL, debug=debug)
  # jireh's implementation
  # pdag_ppc <- phsl:::orient(x=toy_data, local.structure=temp_skeleton[[1]], 
  #                           score='bic-g', blacklist=blacklist, whitelist=c())
  nodes <- names(data)
  # convert arcs to adjacency matrix
  Amat <- bnlearn:::arcs2amat(cpdag$arcs, nodes)
  return(Amat)
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


# bnlearns PC alg, for debugging purposes
# bnlearn:::bnlearn(x=toy_data, test = 'cor', alpha = 0.05, method = "pc.stable", max.sx = 2, debug = TRUE, undirected = TRUE) 


# # EXAMPLE 1
# # true dag: x2 <- x1 -> x3 -> x4
# set.seed(2021)
# true_amat <- rbind(c(0,1,1,0), c(0,0,0,0), c(0,0,0,1), c(0,0,0,0))
# colnames(true_amat) <- c("x1","x2","x3","x4")
# rownames(true_amat) <- c("x1","x2","x3","x4")

# formulas used to generate network
# x1 <- runif(500)
# x2 <- x1^3 + x1*0.5 + rnorm(500)/100
# x3 <- x1^3 + rnorm(500)/100
# x4 <- x3 + rnorm(500)/100

# toy_data <- read.csv("network1.csv")
# cor(toy_data)
# learn skeleton
# temp_skeleton <- GAM_edgedel(toy_data, alpha=0.05, standardized=FALSE)


# # output deleted edge
# del_edges <- temp_skeleton[[2]]
# 
# # convert skeleton to pdag
# # blacklist is removed edges via GAM
# true_bn = empty.graph(colnames(true_amat))
# amat(true_bn) = true_amat
# pc_skeleton = empty.graph(colnames(true_amat))
# amat(pc_skeleton) = get_cpdag(toy_data, temp_skeleton[[3]])
# gam_skeleton = empty.graph(colnames(true_amat))
# amat(gam_skeleton) = get_cpdag(toy_data, temp_skeleton[[1]], blacklist=del_edges[,c(1,2)])
# 
# plot_network(true_bn, pc_skeleton, gam_skeleton)
# 
# 
# 
# ### EXAMPLE 2
# # same network as network 3, but with node 6 removed
# true_amat <- rbind(c(0,1,0,0,0), c(0,0,1,1,0),c(0,0,0,0,0),
#                    c(0,0,0,0,0), c(0,0,0,1,0))
# colnames(true_amat) <- c("x1","x2","x3","x4","x5")
# rownames(true_amat) <- c("x1","x2","x3","x4","x5")
# 
# source("datagen.R")
# # declare causal network here, which contains relations
# true_bn = empty.graph(colnames(true_amat))
# amat(true_bn) = true_amat
# 
# # ninv.arcs <- rbind(c("x2", "x4"))
# # test = DAGdatagen(n = 500, dag.bn = true_bn, ninv.perct = .3, ninv.arcs = ninv.arcs, ninv.method = "nl3")
# # toy_data = as.data.frame(test$DAGdata)
# # plot(toy_data$x4~toy_data$x5)
# toy_data <- read.csv("network2.csv")
# 
# # run PC and GAM regression
# temp_skeleton <- GAM_edgedel(toy_data, alpha=0.05, standardized=FALSE)
# 
# # convert skeleton to pdag
# # blacklist is removed edges via GAM
# true_bn = empty.graph(colnames(true_amat))
# amat(true_bn) = true_amat
# pc_skeleton = empty.graph(colnames(true_amat))
# amat(pc_skeleton) = get_pdag(toy_data, temp_skeleton[[3]])
# gam_skeleton = empty.graph(colnames(true_amat))
# amat(gam_skeleton) = get_pdag(toy_data, temp_skeleton[[1]], blacklist=del_edges[,c(1,2)])
# 
# plot_network(true_bn, pc_skeleton, gam_skeleton)
# 
# 
# 
# 
# # EXAMPLE 3
# # true dag: x1->x2->x3, x2->x4<-x5, x4->x6
# true_amat <- rbind(c(0,1,0,0,0,0), c(0,0,1,1,0,0),c(0,0,0,0,0,0), 
#                    c(0,0,0,0,0,1), c(0,0,0,1,0,0), c(0,0,0,0,0,0))
# colnames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
# rownames(true_amat) <- c("x1","x2","x3","x4","x5","x6")
# 
# # generate data
# set.seed(2022)
# source("datagen.R")
# # declare causal network here, which contains relations
# true_bn = empty.graph(colnames(true_amat))
# amat(true_bn) = true_amat
# 
# # generate data that retains causal network structure
# # need to specify which arcs are non-linear
# set.seed(2022)
# ninv.arcs <- rbind(c("x4", "x2"))
# test = DAGdatagen(n = 1000, dag.bn = true_bn, ninv.perct = .2, 
#                   ninv.arcs = ninv.arcs, ninv.method = "nl3", se = 0.3)
# toy_data = as.data.frame(test$DAGdata)
# 
# # toy_data <- read.csv("network3.csv")
# cor(toy_data)
# par(mfrow = c(1,1))
# plot(toy_data$x4~toy_data$x2)
# 
# # learn skeleton
# temp_skeleton <- GAM_edgedel(toy_data, alpha=0.05, standardized = FALSE)
# # deleted edges via GAM regression
# del_edges <- temp_skeleton[[2]]
# 
# # creating CPDAGs from skeleton
# true_bn = empty.graph(colnames(true_amat))
# amat(true_bn) = true_amat
# pc_skeleton = empty.graph(colnames(true_amat))
# amat(pc_skeleton) = get_cpdag(toy_data, temp_skeleton[[3]])
# gam_skeleton = empty.graph(colnames(true_amat))
# amat(gam_skeleton) = get_cpdag(toy_data, temp_skeleton[[1]], blacklist=del_edges[,c(1,2)])
# 
# # plot the true network and 2 estimated CPDAGs
# plot_network(true_bn, pc_skeleton, gam_skeleton)


