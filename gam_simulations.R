library(bnlearn)
library(RBGL)
library(pcalg)
library(mgcv)
library(devtools)
library(phsl)
library(purrr)
setwd("~/Documents/UCLA/Research/GAM")
source("datagen.R")
source("gam_pPC.R")


# function to simulate networks
# input: network ref number, number of simulations to do, computing cluster for PC alg, invertible edge function
# output: table containing results on edge detection & runtime
simulate_network <- function(network=1, n_simulations=1, cluster=NULL, inv.method="sl"){
  # declare network & inverse arcs
  if(network==1){
    # EXAMPLE 1
    true_amat <- rbind(c(0,1,1,0), c(0,0,0,0), c(0,0,0,1), c(0,0,0,0))
    colnames(true_amat) <- c("x1","x2","x3","x4")
    rownames(true_amat) <- c("x1","x2","x3","x4")
    ninv.arcs <- rbind(c("x1", "x2"), c("x1", "x3"))
  }else if(network==2){
    # EXAMPLE 2
    # same network as network 3, but with node 6 removed
    true_amat <- rbind(c(0,1,0,0,0), c(0,0,1,1,0),c(0,0,0,0,0), 
                        c(0,0,0,0,0), c(0,0,0,1,0))
    colnames(true_amat) <- c("x1","x2","x3","x4","x5")
    rownames(true_amat) <- c("x1","x2","x3","x4","x5")
    # ninv.arcs <- rbind(c("x2", "x3"))
    ninv.arcs <- c("x2", "x4")
  }else if(network==3){
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
  }else if(network=="asia2"){
    asia <- loadRData("asia.rda")
    true_amat <- amat(asia)
    ninv.arcs <- rbind(c("tub", "either"),c("smoke","lung"), c("smoke","bronc"))
  }else if(network=="sachs"){
    sachs <- loadRData("sachs.rda")
    true_amat <- amat(sachs)
    ninv.arcs <- rbind(c("PKC", "Mek"),c("PIP3","PIP2"), c("PKA","P38"), c("PKA","Erk"))
  }
  else if(network=="child"){
    child <- loadRData("child.rda")
    true_amat <- amat(child)
    ninv.arcs <- rbind(c("BirthAsphyxia", "Disease"),c("LungParench","ChestXray"), 
                       c("CardiacMixing","HypoxiaInO2"), c("Disease","Sick"))
    # c("Disease","LungParench") c("BirthAsphyxia", "Disease")
  }else if(network=="mildew"){
    mildew <- loadRData("mildew.rda")
    true_amat <- amat(mildew)
    ninv.arcs <- rbind(c("temp_1", "foto_1"),c("nedboer_1","mikro_1"), 
                       c("lai_2","mikro_2"), c("meldug_3","lai_3"),
                       c("meldug_4","lai_4"), c("dm_2","dm_3"),
                       c("temp_2","foto_2"))
    # c("Disease","LungParench") c("BirthAsphyxia", "Disease")
  }
  
  
  true_bn = empty.graph(colnames(true_amat))
  amat(true_bn) = true_amat
  
  #create data frame
  results <- as.data.frame(matrix(ncol = 10, nrow = 0))

  # formatting folder name according to invertible func type
  if(inv.method=="sl"){
    folder_name = paste("network", network, sep="")
  }else{
    folder_name = paste("network", network, "_random", sep="")
  }
  
  # generate datasets and learn networks
  for(i in c(1:n_simulations)){
    gc()
    # generate dataset
    test = DAGdatagen(n = 1000, dag.bn = true_bn, inv.method = inv.method, ninv.arcs = ninv.arcs, ninv.method = "nl3")
    toy_data = as.data.frame(test$DAGdata)
    # save dataset
    saveRDS(toy_data, file = paste(folder_name, "/example", i,".Rds", sep=""))
    # learn skeleton
    start_time <- Sys.time()
    temp_skeleton <- GAM_edgedel(toy_data, alpha=0.05, standardized=FALSE, cluster=cluster)
    end_time <- Sys.time()
  
    
    # save deleted edges
    del_edges <- temp_skeleton[[2]]
    if(is.null(del_edges) == FALSE){
      del_edges <- as.matrix(del_edges[,c(1,2)])
    }
    # convert PC skeleton to cpdag
    pc_cpdag = bnlearn::empty.graph(colnames(true_amat))
    amat(pc_cpdag) = get_cpdag(toy_data, temp_skeleton[[3]])
    # convert GAM skeleton to cpdag
    gam_cpdag = bnlearn::empty.graph(colnames(true_amat))
    amat(gam_cpdag) = get_cpdag(toy_data, temp_skeleton[[1]], blacklist=del_edges) # make sure deleted edges aren't added back
    
    # save to results
    total_time <- as.numeric(end_time - start_time, unit = "secs")
    gam_time <- attr(temp_skeleton[[1]], "gam_time")
    # compare PC/PC-GAM skeletons w/ true skeleton
    pc_bn_diff <- skeleton_eval(true_bn, pc_cpdag)
    gam_bn_diff <- skeleton_eval(true_bn, gam_cpdag)
    gam_del_results <- (pc_bn_diff - gam_bn_diff)[c(2,1)]
    # save results into larger dataframe
    results <- rbind(results, c(pc_bn_diff, gam_bn_diff, gam_del_results, total_time, gam_time))
  }
  names(results) <- c('pc_truepos', 'pc_falsepos', 'pc_falseneg', 'gam_truepos', 'gam_falsepos', 'gam_falseneg', 'gam_delete_tp', 'gam_delete_fn', 'total_time', 'gam_time')
  write.csv(results, file= paste(folder_name, "/results.csv", sep=""))
  
  return(results)
}

# perform simulations
results <- simulate_network(n_simulations = 300, network="sachs", inv.method="sl")





