
source("generate_data.R")
source("estimate_cpdag.R")
source("orient_edges.R")
source("delete_edges.R")
source("finalize_dag.R")
source('util_functions.R')



#################################
# specify arguments
#################################
args = c("mildew", 1000, "random", 0.05, "55")
network = as.character(args[1])
n_samples = as.integer(args[2])
data_fx = as.character(args[3])
alpha = as.double(args[4])
approach = as.character(args[5])



########################################################
# load in dag and cpdag structures
########################################################
true.dag <- get_true_dag(network)
true.cpdag <- true.dag
if(data_fx == "linear") true.cpdag <- cpdag(true.dag)



#############################################################################
# determine data generating mechanism; overrides passed arguments
#############################################################################
if(data_fx == "linear"){
  ninv.perct <- 0
  nonlinear.perct <- 0
}else if(data_fx %in% c("cubic", "quadratic", "gp", "ninv")){
  ninv.perct <- 1
  nonlinear.perct <- 1
}else{
  ninv.perct <- 0
  nonlinear.perct <- 1
}



########################################################
# generate data
########################################################
rand_dag = DAGdatagen(n = n_samples, dag_amat = amat(true.cpdag), nodeN = length(nodes(true.cpdag)),
                      ninv.perct = ninv.perct, nonlinear.perct = nonlinear.perct,
                      inv.method=data_fx, ninv.method=data_fx, se=1, labels=nodes(true.cpdag))
example_data = rand_dag$DAGdata
example_data  <- apply(example_data, 2, function(x) if(is.numeric(x)){scale(x, center=TRUE, scale=TRUE)} else x)



########################################################
# run algorithm
########################################################
# 1. estimate skeleton
init_graph <- estimate_cpdag(data = example_data, learning_alg = "PC", standardized = TRUE,
                                     pc_params = list(alpha=alpha, alpha_dense=0.25, test='cor'))
init_amat <- init_graph$amat
# 2. edge orientation on cpdag
oriented_graph <- orient_edge(data = example_data, curr_amat=init_amat, alpha=0.05,
                                 train_ratio = 0.5, ll_test_approach = approach)
oriented_amat <- oriented_graph$amat
# 3. edge deletion
pruned_graph <- prune_dag(data = example_data, curr_amat=oriented_amat, alpha=0.00001)
pruned_amat <- pruned_graph$amat
# 4. extend CPDAG to DAG
print("extending to dag")
final_graph <- finalize_dag(data = example_data, curr_amat=pruned_amat, alpha=0.05, train_ratio = 0.5, ll_test_approach = approach)



########################################################
# evaluate learned dag(s)
########################################################
learning_info <- c(network, ncol(example_data), nrow(directed.arcs(true.cpdag)), nrow(undirected.arcs(true.cpdag))/2, n_samples, data_fx)
learning_res <- get_f1(amat(true.cpdag), list(init_amat, oriented_amat, pruned_amat, final_graph$amat))
learning_res <- c(learning_res[1, ], init_graph$total_time,
                  learning_res[2, ], oriented_graph$total_time,
                  learning_res[3, ], pruned_graph$total_time,
                  learning_res[4, ], final_graph$total_time)



########################################################
# display results
########################################################
metrics <- c("shd", "f1", "tp", "fp", "fn", "wrong_dir", "wrong_dir_undir", "n_undir", "time")
steps <- c("init", "orient", "delete", "final")
main_metrics <- paste(rep(metrics, length(steps)), rep(steps, each=length(metrics)), sep="_")
df_names <- c("network", "n_nodes", "n_dir_edges", "n_undir_edges", "n_samples", "data_fx", main_metrics)
alg_res_df <- as.data.frame(matrix(ncol = length(df_names), nrow = 0))
alg_res_df <- rbind(alg_res_df, c(learning_info, learning_res))
names(alg_res_df) <- df_names


print(alg_res_df)
  
