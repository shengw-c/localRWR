library(memoise)

#' Degree normalization of adjacency matrix
#' 
#' @param graph igraph object
#' @param normalize_strategy strategy of normalization, one of ["Row", "Column", "Laplacian"]
#' @param weighted whether it's a weighted graph, default TRUE
#' 
#' @return Normalized matrix
#' 
#' @example 
#' normal_W(graph, normalize_strategy = "Row", weighted = TRUE)
#' 
normal_W <- function(graph, normalize_strategy, weighted = TRUE) {
  
  message(paste0(normalize_strategy, "-normalization..."))
  
  ifweight = ifelse(weighted==TRUE, "weighted", "unweighted")
  message(paste0("Using ", ifweight, " graph..."))
  
  ###For row or column wise normalization
  if (normalize_strategy %in% c("Row", "Column")){
    ##whether using edge weights to generate adjacency matrix
    if (weighted) {
      adj_matrix <- as_adjacency_matrix(graph, sparse = FALSE, attr = "weight")
    }else{
      adj_matrix <- as_adjacency_matrix(graph, sparse = FALSE)
    }
    
    ##calculate row or column sum based on the normalization strategy
    ##the results should be identical for undirected graph
    D = if (normalize_strategy=='Row') diag(Cpp_rowSums(adj_matrix)) else diag(Cpp_colSums(adj_matrix))
    
    ##Check for Zero-Degree Nodes
    if (any(diag(D) == 0)) {
      warning("Zero-degree nodes present. Consider removing or adjusting them.")
    }
    
    ##Inverse of Degree Matrix (D^-1)
    D_inv = Arma_Inv(D)
    
    ##row-wise normalization: W = D^-1A
    ##column-wise normalization: W = AD^-1
    W = if (normalize_strategy=='Row') D_inv %*% adj_matrix else adj_matrix %*% D_inv
  }else{
    ###using igraph function laplacian_matrix to generate the normalized laplacian matrix
    message("Normalized Laplacian normalization...")
    
    if (weighted) {
      message("Using weighted graph...")
      W_Lnorm <- laplacian_matrix(graph, sparse = FALSE, weights=NULL, normalized = TRUE)
    }else{
      message("Using unweighted graph...")
      W_Lnorm <- laplacian_matrix(graph, sparse = FALSE, weights = NA, normalized = TRUE)
    }
    
    ##zero diagonal numbers
    I <- diag(nrow(W_Lnorm))  # Create the identity matrix
    W <- I - W_Lnorm
  }
  return(W)
}

###create a memory cache object to cache the normalized adjacency matrix, to save time
###default 3GB, adjust as needed
if (!exists("cache_mem_gb")) {
  cache_mem_gb = 3 * 1024^3
}
print(paste0("Using ", cache_mem_gb / 1024^3, "GB to store cache..."))
cm <- cachem::cache_mem(max_size = cache_mem_gb, max_age = 12*60*60) ## default memory for 12 hours

###wrap normal_W function with memoise (memoise package), make sure the output of this function could be 
###cached into above cache object, if appropriate
normal_W_cache = memoise(normal_W, cache = cm)

#' Local implementation of network propagation
#' 
#' @description Local implementation of network propagation, including both random walk with restart and heat diffusion. But due to the complexity of calculation, HD is too slow at this point, so only RWR is practical right now
#' 
#' @param graph igraph object
#' @param init_scores initial scores, should be with identical order with the vertices
#' @param int_scores intermediate scores, only needed if you need to manipulate the scores during the process. This should also be with identical order with the vertices
#' @param kernel random walk with restart (RWR) or heat diffusion (HD), only RWR is practical at this point
#' @param normalize_strategy strategy of normalization, one of ["Row", "Column", "Laplacian"]
#' @param weighted whether it's a weighted graph, default TRUE
#' @param param restart probability for RWR and spreading parameter for HD
#' @param niter total run of iteration, default 1000, but can be any positive integer
#' @param delta threshold of convergence
#' 
#' @return a list of two vectors:
#'         vector: scores after network propagation
#'         isconverged: whether the network propagation achieves convergence within given iteration
#'
#' @export
np_local <- function(graph,
                     init_scores,
                     int_scores = NULL,
                     kernel = c("RWR", "HD"), 
                     normalize_strategy = c("Row", "Column", "Laplacian"),
                     weighted = TRUE,
                     param = 0.15,
                     niter = 1000,
                     delta = 1e-6) {
  
  ###to make it run faster, 
  ###1. only allows Row-normalization for adjacency matrix
  ##normalize_strategy="Row"
  ###2. Only support random walk with restart
  kernel = "RWR"
  
  ##Ensure scores are aligned with nodes and scaled
  if (length(init_scores) != vcount(graph)) {
    stop("Number of scores must match the number of nodes in the graph.")
  }
  
  W = normal_W_cache(graph = graph, normalize_strategy = normalize_strategy, weighted = weighted)
  
  if (kernel=='RWR') {
    rwr_cpp(adj_mat = W, 
            init_scores = init_scores, 
            int_scores = int_scores,
            restart_prop = param, 
            num_iter = niter,
            delta = delta)
  }else{
    if (! sum(init_scores) == 1){
      F = init_scores / sum(init_scores)
    }else{
      F = init_scores
    }
    sparse_transition_matrix <- as(transition_matrix, "sparseMatrix")  # Convert to sparse format
    Eigen_expm(-param * sparse_transition_matrix) %*% F
  }
}