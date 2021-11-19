###########################################
########## Backbone of a network ##########
###########################################

## Choose different backbone levels: alpha = 10^-2, 10^-3, 10^-4, etc.
## The backbone is based on "Xu, M. and Liang, S., Input-output networks offer new insights 
## of economic structure".
## Alternative is to use a fixed threshold.
## If both alpha and a fixed threshold are given, the fixed threshold is used.
## Weights can be kept or removed.
## Direction can be kept or removed.

netbackbone <- function(M, alpha = 0.05, thres.fixed = NA, weighted = TRUE, directed = TRUE){
  ## if using the fixed threshold
  n <- dim(M)[1]
  if (is.na(thres.fixed) == FALSE){
    thres <- thres.fixed
    M[M < thres] <- 0
  }
  
  ## if using filter method
  if (missing(thres.fixed)|is.na(thres.fixed)){
    ## normalized out-weight and in-weight
    w_norm_out <- M/rowSums(M) ## produce NaNs for the nodes with no outgoing edges
    w_norm_in <- t(t(M)/colSums(M)) ## produce NaNs for the nodes with no incoming edges
    ## adjacency matrix ignoring weight
    A <- M
    A[A > 0] <- 1
    A_out <- rowSums(A)
    A_in <- colSums(A)
    ## alpha matrices for filtration
    alpha_out <- (1 - w_norm_out)^(A_out - 1) ## A_out = 1 always produces alpha_out = 1
    alpha_in <- t((t(1 - w_norm_in))^(A_in - 1)) ## A_in = 1 always produces alpha_in = 1
    alpha_out[alpha_out < alpha] <- 0
    alpha_out[alpha_out > alpha] <- 1
    alpha_out <- 1 - alpha_out
    alpha_in[alpha_in < alpha] <- 0
    alpha_in[alpha_in > alpha] <- 1
    alpha_in <- 1 - alpha_in
    ## indicate which edges to be kept
    alpha_tot <- alpha_out + alpha_in
    alpha_tot[alpha_tot > 0] <- 1
    M <- M*alpha_tot
    ## assign 0 to NaNs, as they have no incoming edge, outgoing edge or both
    M[is.nan(M)] <- 0
    ## threshold corresponding to alpha
    thres <- min(M[M != 0], na.rm = TRUE)
  }
  if (!weighted){
    M[M > 0] <- 1
  }
  if (!directed){
    M <- M + t(M)
    M <- M - diag(diag(M), n)/2
    warning("The resulting threshold is before merging the directed edges!")
  }
  return(list(AdjacencyMatrix = M, Threshold = thres))
}

##-------------------------------------------------------------------------------------

##########################################
#### Degree and strength distribution ####
##########################################

gen_hist <- function(data = mriot_data, yearlist, mode = c("degree", "strength")) {
  n <- length(yearlist)
  m <- dim(data)[1]
  if (mode == "degree") {
    ## Generate unweighted network with threshold at 0
    data_uw0 <- array(NA, dim = dim(data))
    data_uw0[] <- apply(data, 3, function(M) {
      temp <- netbackbone(M, thres.fixed = 0, weighted = FALSE, directed = TRUE)
      return(temp$AdjacencyMatrix)})
    indeg <- apply(data_uw0, 3, colSums)
    outdeg <- apply(data_uw0, 3, rowSums)
    totdeg <- indeg + outdeg
    colnames(indeg) <- colnames(outdeg) <- colnames(totdeg) <- year_list
    temp_df <- data.frame(year = as.vector(rep(yearlist, each = m, len = n*m*3)),
                          type = as.vector(rep(c("in-degree", "out-degree", "total degree"), 
                                               each = n*m, len = n*m*3)), 
                          degree = as.vector(cbind(indeg[,yearlist], outdeg[,yearlist], totdeg[,yearlist])))
    result <- ggplot(temp_df, aes(x = degree)) + 
      geom_histogram(bins = 30, fill = "white" ,color = "black") + 
      theme(strip.background = element_rect(fill = "grey60", color = "grey60")) + 
      facet_grid(year~type, scales = "free_x")
  } else if (mode == "strength") {
    instr <- apply(data, 3, colSums)
    outstr <- apply(data, 3, rowSums)
    totstr <- instr + outstr
    colnames(instr) <- colnames(outstr) <- colnames(totstr) <- year_list
    temp_df <- data.frame(year = as.vector(rep(yearlist, each = m, len = n*m*3)),
                          type = as.vector(rep(c("in-strength", "out-strength", "total strength"), 
                                               each = n*m, len = n*m*3)), 
                          strength = log(as.vector(cbind(instr[,yearlist], outstr[,yearlist], 
                                                         totstr[,yearlist])) + 1))
    result <- ggplot(temp_df, aes(x = strength)) + 
      geom_histogram(binwidth = 0.5, fill = "white" ,color = "black") + 
      theme(strip.background = element_rect(fill = "grey60", color = "grey60")) + 
      facet_grid(year~type)
  }
  return(result)
}

##-------------------------------------------------------------------------------------

#######################################
########## tail distribution ##########
#######################################

tail_dist <- function(adj, 
                      type = c("out-strength", 
                               "in-strength", 
                               "total-strength", 
                               "out-degree", 
                               "in-degree", 
                               "total-degree"))
  {
  ## load package
  require(ggplot2)
  require(poweRlaw)
  
  ## build unweighted adjacency matrix
  unweight.adj <- adj
  unweight.adj[unweight.adj > 0] <- 1
  
  if (type == "out-strength"){
    temp <- rowSums(adj)
  } else if (type == "in-strength"){
    temp <- colSums(adj)
  } else if (type == "total-strength"){
    temp <- rowSums(adj) + colSums(adj)
  } else if (type == "out-degree"){
    temp <- rowSums(unweight.adj)
  } else if (type == "in-degree"){
    temp <- colSums(unweight.adj)
  } else {
    temp <- rowSums(unweight.adj) + colSums(unweight.adj)
  }

  ## test power-law
  ## continuous for strength; discrete for degree
  if ((type == "out-strength") | (type == "in-strength") | (type == "total-strength")){
    temp_input <- temp + 1e-4
    temp_max <- max(temp_input) + 1
    temp_pl <- conpl$new(temp_input)
  } else {
    temp_input <- temp + 1
    temp_max <- max(temp_input)
    temp_pl <- displ$new(temp_input)
  }
  est_xmin <- estimate_xmin(temp_pl, xmax = temp_max)
  temp_pl$xmin <- est_xmin
  bs_pvalue <- bootstrap_p(temp_pl, xmax = temp_max)
  
  output <- list("res" = temp, "pvalue" = bs_pvalue$p, "alpha" = est_xmin$pars, "xmin" = est_xmin$xmin)
  
  return(output)
}

##-------------------------------------------------------------------------------------

#######################################
########## weighted PageRank ##########
#######################################

## function for computing weighted PageRank
## gamma: damping factor, default value = 0.85
## theta: tuning parameter, default value = 1
## prior.info: node-specific prior information
## default setting = uniformly over all the nodes
wpr <- function(adj, 
                gamma = 0.85, 
                theta = 1, 
                prior.info = rep(1/dim(adj)[1],dim(adj)[1])
)
{
  ## load package
  require(rARPACK)
  
  ## regularity conditions
  if (dim(adj)[1]!=dim(adj)[2]){
    stop("The adjacency matrix is not a square matrix!")
  }
  if ((gamma < 0) | (gamma > 1)){
    stop("The damping factor is not between 0 and 1!")
  }
  if ((theta < 0) | (theta > 1)){
    stop("The tuning parameter is not between 0 and 1!")
  } 
  if (length(prior.info) != dim(adj)[1]){
    stop("The dimension of the prior information is incorrectl!")
  }
  if ((sum(prior.info) == 0) | any(prior.info < 0)){
    stop("The prior information is invalid!")
  }
  if (sum(prior.info) != 1){
    prior.info <- prior.info/sum(prior.info)
    warning("The prior information is not normalized!")
  }
  
  ## get the unweighted adjacency matrix
  unweight.adj <- adj
  unweight.adj[unweight.adj > 0] <- 1
  
  ## construct M and M.star matrix
  n <- dim(adj)[1]
  sink.node <- which(rowSums(adj) == 0)
  M <- theta*t(adj/rowSums(adj)) + (1 - theta)*t(unweight.adj/(rowSums(unweight.adj)))
  M[, sink.node] <- prior.info
  B <- matrix(rep(prior.info, n), nrow = n, ncol = n)
  M.star <- gamma*M + (1 - gamma)*B
  
  ## rARPACK cannot solve solve matrices of 2-by-2
  if (dim(adj)[1] == 2){
    eig_sol <- eigen(M.star)
    eigen_v <- eig_sol$vectors[,1]
    eigen_vstd <- abs(eigen_v)/sum(abs(eigen_v))
    name_v <- c(1:n)
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","WPR")
    return(myres)
  }
  
  ## use rARPACK to solve large-scale matrix
  if (dim(adj)[1] > 2){
    eig_sol <-eigs(M.star, k = 1, which = "LM", mattype = "matrix")
    eigen_v <- Re(eig_sol$vectors)
    eigen_vstd <- abs(eigen_v) / sum(abs(eigen_v))
    name_v <- c(1:n)
    myres <- cbind(name_v, eigen_vstd)
    colnames(myres) <- c("vertex","WPR")
    return(myres)
  }
}
