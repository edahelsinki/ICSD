#' THIS SOURCE CODE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND, AND ITS
#' AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH (JMLR) AND JMLR’S
#' PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES, INCLUDING BUT
#' NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#' PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES
#' ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE
#' AUTHOR NOR JMLR, NOR JMLR’S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR
#' DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality
#' of the foregoing, neither the author, nor JMLR, nor JMLR’s publishers and
#' distributors, warrant that the Source Code will be error-free, will operate
#' without interruption, or will meet the needs of the user.

#' Extract a list of parents for each variable in a graph.
#'
#' @param adj.matrix A numerical adjacency matrix representing the graph.
#' @return A named list with an entry for each variable containing a vector of
#'   the names of the variable's parents in the graph.
get.parents <- function (adj.matrix) {
  parents <- list()
  varnames <- colnames(adj.matrix)

  sapply(seq_along(colnames(adj.matrix)), function (column) parents[varnames[column]] <<- list(varnames[which(adj.matrix[, column] == 1)]))

  parents[tryCatch(topo_sort(graph_from_adjacency_matrix(adj.matrix)), warning=function (w) seq_along(colnames(adj.matrix)))]
}


#' Compute the adjusted r-squared of a given data vector and predictions.
#'
#' @param data A vector of the true values.
#' @param preds A vector of the predicted values.
adj.r.squared <- function (data, preds, n, p) {
  ss.res <- sum((data - preds)^2)
  ss.tot <- sum((data - mean(data))^2)

  r2 <- 1 - (ss.res / ss.tot)
  1 - (1 - r2) * (n - 1) / (n - p - 1)
}


#' Return the mean of adjusted r-squared values of a given model represented by
#' an adjacency matrix using assumptions of Gaussianity and linearity.
#'
#' @param data A data frame containing both the training and validation data.
#' @param training.idx A vector of row indices of the training data.
#' @param adj.matrix A numerical adjacency matrix.
#' @param variables The variables to estimate or {NULL} if all variables are
#'   included.
#' @param train A boolean indicating whether the mean r-squared should be
#'   computed for the training data, too.
#' @return Mean of r-squared values of the given model.
calc.mean.R2 <- function (data, training.idx, adj.matrix, variables=NULL, train=F) {
  if (any(adj.matrix == 1 & t(adj.matrix) == 1)) adj.matrix <- get.adj.matrix(pdag2dag(getGraph(adj.matrix))$graph)

  data <- as.data.frame(scale(data, center=F))
  parents <- get.parents(adj.matrix)
  vars <- if (is.null(variables)) names(parents) else variables

  r2s <- as.data.frame(sapply(vars, function (column) {
    if (length(parents[[column]]) == 0) c(validation=0, train=0)
    else {
      fit.lm <- lm(as.formula(paste(column, '~', paste(unlist(parents[column]), collapse = '+'))), data=data[training.idx, ])
      c(validation=adj.r.squared(data[-training.idx, column],  predict.lm(fit.lm, newdata=data[-training.idx, ]),
                                 nrow(data) - length(training.idx), length(parents[[column]])),
        train=summary(fit.lm)$adj.r)
    }
  }))

  fact <- ncol(data) / (ncol(data) - 1)

  if (is.null(variables) && train) round(rowMeans(r2s) * fact, 5)
  else if (is.null(variables)) round(mean(unlist(r2s['validation', ])) * fact, 5)
  else unlist(r2s['validation', , drop=F]) * fact
}


#' Return an estimate of the log-likelihood of a given model represented by an
#' adjacency matrix using assumptions of Gaussianity and linearity.
#'
#' The estimation for a single variable is computed as the mean squared error of
#' a linear regression fit on the variable's parents, divided by 2. For a root
#' variable, i.e. a variable with no parents, the training data mean of the
#' variable is used as an estimate for computing the error. Finally, the
#' estimates are subtracted from each other.
#'
#' @param data A data frame containing both the training and validation data.
#' @param training.idx A vector of row indices of the training data.
#' @param adj.matrix A numerical adjacency matrix.
#' @param variables The variables to estimate.
#' @return An estimate of the log-likelihood of the given model.
calc.loglik <- function (data, training.idx, adj.matrix, variables=NULL) {
  if (any(adj.matrix == 1 & t(adj.matrix) == 1)) adj.matrix <- get.adj.matrix(pdag2dag(getGraph(adj.matrix))$graph)

  data <- as.data.frame(scale(data, center=F))
  parents <- get.parents(adj.matrix)
  vars <- if (is.null(variables)) names(parents) else variables

  log.liks <- sapply(vars, function (column) {
    if (length(parents[[column]]) == 0) -mean((data[-training.idx, column] - mean(data[training.idx, column]))**2) / 2
    else {
      fit.lm <- lm(as.formula(paste(column, '~', paste(unlist(parents[column]), collapse = '+'))), data=data[training.idx, ])
      -mean((data[-training.idx, column] - predict.lm(fit.lm, newdata=data[-training.idx, ]))**2) / 2
    }
  })
  
  if (is.null(variables)) round(sum(log.liks), 5) else setNames(as.list(log.liks), variables)
}


#' Computes a set of changes that would improve the score of the model by
#' testing additions, deletions, and reversals of edges in the model.
#'
#' @param data A data frame containing both the training and validation data.
#' @param training.idx A vector of row indices of the training data.
#' @param adj.matrix A numerical adjacency matrix.
#' @param calc.score A function for computing the score of a model.
#' @param all A boolean indicating whether all edits whether improvements or not
#'   are included.
#' @param ucert.threshold If a reversal of an edge would change the score at
#'   most by this value, the change is determined uncertain, regardless of
#'   whether the change would be positive or negative. By default 0.001.
#' @return A named list with components `Add`, `Del`, `Rev`, and `Ucert`, each
#'   containing a matrix of improvements, defined by the row and column indices
#'   and the change in score.
calc.improvements <- function (data, training.idx, adj.m, calc.score, all=F, ucert.threshold=0.001) {
  varnames <- colnames(adj.m)
  curr.score <- calc.score(data, training.idx, adj.m, varnames)
  add.scores <- del.scores <- rev.scores <- ucert.scores <- matrix(0, nrow(adj.m), ncol(adj.m))

  missing.edges <- which(adj.m == 0 & t(adj.m) == 0 & row(adj.m) != col(adj.m), arr.ind=T)
  if (nrow(missing.edges) > 0)
    apply(missing.edges, 1, function (cell) {
      adj.m[cell[1], cell[2]] <- 1
      score.diff <- round(calc.score(data, training.idx, adj.m, varnames[cell[2]])[[varnames[cell[2]]]] - curr.score[[varnames[cell[2]]]], 5)
      adj.m[cell[1], cell[2]] <- 0

      if (score.diff > 0 || all) add.scores[cell[1], cell[2]] <<- score.diff
    })

  dir.edges <- which(adj.m == 1 & t(adj.m) == 0, arr.ind=T)
  if (nrow(dir.edges) > 0)
    apply(dir.edges, 1, function(cell) {
      adj.m[cell[1], cell[2]] <- 0
      del.score.diff <- round(calc.score(data, training.idx, adj.m, varnames[cell[2]])[[varnames[cell[2]]]] - curr.score[[varnames[cell[2]]]], 5)
  
      adj.m[cell[2], cell[1]] <- 1
      mod.score <- calc.score(data, training.idx, adj.m, varnames[c(cell[1], cell[2])])
      rev.score.diff <- round(mod.score[[varnames[cell[1]]]] + mod.score[[varnames[cell[2]]]] - curr.score[[varnames[cell[1]]]] - curr.score[[varnames[cell[2]]]], 5)
      adj.m[cell[1], cell[2]] <- 1
      adj.m[cell[2], cell[1]] <- 0

      add.del <- del.score.diff > 0 && del.score.diff > rev.score.diff

      if (add.del || all) del.scores[cell[1], cell[2]] <<- del.score.diff
      if ((!add.del && rev.score.diff > ucert.threshold) || all) rev.scores[cell[1], cell[2]] <<- rev.score.diff
      else if (!add.del && rev.score.diff != 0 && abs(rev.score.diff) <= ucert.threshold && !all) ucert.scores[cell[1], cell[2]] <<- rev.score.diff
    })

  list(Add=add.scores, Del=del.scores, Rev=rev.scores, Ucert=ucert.scores)
}


#' Check a graph represented by an adjacency matrix for cycles.
#'
#' Cycles are found using depth-first-search of the graph.
#'
#' @param adj.m An adjacency matrix.
#' @return Logical indicating whether the graph is acyclic.
is.acyclic <- function (adj.m) {
  visited <<- recursion.stack <<- rep(F, ncol(adj.m))

  is.cyclic.util <- function (column) {
    visited[column] <<- recursion.stack[column] <<- T

    for (neighbour in which(adj.m[column, ] == 1))
      if ((!visited[neighbour] && is.cyclic.util(neighbour)) || recursion.stack[neighbour]) return(T)

    recursion.stack[column] <<- F
    F
  }

  for (column in seq_len(ncol(adj.m))) {
    if (!visited[column] && is.cyclic.util(column)) return(F)
  }

  T
}


#' Count the number of v-structures in a (partially) oriented graph.
#'
#' @param adj.matrix A numerical adjacency matrix.
count.v.structures <- function (orig.adj.matrix) {
  adj.matrix <- orig.adj.matrix
  adj.matrix[which(adj.matrix == t(adj.matrix))] <- 0
  possible.colliders <- which(apply(adj.matrix, 2, sum) > 1)
  if (length(possible.colliders) == 0) return(0)

  sum(apply(adj.matrix[, possible.colliders, drop=F], 2, function (column) {
    possible.nodes <- combn(which(column == 1), 2)
    sum(apply(possible.nodes, 2, function (idx) if (orig.adj.matrix[idx[1], idx[2]] == 1 || orig.adj.matrix[idx[2], idx[1]] == 1) 0 else 1))
  }))
}


#' Return given graph object as an adjacency matrix.
#' 
#' @param graph.nel A graph object, either graphNEL or graphAM.
#' @return A numerical adjacency matrix.
get.adj.matrix <- function (graph) {
  graph.nel <- as(graph, 'graphNEL')
  nodes <- graph.nel@nodes
  edges <- graph.nel@edgeL

  adj.m <- matrix(0, length(nodes), length(nodes), dimnames=list(nodes, nodes))
  sapply(nodes, function (node) if (length(edges[[node]]$edges) > 0) adj.m[node, edges[[node]]$edges] <<- 1)

  adj.m
}