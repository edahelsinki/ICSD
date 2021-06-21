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

#' Recurse over all possible combinations of orientations of undirected edges in
#' an adjacency matrix and return a list of all valid DAGS.
#'
#' @param adj.m A numerical adjacency matrix.
#' @param udir.edges A matrix of indices of undirected edges in the adjacency
#'   matrix.
#' @param n.v.structures A numerical, number of v-structures in the original
#'   CPDAG.
#' @return A list of valid DAGs as adjacency matrices.
adj.matrix2dags <- function (adj.m, udir.edges, n.v.structures) {
  dags <- NULL
  if (length(udir.edges) == 0 && is.acyclic(adj.m)) dags <- list(adj.m)
  else if (length(udir.edges) != 0) {
    adj.m[udir.edges[1, 1], udir.edges[1, 2]] <- 0
    if (count.v.structures(adj.m) == n.v.structures)
      dags <- rbind(dags, adj.matrix2dags(adj.m, udir.edges[-1, , drop=F], n.v.structures))

    adj.m[udir.edges[1, 1], udir.edges[1, 2]] <- 1
    adj.m[udir.edges[1, 2], udir.edges[1, 1]] <- 0
    if (count.v.structures(adj.m) == n.v.structures)
      dags <- rbind(dags, adj.matrix2dags(adj.m, udir.edges[-1, , drop=F], n.v.structures))
  }

  if (is.null(dags)) return(dags)
  as.matrix(dags)
}


#' A set of functions to initialize an algorithm with the set of given
#' arguments.
#'
#' @param algfunction Function to call to initialize the algorithm.
#' @param ... Parameters to pass to specific algorithms.
#' @seealso [pcalg], [fcialg], [gesalg], [kpcalg], [fhcalg], [lingamalg].
alg <- function (algfunction, ...) algfunction(...)


#' Call the correct method to run the given algorithm based on its class.
#'
#' @param alg A list structure whose class determines the correct method to
#'   call, the result of one of [pcalg], [fcialg], [gesalg], [kpcalg], [fhcalg],
#'   [lingamalg].
#' @param data The data set on which the algorithm is run, both training and
#'   validation data.
#' @param training.idx A vector of row indices of the training data.
#' @param calc.score A function used to compute scores for each model.
#' @return A matrix where each row represents one model. Columns are\n `name` -
#'   a unique identifier of the model,\n `score` - the score of the model,\n
#'   `valid` - The logical TRUE,\n `adj.m` - adjacency matrix representation of the models
#'   underlying DAG.
#' @seealso [pcalg::pc], [pcalg::fci], [pcalg::ges], [kpcalg::kpc], [SELF::fhc],
#'   [pcalg::lingam]
run <- function (alg, data, training.idx, calc.score) UseMethod('run')


#' Initialize a PC algorithm.
#'
#' @param vars A numerical vector containing the indices of variables in the
#'   data set to include when running the algorithm.
#' @param alpha A number in the range `(0,1)`, determines the significance level
#'   used for performing conditional independence tests.
#' @param suffStat A list with two named elements, `C` for the correlation
#'   matrix and `n` for the sample size.
#' @param correlation A character string, name of the correlation coefficient
#'   used.
#' @return A list structure of class `pcalg` which contains the given
#'   parameters.
pcalg <- function (vars, alpha, suffStat, correlation) {
  structure(list(alpha=alpha, suffStat=suffStat, vars=vars, name=paste0('PC.', alpha, '.', correlation, '(', paste(vars, collapse=','), ')')), class='pcalg')
}


#' @rdname run
run.pcalg <- function (alg, data, training.idx, calc.score) {
  cat('Running PC...\n')
  varnames <- colnames(data)[alg$vars]

  pc.fit <<- pc(
    alg$suffStat, indepTest=gaussCItest, u2pd='retry',
    labels=varnames, alpha=alg$alpha, skel.method='stable'
  )

  adj.m <- get.adj.matrix(pc.fit@graph)
  udir.edges <- which(adj.m == 1 & t(adj.m) == 1, arr.ind=T)
  udir.edges <- udir.edges[which(udir.edges[, 'row'] < udir.edges[, 'col']), , drop=F]
  all.dags <- adj.matrix2dags(adj.m, udir.edges, count.v.structures(adj.m))

  if (is.null(all.dags)) return(NULL)
  adj.ms <- matrix(nrow=nrow(all.dags), ncol=5)
  adj.m <- matrix(0, ncol(data), ncol(data), dimnames=list(colnames(data), colnames(data)))

  sapply(1:nrow(adj.ms), function (i) {

    adj.m[varnames, varnames] <- all.dags[i, ][[1]]
    adj.ms[i, ] <<- list(name=paste0(alg$name, '_', i), score=calc.score(data, training.idx, adj.m),
                         valid=T, adj.m=adj.m, cpdag=ifelse(dag2cpdag(adj.m), 1, 0))
    adj.ms <<- matrix(adj.ms, ncol=5, dimnames=list(NULL, c('name', 'score', 'valid', 'adj.m', 'cpdag')))
  })

  if (nrow(adj.ms) > 5) adj.ms <- adj.ms[sort.list(unlist(adj.ms[, 'score']), decreasing=T), ]
  adj.ms <- cbind(adj.ms, cpdag.score=rep(calc.score(data, training.idx, get.adj.matrix(pdag2dag(pc.fit@graph)$graph)), nrow(adj.ms)))
  adj.ms[1:(min(5, nrow(adj.ms))), ]
}


#' Initialize a GES algorithm.
#'
#' @param vars A numerical vector containing the indices of variables in the
#'   data set to include when running the algorithm.
#' @param score An instance of a class derived from [pcalg::Score], determining
#'   the score used for the greedy search.
#' @return A list structure of class `gesalg` which contains the given
#'   parameters.
gesalg <- function (vars, score) structure(list(score=score, vars=vars, name=paste0('GES', '(', paste(vars, collapse=','), ')')), class='gesalg')


#' @rdname run
run.gesalg <- function (alg, data, training.idx, calc.score) {
  cat('Running GES...\n')

  ges.fit <- ges(alg$score, labels=colnames(data)[alg$vars], phase=c('forward', 'backward', 'turning', iterate=T))

  ges.graph <- matrix(0, ncol(data), ncol(data), dimnames=list(colnames(data), colnames(data)))
  lapply(ges.fit$essgraph$.nodes, function (node)
    if (length(ges.fit$essgraph$.in.edges[[node]]) > 0) ges.graph[sapply(ges.fit$essgraph$.in.edges[[node]], function (var) alg$vars[var]), node] <<- 1
  )

  udir.edges <- which(ges.graph == 1 & t(ges.graph) == 1, arr.ind=T)
  udir.edges <- udir.edges[which(udir.edges[, 'row'] < udir.edges[, 'col']), , drop=F]

  all.dags <- adj.matrix2dags(ges.graph, udir.edges, count.v.structures(ges.graph))

  if (is.null(all.dags)) return(NULL)

  adj.ms <- matrix(nrow=nrow(all.dags), ncol=5)

  sapply(1:nrow(all.dags), function (i) {
    adj.ms[i, ] <<- list(
      name=paste0(alg$name, '_', i), score=calc.score(data, training.idx, all.dags[i, ][[1]]),
      valid=T, adj.m=all.dags[i, ][[1]], cpdag=ges.graph)
    adj.ms <<- matrix(adj.ms, ncol=5, dimnames=list(NULL, c('name', 'score', 'valid', 'adj.m', 'cpdag')))
  })

  if (nrow(adj.ms) > 5) adj.ms <- adj.ms[sort.list(unlist(adj.ms[, 'score']), decreasing=T), ]
  adj.ms <- cbind(adj.ms, cpdag.score=rep(calc.score(data, training.idx, get.adj.matrix(pdag2dag(getGraph(ges.graph))$graph)), nrow(adj.ms)))
  adj.ms[1:(min(5, nrow(adj.ms))), ]
}


#' Initialize a LiNGAM algorithm.
#'
#' @param vars A numerical vector containing the indices of variables in the
#'   data set to include when running the algorithm.
#' @return A list structure of class `lingamalg` which contains the given
#'   parameters.
lingamalg <- function (vars) structure(list(vars=vars, name=paste0('LiNGAM', '(', paste(vars, collapse=','), ')')), class='lingamalg')


#' @rdname run
run.lingamalg <- function (alg, data, training.idx, calc.score) {
  cat('Running LiNGAM...\n')

  tryCatch({
    lingam.fit <- lingam(data[training.idx, alg$vars])

    adj.m <- matrix(0, ncol(data), ncol(data), dimnames=list(colnames(data), colnames(data)))
    adj.m[colnames(data)[alg$vars], colnames(data)[alg$vars]] <- t(lingam.fit$Bpruned)
    adj.m[adj.m != 0] <- 1

    score <- calc.score(data, training.idx, adj.m)

    list(name=alg$name, score=score, valid=T, adj.m=adj.m,
      cpdag=ifelse(dag2cpdag(adj.m), 1, 0), cpdag.score=score)
  }, error = function (e) {
    cat(paste('ERROR:', e))
    NULL
  })
}


#' Return a hard coded list of initial algorithms.
#'
#' @param data A data frame of the training data on which the algorithms will be
#'   run on.
#' @return A list of algorithms, list structures returned by [alg].
get.initial.alg <- function (data) {
  vars <- 1:ncol(data)
  corr <- cor(data)
  suffStat <- list(C=corr, n=nrow(data))
  gauss.score <- new('GaussL0penObsScore', data=data)

  list(
    alg(pcalg, vars, alpha=0.1, suffStat=suffStat, correlation='pearson'),
    alg(pcalg, vars, alpha=0.01, suffStat=suffStat, correlation='pearson'),
    alg(gesalg, vars, score=gauss.score),
    alg(lingamalg, vars)
  )
}


#' Run all of the given algorithms on training data.
#'
#' @param algorithms A list of algorithms, list structures returned by [alg].
#' @param data Training and validation data as a data frame.
#' @param training.ratio The ratio of the data used for training, 0.5 by
#'   default.
#' @param calc.score The function used to compute scores for the models,
#'   [calc.mean.R2] by default.
#' @return A matrix where each row represents one model. Columns are\n `name` -
#'   a unique identifier of the model,\n `score` - the score of the model,\n
#'   `valid` - The logical TRUE,\n `adj.m` - adjacency matrix representation of the models
#'   underlying DAG.
run.all.algorithms <- function (algorithms, data, training.ratio=0.5, calc.score=calc.mean.R2) {
  training.idx <- 1:(nrow(data) * training.ratio)

  models <- c()
  sapply(algorithms, function (algorithm) models <<- rbind(models, run(algorithm, data, training.idx, calc.score)))

  if (nrow(models) > 1) models <- models[sort.list(unlist(models[, 'score']), decreasing=T), ]

  models
}

