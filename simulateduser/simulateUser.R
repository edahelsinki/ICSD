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

setwd('../code/')
source('global.R')
setwd('../simulateduser/')

source('randomGraphs.R')


find.map <- function (calc.posterior) {
  map <- list(p=calc.posterior(), type=NULL, cell=NULL)

  missing.edges <- which(current$edit.mtx == t(current$edit.mtx) & row(current$edit.mtx) != col(current$edit.mtx), arr.ind=T)
  if (nrow(missing.edges) > 0)
    apply(missing.edges, 1, function (cell) {
      current$edit.mtx[cell[1], cell[2]] <<- 1
      if (calc.posterior() > map$p) map <<- list(p=calc.posterior(), type='Add', cell=cell)
      current$edit.mtx[cell[1], cell[2]] <<- 0
    })

  dir.edges <- which(current$edit.mtx == 1 & t(current$edit.mtx) == 0, arr.ind=T)
  if (nrow(dir.edges) > 0)
    apply(dir.edges, 1, function (cell) {
      current$edit.mtx[cell[1], cell[2]] <<- 0
      if (calc.posterior() > map$p) map <<- list(p=calc.posterior(), type='Del', cell=cell)
      current$edit.mtx[cell[2], cell[1]] <<- 1
      if (calc.posterior() > map$p) map <<- list(p=calc.posterior(), type='Rev', cell=cell)
      current$edit.mtx[cell[1], cell[2]] <<- 1
      current$edit.mtx[cell[2], cell[1]] <<- 0
    })

  map
}


simulate.user.prior <- function(k, true.graph, priors) {

  training.idx <- 1:(nrow(current$data) * current$training.ratio)
  
  get.probs <- function (model=current$edit.mtx)
    c(priors[model == 1], 1 - (priors + t(priors))[model == t(model) & col(priors) > row(priors)])
  
  calc.prior <- function (model=current$edit.mtx) sum(log(get.probs(model)))
  calc.posterior <- function (model=current$edit.mtx) calc.prior(model) + calc.loglik(current$data, training.idx, model)
  
  start.score <- current$calc.score(current$data, training.idx, current$edit.mtx)
  start.shd <- calc.shd(current$edit.mtx, true.graph)
  
  while (T) {
    map <- find.map(calc.posterior)

    # No edit would improve the posterior
    if (is.null(map$type)) break

    # Perform the improvement action
    if (map$type == 'Rev') reverse.edge(map$cell[2], map$cell[1])
    else if (map$type == 'Add') add.edge(map$cell[2], map$cell[1])
    else if (map$type == 'Del') delete.edge(map$cell[1], map$cell[2])
  }
  
  return(list(start.score=start.score, start.shd=start.shd))
}


simulate.compare.prior <- function (data, k, true.graph, ratio.unknown, multi=F) {
  priors <- true.graph * k
  priors[priors == 0] <- (1 - k) / 2
  diag(priors) <- 0
  n.nodes <- ncol(true.graph)
  known.edges <- which(col(priors) > row(priors))

  if (ratio.unknown > 0) {
    unknown.edges <- sample(which(col(priors) > row(priors)), ratio.unknown * (n.nodes^2 - n.nodes) / 2)
    known.edges <- known.edges[which(!known.edges %in% unknown.edges)]
    priors[unknown.edges] <- 1/3
    priors[t(priors) == 1/3] <- 1/3
  }

  current <<- initcurrent(data=data)
  if (multi) {
    final.graphs <- list()

    # Start from an empty graph
    current$edit.mtx <<- matrix(0, nrow(true.graph), ncol(true.graph), dimnames=list(rownames(true.graph), colnames(true.graph)))
    res <- simulate.user.prior(k, true.graph, priors)
    results <- c(k=k, unknown=ratio.unknown, compare.graphs(current$edit.mtx, true.graph, print=F), edits=nrow(current$edits), unlist(res))
    final.graphs$empty <- current$edit.mtx

    # Start from the true graph
    current$edit.mtx <<- true.graph
    res <- simulate.user.prior(k, true.graph, priors)
    results <- rbind(results, c(k=k, unknown=ratio.unknown, compare.graphs(current$edit.mtx, true.graph, print=F), edits=nrow(current$edits), unlist(res)))
    final.graphs$true <- current$edit.mtx

    names <- unlist(current$models[, 'name'])
    short.names <- c('GES', 'PC.0.01', 'PC.0.1', 'LiNGAM')
    best <- sapply(short.names, function (str) which(startsWith(names, str))[1])

    # Start from the top result for each of the default algorithms
    for (i in which(!is.na(best))) {
      change.model(names[best[i]])
      res <- simulate.user.prior(k, true.graph, priors)
      results <- rbind(results, c(k=k, unknown=ratio.unknown, compare.graphs(current$edit.mtx, true.graph, print=F), edits=nrow(current$edits), unlist(res)))
      final.graphs[[short.names[i]]] <- current$edit.mtx
    }

    rownames(results) <- c('empty', 'true', names(best[!is.na(best)]))
    return(list(results=results, distance.mtx=calc.distances(final.graphs)))
  }

  results <- simulate.user.prior(k, true.graph, priors)
  return(list(results=c(k=k, unknown=ratio.unknown, compare.graphs(current$edit.mtx, true.graph, print=F), edits=nrow(current$edits), unlist(results))))
}

