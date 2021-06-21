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

create.random.graph <- function (n.nodes=4, n.sample=1000, edge.prob=.3, noise.min=.05, noise.max=2, rchoices=c('rnorm', 'runif')) {
  nodes <- LETTERS[1:n.nodes]
  edges <- rbinom(sum(seq(n.nodes - 1)), 1, edge.prob)
  rfuncs <- sample(rchoices, n.nodes, replace=T)
  params <- runif(n.nodes, min=noise.min, max=noise.max)

  adj.m <- matrix(0, n.nodes, n.nodes, dimnames=list(nodes, nodes))

  curr.idx <- 1
  curr.size <- n.nodes - 1
  pars <- setNames(lapply(seq_along(nodes), function (i) {
    if (curr.size == 0) return(c())

    poss.edges <- edges[curr.idx:(curr.idx + curr.size - 1)]

    curr.idx <<- curr.idx + curr.size
    curr.size <<- curr.size - 1

    adj.m[which(poss.edges == 1), n.nodes + 1 - i] <<- 1
    LETTERS[which(poss.edges == 1)]
  }), nodes[n.nodes:1])

  relations <- setNames(lapply(seq_along(nodes), function (node) { 
    if (length(pars[[nodes[node]]]) == 0) list(rfunc=rfuncs[node], param=params[node])
    else list(rfunc=rfuncs[node], param=params[node], par=pars[[nodes[node]]], coef=runif(length(pars[[nodes[node]]]), -3, 3))
  }), nodes)

  data <- create.linear.mixed.data(n.sample, relations)

  return(list(data=data, true.graph=adj.m))
}


calc.shd <- function (model, true) 
  return(sum(model == 1 & true == 0 & t(true) == 0)
         + sum(model == 0 & t(model) == 0 & (true == 1 | t(true) == 1)) / 2
         + sum(model == 1 & t(true) == 1))


compare.graphs <- function (target, true, print=T) {
  tp <- sum(target == 1 & true == 1)
  fp <- sum(target == 1 & true == 0)
  tn <- sum(target == 0 & true == 0) - ncol(target)
  fn <- sum(target == 0 & true == 1)

  exact <- c(tp=tp, fp=fp, tn=tn, fn=fn)

  tp <- sum(target == 1 & (true == 1 | t(true) == 1))
  fp <- sum(target == 1 & true == 0 & t(true) == 0)
  tn <- (sum(target == 0 & t(target) == 0 & true == 0 & t(true) == 0) - ncol(target)) / 2
  fn <- sum(target == 0 & t(target) == 0 & (true == 1 | t(true) == 1)) / 2

  skeleton <- c(tp=tp, fp=fp, tn=tn, fn=fn)

  incorrect.orientation <- sum(target == 1 & t(true) == 1)

  training.idx <- 1:(nrow(current$data) * current$training.ratio)
  target.score <- current$calc.score(current$data, training.idx, target)
  true.score <- current$calc.score(current$data, training.idx, true)

  c(exact=exact, skeleton=skeleton, target.score=target.score, true.score=true.score, shd=fn + fp + incorrect.orientation)
}


calc.distances <- function (graphs)
  unlist(apply(combn(names(graphs), 2), 2, function (colm)
    setNames(list(calc.shd(graphs[[colm[1]]], graphs[[colm[2]]])), sprintf('%s-%s', colm[1], colm[2]))))
