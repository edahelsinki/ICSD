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

source('simulateUser.R')


args <- commandArgs(trailingOnly=TRUE)
k <- as.double(args[1])                 # Arguments used: 0.33 0.34 0.37 0.4 0.5

n.graphs <- 100
n <- 1000
unknown <- 0.33
noise <- 0.01
nodes <- 7

if (k == 0.33) k <- 1/3 # Uninformative prior

all.results <- list()
all.distances <- c()

sapply(seq(n.graphs), function (i) {
  true.process <- create.random.graph(n.nodes=nodes, n.sample=n, noise.min=noise, noise.max=noise)
  res <- simulate.compare.prior(true.process$data, k, true.process$true.graph, unknown, multi=T)

  all.results <<- rbind(all.results, cbind(nodes=nodes, n=n, noise=noise, res$results))
  all.distances <<- rbind(all.distances, c(nodes=nodes, noise=noise, unknown=unknown, n=n, k=k, res$distance.mtx))
})

# Save results
saveRDS(all.results, sprintf('multi/k%.2f.unknown%.2f.noise%.2f.nodes%d.Rds', k, unknown, noise, nodes))

# Save distance results
saveRDS(all.distances, sprintf('distances/k%.2f.unknown%.2f.noise%.2f.nodes%d.Rds', k, unknown, noise, nodes))
