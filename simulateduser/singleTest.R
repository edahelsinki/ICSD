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
nodes <- as.integer(args[1])            # Arguments used: 5 10
unknown <- as.double(args[2])           # Arguments used: 0 0.33
k <- as.double(args[3])                 # Arguments used: 0.33 0.34 0.37 0.4 0.5
n <- as.integer(args[4])                # Arguments used: 50 100000

if (k == 0.33) k <- 1/3 # Uninformative prior
n.graphs <- 100
noise <- 0.01

all.results <- c()

sapply(seq(n.graphs), function (i) {
  true.process <- create.random.graph(n.nodes=nodes, n.sample=n, noise.min=noise, noise.max=noise)
  all.results <<- rbind(all.results, c(nodes=nodes, n=n, noise=noise, simulate.compare.prior(true.process$data, k, true.process$true.graph, unknown)$results))
})

# Save results
saveRDS(all.results, sprintf('single/nodes%d.u%.2f.k%.2f.n%d.Rds', nodes, unknown, k, n))
