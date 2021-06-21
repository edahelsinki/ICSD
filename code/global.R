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

library(igraph)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(pcalg)

source('graphComp.R')
source('csdAlgorithms.R')


#' Initialize the session state.
#'
#' Runs the default set of algorithms on the given data set, selecting the one
#' with the highest score as the edit matrix.
#'
#' @param data A data frame.
#' @return The initialized state as a named list:\n data - A data frame
#'   containing the specified data set. \n algorithms - A list of algorithms,
#'   the result of [get.initial.alg()]. \n models - A matrix, result of
#'   [run.all.algorithms()]. \n calc.score - The function [calc.mean.R2()]. \n
#'   edit.mtx - A numeric adjacency matrix, chosen to be the adjacency matrix of
#'   the model with highest score. \n training.ratio - 0.5.
initcurrent <- function (data=NULL) {
  algorithms <- get.initial.alg(data[(1:(nrow(data) * 0.5)), ])
  models <- run.all.algorithms(algorithms, data, training.ratio=0.5)

  list(
    data=data,
    algorithms=algorithms,
    models=models,
    calc.score=calc.mean.R2,
    edit.mtx=models[[1, 'adj.m']],
    training.ratio=0.5
  )
}


#' Run selected algorithms on the current data, using the selected training
#' ratio.
#'
#' Updates the current edit and summary matrices according to the results.
#'
#' @param scoring.metric A string, name of the scoring function to use in the
#'   evaluation of the results.
#' @return The character string 'Run selected algorithms' to add to the action
#'   log.
run.algorithms <- function (scoring.metric) {
  current$calc.score <<- eval(parse(text=scoring.metric))
  current$models <<- run.all.algorithms(current$algorithms, current$data, current$training.ratio, current$calc.score)
  current$edit.mtx <<- current$models[[1, 'adj.m']]

  'Run selected algorithms'
}


#' Return models.
get.models <- function () current$models


#' Add an edge from `y` to `x` to the current edit matrix.
#'
#' @param x The end node of the edge.
#' @param y The start node of the edge.
#' @return A character string describing the action to add to the action log.
add.edge <- function (x, y) {
  current$edit.mtx[y, x] <<- 1
  cols <- colnames(current$data)

  sprintf('Add edge %s -> %s', cols[y], cols[x])
}


#' Reverse the edge from `y` to `x` in the current edit matrix.
#'
#' @param x The end node of the edge.
#' @param y The start node of the edge.
#' @return A character string describing the action to add to the action log.
reverse.edge <- function (x, y) {
  current$edit.mtx[y, x] <<- 0
  current$edit.mtx[x, y] <<- 1
  cols <- colnames(current$data)

  sprintf('Reverse edge to %s -> %s', cols[x], cols[y])
}


#' Delete the edge from `x` to `y` from the current edit matrix.
#'
#' @param y The start node of the edge.
#' @param x The end node of the edge.
#' @return A character string describing the action to add to the action log.
delete.edge <- function (x, y) {
  current$edit.mtx[x, y] <<- 0
  cols <- colnames(current$data)

  sprintf('Remove edge %s -> %s', cols[x], cols[y])
}


#' Set the selected model as the current edit matrix.
#'
#' @param name A character string, name of the selected model.
#' @return A character string describing the action performed to add to the
#'   action log.
change.model <- function (name) {
  models <- get.models()
  current$edit.mtx <<- models[[which(models[, 'name'] == name), 'adj.m']]
  sprintf('Select model %s', name)
}


options(stringsAsFactors = FALSE)