# Experiments for _Interactive Causal Structure Discovery in Earth System Sciences_

This directory contains the code needed to run the experiments and use cases in the paper
> **Laila Melkas, Rafael Savvides, Suyog Chandramouli, Jarmo Mäkelä, Tuomo Nieminen, Ivan Mammarella, and Kai Puolamäki**
> *Interactive Causal Structure Discovery in Earth System Sciences*.
> In ACM SIGKDD Workshop on Causal Discovery, 2021, PMLR.


#### Requirements
* R
* R packages: ggplot2, igraph, RColorBrewer, data.table, pcalg
* Python 3
* Python packages: numpy, pandas

#### Data
FLUXNET2015 Hyytiälä, measurements between January 2013 and December 2015, inclusive. The data is licensed under a Creative Commons Attribution 4.0 International License.

**Reference for the data**: Ivan Mammarella. Drought 2018 Fluxdata Preview Selection, Hyytiälä, 1995-12-31–2018-12-31, 2020. URL: https://hdl.handle.net/11676/EBmVEuoJaOmOw8QmUyyh6G-n.

## Description of folders

Folder `code` contains the code base for applying the suggested approach to perform interactive causal discovery on a given data set.
The code used to process the full FLUXNET Hyytiälä data set and files containing the processed data are located in the `data` folder.
Use cases presented in the paper can be rerun by running the Rmarkdown file `KDD_use_cases.Rmd` in folder `usecases`.
Scripts for running the experiments with synthetic data and simulated user are located in the folder `simulateduser`.
All of the figures created in the experiments and use cases are saved in folder `figures`.
Below is a more detailed description of how to run the simulated user experiments.

## Synthetic data and simulated user

All of the files referenced below are located in the folder `simulateduser`.
Folders `single`, `multi`, and `distance` already contain the results from executed simulations.
The simulations can be re-executed by running script `singleTest.R` or `multiTest.R` depending on the experiment with appropriate parameters.
Code to inspect and plot the results is located in the file `KDD_results.Rmd`.
Below is a more detailed list of files and their contents within the folder.

#### Structure
* `randomGraphs.R`: Code to create random graphs and compare two graphs against each other in terms of structural Hamming distance, differences in structure, and scores.
* `simulateUser.R`: Code to simulate a user with parameters `k` for confidence in knowledge and `unknown` for ratio of possible edges with flat priors.
* `singleTest.R`: Code to run Experiment 1. Results are saved in folder `single`. Parameters to pass when running the script consist of number of nodes (5 and 10), sample size (50 and 100000), ratio of unknown information (0 and 0.33), and level of confidence in known information (0.33, 0.34, 0.37, 0.4, 0.5) – parameters used in the paper in parentheses. **NB:** the script may take a couple of days to run with large sample size and uniform prior (`k = 0.33`).
* `multiTest.R`: Code to run Experiment 2. Results similar in format to Experiment 1 are saved in folder `multi` and pairwise distances between the final models when starting from different initial models are stored in folder `distance`. When running the script, the parameter `k = [0.33, 0.34, 0.37, 0.4, 0.5]` is passed to control the level of confidence in known information.
* `KDD_results.Rmd`: RMarkdown file to plot and print results. Figures created are saved in the folder `figures`.
