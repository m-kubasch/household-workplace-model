# household-workplace-model


The content of this repository is associated to the preprint "Large population limit for a multilayer *SIR* model including households and workplaces" (M. Kubasch, 2023 - available on [ArXiv](https://arxiv.org/abs/2305.17064)), whose abstract yields some context:

> We study a multilayer *SIR* model with two levels of mixing, namely a global level which is uniformly mixing, and a local level with two layers distinguishing household and workplace contacts, respectively. We establish the large population convergence of the corresponding stochastic process. For this purpose, we use an individual-based  model whose state space explicitly takes into account the duration of infectious periods. This allows to deal with the  natural correlation of the epidemic states of individuals whose household and workplace share a common infected. In a general setting where a non-exponential distribution of infectious periods may be considered, convergence to the unique deterministic solution of a measure-valued equation is obtained. In the particular case of exponentially distributed infectious periods, we show that it is possible to further reduce the obtained deterministic limit, leading to a closed, finite dimensional dynamical system capturing the epidemic dynamics. This model reduction subsequently is studied from a numerical point of view. We illustrate that the dynamical system derived from the large population approximation is a pertinent model reduction when compared to simulations of the stochastic process or to an alternative edge-based compartmental model, both in terms of accuracy and computational cost.

It contains both the code developed for this project, as well as the simulated data allowing to reproduce the results. 



## Contents

### *code*

1. **Modules:**

	- `agentbasedmodel`: simulations of the stochastic agent-based model (ABM) using Gillespie's algorithm. 
	- `largepoplimit`: dynamical system arising as the large population limit of the stochastic model.
	- `ebcm`: alternative dynamical system in the family of edge-based compartmental models (EBCM).
	- `comptime`: assessing computation time of either stochastic simulations or numerical solution of the dynamical systems.
	- `ioput`: everything connected to input/output.
	- `graphics`: everything connected to visualizing the data.


2. **Mains:** the following scripts give examples of how to use these modules.

	- `main-simu.py`: simulation of one epidemic trajectory given by the stochastic ABM. The output is only saved if the proportion of infected reaches at least 0.5%, otherwise the trajectroy is discarded (this can be changed by using the *threshold* parameter of `ioput.save_simulation`). The script is designed to be launched using the command line.
	- `main-reduction.py`: numerical solution of the dynamic system arising as the large population limit of the stochastic model, for all parameter sets considered in the article.
	- `main-ebcm.py`: numerical solution of the EBCM. The script is designed to be launched using the command line.
	- `main-ctime.py`: script allowing to time the execution of either one simulation of the stochastic ABM, or the numerical solution of either of the two dynamical systems (large population limit, EBCM).
	- `main-simu-initialcdt.py`: simulation of an epidemic trajectory of the ABM starting from a single infected, up to the point when the proportion of infected reaches a given threshold. If this point is reached, the observed state of the variables of the large population limit is saved. The script is designed to be launched using the command line.
	- `main-infere-initialcdt.py`: numerical solution of the large population limit dynamical system, using an initial condition infered from simulation data.
	- `main-figure2.py`, `main-figure3.py`, `main-figureA1.py`, `main-figureA2.py`, `main-figureSM1.py`, `main-tableSM2.py`: scripts allowing to reproduce the corresponding figures and table of the article. 

### *article-data*

This folder contains the data necessary for reproducing the results presented in the article. All files were generated using the code contained in this repository. A brief overview of the contents:
- *ctime*: computation times (ABM, large population limit, EBCM) obtained for different parameter sets. The data is organized by proportions of infections in the general population and reproduction number, which are enough to distinguish the scenarios considered in the article.
- *ebcm*: outputs of the EBCM.
- *figures*: the title says it all.
- *limit*: outputs of the large population limit dynamical system.
- *simulations*: outputs of the ABM. Roughly organized by main application of the data (comparison of the ABM to the large population limit or to the EBCM, inference of an appropriate initial condition for the large population limit from simulated data).
 

## License

The code content is licensed under the [GNU GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html#license-text), whereas data and figures are licensed under the [CC-BY-4.0 license](https://creativecommons.org/licenses/by/4.0/).


## Contact

Madeleine Kubasch -- madeleine.kubasch@polytechnique.edu
