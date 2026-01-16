[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Calcium-based plasticity model and irregular, bursty activity patterns
==============================

Synaptic plasticity—the activity-dependent change in the strength of connections between neurons—is widely believed to underlie learning and memory. It is sensitive not only to the rate of neuronal firing but also to the precise timing of presynaptic and postsynaptic spikes. In experimental settings, plasticity is often studied using regular stimulation patterns, where action potentials are evenly spaced and spike timing is fixed. However, such artificial differs from the irregular and bursty firing patterns observed _in vivo_. In this repository, we investigate how more naturalistic firing patterns influence synaptic plasticity.

The Python code provided here implements the calculations for synaptic strength changes induced by irregular spike pairs. For related implementations and other publications of the calcium-based plasticity model, see [this repository](https://github.com/mgraupe/CalciumBasedPlasticityModel).

These scripts are associated with the following manuscript, currently submitted for publication:

**Yulia Dembitskaya, Silvana Valtcheva, Yihui Cui, Zhiwei Zheng, Srdjan Ostojic, Michael Graupner, and Laurent Venance (2026).** _Irregular and bursty neuronal firing extend spike-timing dependent plasticity window to behavioral timescales._

## Features
- The code implements the fit of the calcium-based plasticity model to experimental data obtained with regular spike-pair stimulation protocol. 
- It then uses the obtained parameter set to predict plasticity results for irregular spike-pairs and burts. 
- Includes scripts for reproducing figures from the associated manuscript

## Script  Reference

The below table describe the indidvual scripts and their purpose. 

| Script                                         | Purpose                                                                                                                                                                 | Input                                                                                                                                                        | Output                                                                                                                                                                                  | Notes                                                      |
|------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------|
| [`calcium_model_fit.py`](calcium_model_fit.py) | Fits the calcium-based plasticity model to the experimental of 1, 3, 5 and 10 Hz regular spike-pair stimulation.                                                        | Experimental data for the four stimulation protocols (in [experimental_data](experimental_data/)                                                             | Saves solutions in the pickled [`solutions.p`](solutions.p) file                                                                                                                        | Parameter fitting script                                   |
| [`checkOutcomes.py`](checkOutcomes.py)         | Generates figures from regular and irregular spike-pair stimulatoin containing regular exp. data and mdoel fit/predictions                                              | The paramter sets for which the figure is generated is specified in the script. Note that the publicaiton is based on the parameter set named `VenancesBin0` | PDF/PNG matplotlib figures for [regular](FigsSimResults/regularDataFitVenance_VenancesBin0.pdf) and [irregular stim.](FigsSimResults/regular-irregular-DataFitVenance_VenancesBin0.pdf) | Plotting model output                                      |
| [`fig_ModelRegularIrregularFrequencies.py`](fig_ModelRegularIrregularFrequencies.py)  | Generates and plots the model prediction for irregular spike-pair data based the parameter set obtained from fitting regular spike-pair data | Uses the fitted parameter set `VenancesBin0`                                                                                                                              | PDF/PNG Figure with the model plasticity outcome for regular and irregular spike-pairs <br/>   <br/>![Analysis pipline](publicationFigures/fig_regularIrregularFrequencies_VenancesBin0_v3.png)                                                                                              | Computation and visualization of irregular plasticity data |

## Requires

All scripts are running with **Python 3**.
Standard python packages such as **numpy**, **scipy**, **pylab**, **time**, **os**,  **sys** and **matplotlib** are required.

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0).
Note that the software is provided "as is", without warranty of any kind, express or implied.
If you use the code or data, please cite us!.


