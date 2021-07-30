## Dynamic linear model (DLM)
### File descriptions
1. ``dlm_functions.py`` contains the key functions used for the DLM model. This is modified from Liu et al., (2019).
2. ``dlm_ndvi_analysis.py`` contains the main function for DLM implementation. Four parameters should be passed to this function, including the spin up time (5 years), spin up repeat times (2 times), delta values (0.98), and model experiment id. In this example, we only provide the setup, i.e., detault model (temperature and precipitation), reduced model (only precipitation), extended model (temperature, precipitation, and radiation), and randomized experiment (using the default model with different random seeds).
3. ``get_component_from_DLM.R`` generates the component (contribution) of each independent variable. These components will be further used for overshoot identification. 

**References**

Liu, Y., Kumar, M., Katul, G.G., Porporato, A., 2019. Reduced resilience as an early warning signal of forest mortality. Nat. Clim. Chang. 9, 880–885. https://doi.org/10.1038/s41558-019-0583-9
