## Overshoot identification
### File description
**``identify_overshoot.R``**  contains functions for overshoot drought identification. Some of the major functions are described below.

>``drought_identify`` functions used to identify drought event for each pixel. 
>
>**Inputs** vegetation index (VI) time series, SPEI, precipitation sensitivity, temeprature sensitivity, precipitation contribution, temeprature contribution from DLM, and mean vegetation index.
>
>**Outputs** time series of drought event id number, e.g., 0 0 0 1 1 1 0 0 2 2 2 2 .... represents first drought happended between April to June, second drought happened in September to December.   
