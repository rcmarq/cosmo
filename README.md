# cosmo

`cosmo` is a suite of add-on functions to the original <a href="">Double Spike Toolbox</a>. These functions allow for the optimization of a double spike suited for analyzing samples with mass-independent isotope effects (i.e., radiogenic excess, nucleosynthetic anomalies). These scripts can simply be copied and pasted to an existing download of the Double Spike Toolbox that has the 'fixed voltage' error model (v. 1.02 or later). 

### New functionalities
The six additional (6) functions are described below:

1. `calcratiocovIN.m`<br>
Calculates the covariance matrix of the ratios that propagates the additional uncertainty from internal normalization. Ideally added to the 'private' folder of the original Double Spike Toolbox. 

2. `shake.m`<br>
Sets up the parameters necessary for performing calculations with the `errorwsplit.m` and all dependent functions.

3. `errorwsplit.m`<br>
Similar to the original `errorestimate.m` function, but considers the splitting of the sample between spiked and unspiked measurements.

4. `cosmo.m`<br>
Similar to the original `cocktail.m` function, but considers the splitting of the sample between spiked and unspiked measurements.

5. `errorcurve3.m`<br>
Error curve that shows uncertainties in &#945; as a function of the sample splitting between spiked and unspiked measurements. The user also has the option to display the uncertainty in the standard/unspiked measurement in &#949; units.

6. `errorsurface.m`<br>
2d projection of error surface on &#945 (or corresponding ppm/amu uncertainty) as a function of the sample splitting and the spike proportion in the spike-sample mixture. 
