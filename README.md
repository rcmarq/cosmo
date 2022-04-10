# cosmo

`cosmo` is a suite of add-on functions to the original <a href="">Double Spike Toolbox</a>. These functions allow for the optimization of a double spike suited for analyzing samples with mass-independent isotope effects (i.e., radiogenic excess, nucleosynthetic anomalies). These scripts can simply be copied and pasted to an existing download of the Double Spike Toolbox that has the 'fixed voltage' error model (v. 1.02 or later). 

## New Functionalities
The five additional (5) functions are described below:
<ol>
<li>`calcratiocovIN.m` - Calculates the covariance matrix of the ratios that propagates the additional uncertainty from internal normalization</li>
<li>`shake.m` - Sets up the parameters necessary for performing calculations with the `errorwsplit.m`,`cosmo.m`, and `errocurve3.m` functions</li>
<li>`errorwsplit.m` - Similar to the original `errorestimate.m` function, but considers the splitting of the sample between spiked and unspoken measurements.</li>
<li>`cosmo.m` - Similar to the original `cocktail.m` function, but considers the splitting of the sample between spiked and unspiked measurements.</li>
<li>`errorcurve3.m` - Error curve that shows uncertainties in $\alpha$ as a function of the sample splitting between spiked and unspiked measurements. The user also has the option to display the uncertainty in the standard/unspiked measurement in $\varepsilon$ units.</li>
</ol>