# mbDriver
mbDriver: identifying driver microbes in microbial communities based on time-series microbiome data

## Framework of mbDriver
![](https://github.com/tanxiaoxiu/mbDriver/blob/master/Framework.png)

## a. Data preprocessing
Observed temporal abundance data are denoised and smoothed, using smoothing splines based on the negative binomial distribution, to obtain estimates of the species abundance curves and their derivatives.
## b. Parameter estimation 
The generalized Lotka-Volterra (gLV) equations and regularized least squares are employed for dynamic modeling and estimation of the growth rates and interaction parameters in the gLV model, respectively. 
## c. Driver prediction
The Driver score index is introduced from a causal perspective for identifying driver microbes in a microbial community. This index quantifies the impact of each microbe on changes in the community’s steady state, and is derived from manipulating the causal graph implied by the gLV equations. A higher score for a microbe indicates its capability to induce more significant changes in the community.
