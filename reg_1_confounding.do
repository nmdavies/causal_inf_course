//Neil Davies 02/10/20
//A simulation of confounding and selection bias

clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}

//Create exposure
egen exposure=rowtotal(u c_*)

//Create outcome
egen outcome=rowtotal(exposure c_* v)

//The true causal effect is equal to 1
reg outcome exposure

//But the regression coefficient is bias and equal to 1.9

//We can adjust for the confounders
reg outcome exposure c_*

//To recover the causal estimate=1.

//However in most applied examples, it is unlikely to be possible measure all the confounders with sufficient precision. 
//For example, suppose we only measure the covariates with error:

//Create mismeasured covariates
forvalues i=1(1)10{
	gen c_measured_`i'=c_`i'+rnormal()
	}
	
//If we adjust for the measured covariates
reg outcome exposure c_measured_*

//The estimated effect is only 1.8, and highly biased. 

//Note this is only with a very modest amount of measurement error - it assumes that we've measured variables that account for 50%
//of the variation in each confounder. This is pretty good compared to many real life studies, where complex covariates like 
//"smoking history" are accounted for via a single binary covariate. 

//Thus - it is generally implausible that measure covariates can sufficiently account for all differences between exposure groups. 
