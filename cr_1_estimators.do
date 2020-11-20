//Neil Davies 19/11/20
//Simulations for lecture 3 estimation using instrumental variables


clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
gen instrument=rnormal()

//Create exposure
egen exposure=rowtotal(u c_* instrument)

//Create outcome
egen outcome=rowtotal(exposure c_* v )


//The true causal effect is equal to 1
reg outcome exposure,ro

//But the regression coefficient is biased and equal to 1.9

//Using ivreg2
ivreg2 outcome (exposure =instrument),ro first endog(exposure)

//*******************************************************
//*************With mutliple instruments ****************
//*******************************************************

clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
forvalues i=1(1)10{
	gen instrument_`i'=rnormal()
	}

//Create exposure
egen exposure=rowtotal(u c_* instrument_*)

//Create outcome
egen outcome=rowtotal(exposure c_* v )


//The true causal effect is equal to 1
reg outcome exposure,ro

//But the regression coefficient is biased and equal to 1.9

//Using ivreg2
ivreg2 outcome (exposure =instrument_*) c_1,ro first endog(exposure)

//*******************************************************
//*************Causal risk difference ************************
//*******************************************************


clear
set obs 100000

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
forvalues i=1(1)10{
	gen instrument_`i'=rnormal()
	}

//Create exposure
egen exposure=rowtotal(c_* instrument_*)

//Create outcome
//Note no error terms are required
egen outcome=rowtotal( c_*)

//Add in the effect of the exposure
set seed 1234
replace outcome=outcome+0.2*exposure

//Generate the potential outcomes
set seed 1234
gen outcome_bin0=rbinomial(1,1/(1 + exp(-outcome)))
set seed 1234
gen outcome_bin1=rbinomial(1,1/(1 + exp(-outcome-0.2)))

//The true risk difference 
gen rd=outcome_bin1-outcome_bin0
sum rd

//The LPM is biased and over estiamtes the effect of the exposure
reg outcome_bin0 exposure

//We can estimate with ivreg2 and robust standard errors
ivreg2 outcome_bin0 (exposure=instrument_*),ro

//The IV estimate is consistent

//We can also use the gernalised method of moments to estimate the same parameter using an additive models
gmm (outcome_bin0- {ey0} - exposure*{psi}), instruments(instrument_*)
estat overid

//Where psi is the estimate of the risk difference

//************************************************************
//*****************Two residual inclusion estimator (2SRI) ***
//************************************************************

reg exposure instrument_*
predict x_hat2, xb
predict res_hat, res

logistic outcome_bin0 x_hat2 res_hat

//Note the standard errors are likely to be underestimated because of the estimation error in the first stage.
//Can fix the standard errors using bootstrapping

cap prog drop resid_incl
prog def resid_incl
cap drop x_hat2 res_hat

reg exposure instrument_*
predict x_hat2, xb
predict res_hat, res
logistic outcome_bin0 x_hat2 res_hat

end 

bootstrap , reps(1000) : resid_incl

//Standard errors are increased by 32% and give the correct coverage


//**************************************************************
//********************Causal risk ratio*************************
//**************************************************************

sum outcome_bin0 outcome_bin1 

tabstat outcome_bin0 outcome_bin1,save

//True causal risk ratio
di r(StatTotal)[1,2]/r(StatTotal)[1,1]
//Causal risk ratio=1.0375

//If we want to estimate 
logistic outcome_bin0 exposure 

//But the regression coefficient is biased and equal to 1.52, upwards biased

//Using the gmm (generalized method of moments function)
//Using an multiplicative function for the causal risk ratio

gmm (outcome_bin0*exp(-1*exposure*{psi}) - {ey0}), instruments(instrument_*)
lincom [psi]_cons, eform 
//causal risk ratio=1.043
estat overid


gmm (outcome_bin0*exp(-exposure*{psi} - {logey0}) - 1), instruments(instrument_*)
lincom [psi]_cons, eform //
//GMM estimate causal risk ratio=1.0436

//**************************************************************
//********************Causal odds ratio*************************
//**************************************************************

//The true causal odds-ratio is

tab outcome_bin0 
/*
outcome_bin |
          0 |      Freq.     Percent        Cum.
------------+-----------------------------------
          0 |     49,919       49.92       49.92
          1 |     50,081       50.08      100.00
------------+-----------------------------------
      Total |    100,000      100.00
*/

tab outcome_bin1 

/*
outcome_bin |
          1 |      Freq.     Percent        Cum.
------------+-----------------------------------
          0 |     48,036       48.04       48.04
          1 |     51,964       51.96      100.00
------------+-----------------------------------
      Total |    100,000      100.00
*/

di (51964/48036)/(50081/49919)
//True 1.0782727

//Estimate causal odds-ratio using GMM:
//Generate interactions
forvalues i=1(1)10{
	gen exposure_iv_`i'=exposure*instrument_`i'
	}
	
//Generate association model
logit outcome_bin0 exposure instrument_* exposure_iv_*
matrix from = e(b)
predict xblog, xb 	

//Causal model with incorrect SEs
gmm (invlogit(xblog - exposure*{psi}) - {ey0}), instruments(instrument_*)
lincom [psi]_cons, eform
matrix from = (from,e(b))

//Causal model with correct SEs
gmm (outcome_bin0 - invlogit({xb:exposure instrument_* exposure_iv_*} + {b0})) ///
	(invlogit({xb:} + {b0} - exposure*{psi}) - {ey0}), ///
	instruments(1:exposure instrument_* exposure_iv_*) ///
	instruments(2:instrument_*) ///
	winitial(unadjusted, independent) from(from) 

lincom [psi]_cons, eform
//causal odds ratio
estat overid


//************************************************************
//*****************Two sample estimator***********************
//************************************************************


clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
forvalues i=1(1)10{
	gen instrument_`i'=rnormal()
	}

//Create exposure
egen exposure=rowtotal(u c_* instrument_*)

//Create outcome
egen outcome=rowtotal(exposure c_* v )

//Suppose the exposure is very expensive to measure and is only recorded in 10% of the sampple
replace exposure=. if _n>_N*0.1

//Can use a two-stage least squares estimator 
reg exposure instrument_*,robust
predict exposure_hat if exposure==.

//The true causal effect is equal to 1
reg outcome exposure,ro

//Estiamte association with the prediction
reg outcome exposure_hat,ro

//Compare to running the IV analysis in the 10000 with exposure measures:
ivreg2 outcome (exposure =instrument_*),ro first endog(exposure)

//************************************************************
//*****************Two exposure estimators********************
//************************************************************


clear
set obs 100000

//Create error terms
gen u_x1=rnormal()
gen u_x2=rnormal()
gen v=rnormal()

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
forvalues i=1(1)10{
	gen instrument_x1_`i'=rnormal()
	}
forvalues i=1(1)10{
	gen instrument_x2_`i'=rnormal()
	}

//Create exposure
egen exposure_x1=rowtotal(u_x1 c_* instrument_x1_*)
egen exposure_x2=rowtotal(u_x2 c_* instrument_x2_*)

//Create outcome
egen outcome=rowtotal(exposure* c_* v )

//Use ivreg2 to estimate the effects of multiple exposures
ivreg2 outcome (exposure* =instrument_*),ro first endog(exposure_*)
