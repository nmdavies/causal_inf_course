//Neil Davies 11/12/20 
//Illustrating the importance of point identifying assumptions

foreach n in 10000 100000 1000000{
	clear
	set seed 2232011
	set obs `n'

	generate z = rbinomial(1,.5)
	generate u = rbinomial(1,.5)
	generate px = .05 + .1*z + .1*u
	generate x = rbinomial(1,px)
	generate p = .1  + .1*x + .1*u
	generate y = rbinomial(1,p)

	//Check frequency of the outcome, exposure and instrument
	tab y
	tab x
	tab z

	//The bounds can be calculated using the following command (Palmer et al. 2011)
	bpbounds y (x = z)

	//Whereas if we are willing to assume a 4th point identifying assumption our results are much more precise:
	ivreg2 y (x=z),ro	
}
*************************************************************************************
*************************CONSTANT TREATMENT EFFECT***********************************
*************************************************************************************
***************************CONTINIOUS OUTCOMES***************************************

//Simulation illustrating a) a constant effect of the exposure on a continious outcome


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

//Generate true potential outcomes
gen y_0=outcome
gen y_1=outcome+1

gen TE=y_1-y_0

sum TE

//Treatement effect the same in everyone, and obviously equal to 1.

*************************************************************************************
***************************BINARY OUTCOMES***************************************

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

gen TE=outcome_bin1-outcome_bin0

sum TE

//Note that while the exposure has a constant effect on the log odds scale. It has a variable effect on the binary outcome.
//For some individuals increasing the exposure causes the outcome, for others it cures the outcome, for most, it has no effect.
//Therefore, while the effects are constant on the log odds scale, they are variable on the binary outcome.
//The only way to have a constant effect on a binary outcome is if the exposure entirely determines the outcome. 



*************************************************************************************
*************************HETEROGENOUS TREATMENT EFFECT***********************************
*************************************************************************************
***************************CONTINIOUS OUTCOMES***************************************

clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Generate treatment effects
gen treatment=rnormal(1,1)

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
gen instrument=rnormal()

//Create exposure
egen exposure=rowtotal(u c_* instrument)

//Create outcome
egen outcome=rowtotal( c_* v )
replace outcome=outcome+treatment*exposure

//The true causal effect is equal to 1
reg outcome exposure,ro

//But the regression coefficient is biased and equal to 1.9

//Using ivreg2
ivreg2 outcome (exposure =instrument),ro first endog(exposure)

//Generate true potential outcomes
gen y_0=outcome
gen y_1=outcome+1*treatment

gen TE=y_1-y_0

sum TE


//Suppose we have a subgroup that are not given the treatment because it is very detrimental, but that this subgroup is not affected by the instrument
clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()

//Generate treatment effects
gen treatment=rnormal(1,1)
//Replace treatment effect equal to -1 in 10% of the sample
replace treatment =-1 if _n<10000

//Create confounders
forvalues i=1(1)10{
	gen c_`i'=rnormal()
	}
	
//Create the instrumental variables
gen instrument=rnormal()

//Create exposure
egen exposure=rowtotal(u c_* instrument)
replace exposure=exposure-instrument if _n<10000

//Create outcome
egen outcome=rowtotal( c_* v )
replace outcome=outcome+treatment*exposure

//The true causal effect is equal to 1
reg outcome exposure,ro

//But the regression coefficient is biased and equal to 1.9

//Using ivreg2
ivreg2 outcome (exposure =instrument),ro first endog(exposure)

//Generate true potential outcomes
gen y_0=outcome
gen y_1=outcome+1*treatment

gen TE=y_1-y_0

sum TE
//The treatment effect in those without the adverse effects: 
sum TE if _n>=10000

//This is exactly the sub-group identified by the IV analysis. 

