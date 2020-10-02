//Neil Davies 03/10/20
//This simulates selection/collider bias

clear
set obs 100000

//Create error terms
gen u=rnormal()
gen v=rnormal()


//Create exposure
gen exposure=u

//Create outcome
gen outcome=exposure+v

//The true causal effect is equal to 1
reg outcome exposure

//However suppose there is a selection process into the study, and we only sample those who have higher values of the outcome and the exposure
gen selection=0.2*(outcome+exposure)+0.8*rnormal()

//Running a regression in the selected sample results in bias
reg outcome exposure if selection>0

//Coefficient, which is equal to 1 above (there are now no confounders), is equal to 0.93 in the selected sample.

//Similarly if the collider is adjusted for, the regression coefficient is biased
reg outcome exposure selection 

//The coefficient is now equal to 0.89.
