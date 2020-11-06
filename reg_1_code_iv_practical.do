//Neil Davies 06/11/20
//This do-file runs through the code you will need for the practical for the principals of instrumental variable analysis

//Practical
//Open dataset
use iv_practical_2020, clear

//1. Describe data
describe
summarize

//2. Generate instrument
bysort physician_id (visit_order):gen prior_prescription=prescribed_cox_2[_n-1]

//3. Test instrument strength
regress prescribed_cox_2 prior_prescription, robust
regress prescribed_cox_2 prior_prescription age female, robust

//4. Test whether the instruments associate with the age and sex
regress age prior_prescription, robust
regress female prior_prescription, robust

//5. Estimate the association of treatment with the outcome with and without adjustment
regress has_gi_event prescribed_cox_2, robust
regress has_gi_event prescribed_cox_2 age female, robust

//6. Manually calculate the Wald estimator
regress prescribed_cox_2 prior_prescription, robust
scalar denominator =_b[prior_prescription]
regress has_gi_event prior_prescription, robust
scalar numerator = _b[prior_prescription]
display scalar(numerator)/scalar(denominator)

//7. Use the ivregress command to estimate the effect of prescribing COX-2s on the outcome
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription), robust

//8. Perform 2SLS manually using two linear regressions
regress prescribed_cox_2 prior_prescription
predict prescribed_cox_2_hat
regress has_gi_event prescribed_cox_2_hat, robust

//9. Create multiple instruments
bysort physician_id (visit_order):gen prior_prescription2= prescribed_cox_2 [_n-2]
bysort physician_id (visit_order):gen prior_prescription3= prescribed_cox_2 [_n-3]

//10. Run instrumental variable regression with multiple instruments
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription prior_prescription2 prior_prescription3), robust

//11. Run the instrumental variable regression seperately for each instrument
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription), robust
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription2), robust
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription3), robust

//12. Investigate the plausibility of the instrumental variable assumptions using the in built diagnostics 
ivregress 2sls has_gi_event (prescribed_cox_2 = prior_prescription*), robust
estat firststage, all
estat endogenous 
estat overid
