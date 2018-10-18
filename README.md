# A Personalized Screening Strategy to Monitor the Development of Chronic Allograft Failure in Renal Transplant Recipients

## Code for the Model
We use a joint model for time to event and longitudinal data to model the evolution of serum creatinine (SCr) and protein-creatinine ratio (PCR) levels over time, and to simulatenouly model their association with the hazard of graft failure. The R package we use for this purpose is called JMbayes (https://cran.r-project.org/web/packages/JMbayes/JMbayes.pdf). The API we use however are currently not hosted on CRAN, and can be found here:
https://github.com/drizopoulos/JMbayes/blob/master/man/mvJointModelBayes.Rd

Now since public access to kidney transplant data is forbidden, you will most probably be not able to run our code. Having said that, the code for the model that we fit can be found in this repository. For the relative risk part, we first create a model using the well known coxph API (https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html). The code for fitting the models can be found here:
https://github.com/anirudhtomer/amctx/blob/master/src/coxAnalysis.R

As you can see we fit multiple models, however the one we use finally is called "coxModel_clinical".

Using the definitions from the coxph model object, we then create a joint model.
https://github.com/anirudhtomer/amctx/blob/master/src/jointAnalysis.R

As you can see we fit multiple models, however the one we use finally is called "mvJoint_pcr_creatinine_tdboth_complex". For the joint model with only SCr as the longitudinal outcome the joint model we use is "mvJoint_creatinine_tdboth_complex".  

## Code for the Personalized Schedules, and Simulation Study
Using the fitted joint model with only SCr as the longitudinal outcome, we simulate evolution of SCr for 625 patients. We also generate the time of graft failure for these simulated patients. The simulation code can be found here (start reading from the bottom to the top):
https://github.com/anirudhtomer/amctx/blob/master/src/Simulation%20Study/simCommon.R

Now for 50 of these patients we also generate personalized schedules, and the fixed schedule of SCr measurements. The code for the fixed schedule can be found here:
https://github.com/anirudhtomer/amctx/blob/master/src/Simulation%20Study/FixedSchedule.R
The code for the personalized schedule can be found here:
https://github.com/anirudhtomer/amctx/blob/master/src/Simulation%20Study/pers_schedule_entr_KL.R

Then we run simulations for 50 test patients to check which schedule conducts how many SCr measurements how much do they overshoot/undershoot the real target time at which the dynamic risk of graft failure is more than 5%. This can be found here:
https://github.com/anirudhtomer/amctx/blob/master/src/Simulation%20Study/pers_schedule_entr_KL.R

Lastly, we generate graphs to see our results, which can be found here:
https://github.com/anirudhtomer/amctx/blob/master/src/Simulation%20Study/produceResults.R
