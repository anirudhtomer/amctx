temp = pDynStopTime(simJointModel_replaced, newdata = patientDsList[[1]], Dt = 1, K = 50, seed = 2015, idVar="amctx")
temp$summary$times
exp(temp$summary$Info)

temp2 = dynInfo(simJointModel_replaced, newdata = patientDsList[[1]], Dt = 0.26, K = 50, seed = 2015, idVar="amctx")
temp3 = dynInfo_mod(simJointModel_replaced, newdata = patientDsList[[1]], Dt = 0.26, K = 50, seed = 2015, idVar="amctx")
