test_that("predict", {
  
  skip_on_cran()
  
  # small test to run fast
  data("Male_Gammarus_Single")
  Male_Gammarus_Single <- Male_Gammarus_Single[Male_Gammarus_Single$replicate == 1, ]
  modelData_MGS <- modelData(Male_Gammarus_Single, time_accumulation = 4)
  fit_MGS <- fitTK(modelData_MGS, iter = 10, chains = 2)
  
  data("Male_Gammarus_seanine_growth")
  modelData_MGSG <- modelData(Male_Gammarus_seanine_growth, time_accumulation = 1.417)
  fit_MGSG <- fitTK(modelData_MGSG, iter = 10, chains=2)

  data("Chiro_Creuzot")
  Chiro_Creuzot <- Chiro_Creuzot[Chiro_Creuzot$replicate == 1,]
  modelData_CC <- modelData(Chiro_Creuzot, time_accumulation = 1.0)
  fit_CC <- fitTK(modelData_CC, iter = 10, chains=2)
  
  # SAME DATA for test
  # Data 4 prediction should respect the exposure routes
  data_4pred_MGS <- data.frame( time = 0:25, expw = 7.08e-05)
  predict_MGS <- predict(fit_MGS, data_4pred_MGS, fixed_init = TRUE)
  plot(fit_MGS)
  plot(predict_MGS)
  
  predict_MGS <- predict(fit_MGS, data_4pred_MGS, fixed_init = FALSE)
  plot(fit_MGS)
  plot(predict_MGS)
  
  ###############
  data_4pred_MGSG <- data.frame(time = sort(c(0:6,1.417)), expw = 15.533)
  predict_MGSG <- predict(fit_MGSG, data_4pred_MGSG)
  plot(fit_MGSG)
  plot(predict_MGSG)
  
  data_4pred_CC <- data.frame( time = seq(0,4,0.5), expw = 22.9, exps = 1315.7, exppw = 16.24)
  predict_CC <- predict(fit_CC, data_4pred_CC)
  plot(fit_CC)
  plot(predict_CC)
  
  expect_true(all(class(plot(predict_MGS)) == c("gg","ggplot")))
  expect_true(all(class(plot(predict_MGSG)) == c("gg", "ggplot")))
  expect_true(all(class(plot(predict_CC)) == c("gg", "ggplot")))
  
})

  