test_that("PriorPost",{
  
  # Very small fit to run super fast
  data("Male_Gammarus_Single")
  modelData_MGS <- modelData(Male_Gammarus_Single, time_accumulation = 4)
  fit_MGS <- fitTK(modelData_MGS, iter = 10, chains = 2)
  data("Male_Gammarus_seanine_growth")
  modelData_MGSG <- modelData(Male_Gammarus_seanine_growth, time_accumulation = 1.417)
  fit_MGSG <- fitTK(modelData_MGSG, iter = 10, chains=2)
  
  test_that("df_PriorPost", {
    
    df_MGS <- df_PriorPost(fit_MGS)
    df_MGSG <- df_PriorPost(fit_MGSG)
    
    # Check class
    expect_true(class(df_MGS) == "data.frame")
    expect_true(class(df_MGSG) == "data.frame")
    
    # Check column names 
    expect_true(all(colnames(df_MGS) == c("parameter", "type", "value")))
    expect_true(all(colnames(df_MGSG) == c("parameter", "type", "value")))
    
    # Check parameter of prior and posterior are equals
    df_post_MGS = df_MGS[df_MGS$type == "posterior",]
    df_prior_MGS = df_MGS[df_MGS$type == "prior",]
    expect_true(all(df_post_MGS$parameter == df_prior_MGS$parameter))
    
    
  })
  
  test_that("plot_PriorPost", {
    
    plt_MGS <- plot_PriorPost(fit_MGS)
    plt_MGSG <- plot_PriorPost(fit_MGSG)
    
    # Check class
    expect_true(all(class(plt_MGS) == c("gg", "ggplot")))
    expect_true(all(class(plt_MGSG) == c("gg", "ggplot")))
    
  })
  
})



