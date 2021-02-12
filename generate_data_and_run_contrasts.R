  library(plyr)
  library(gamlss.dist)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(ARTool)
  library(emmeans)
  library(data.table)
  library(ggm)
  library(lme4)
  library(nimble)
  library(broom)
  library(hash)
  #library(ARTool)
  
  source("artcon_draft.R")
  
  ####### ANALYZE AND LOG DATA #######
  
  # make full factorial ART model using all factors from data
  # e.g. m = art(Y ~ X1*X2*X3, data)
  # RETURNS: list of art models
    # e.g. if between-Ss list("Fixed" = art(Y ~ X1*X2*X3, data))
    # e.g. if within-Ss list("Fixed" = art(Y ~ X1*X2*X3 + (1|S), data), "Random" = art(Y ~ X1*X2*X3 + Error(S)))
  makeFullFactorialARTModel = function(data, factorNames, between){
    # string models first
    
    # between
    # e.g. factors in data are X1, X2, X3 -> Y ~ X1*X2*X3
    strARTFormula = paste("Y ~ ", paste(factorNames, collapse = "*"), sep="")
    
    if(between){
      # turn into formula
      ARTFormulaFixed = as.formula(strARTFormula)
      mFixed = art(ARTFormulaFixed, data)
      return(list("Fixed" = mFixed))
    }
    
    else{
      # random
      # e.g. factors in data are X1, X2, X3 -> Y ~ X1*X2*X3 + (1|S)
      strARTFormulaRandom = paste(strARTFormula, "(1|S)", sep=" + ")
      ARTFormulaRandom = as.formula(strARTFormulaRandom)
      mRandom = art(ARTFormulaRandom, data)
      # fixed
      # e.g. factors in data are X1, X2, X3 -> Y ~ X1*X2*X3 + Error(S)
      strARTFormulaFixed = paste(strARTFormula, "Error(S)", sep=" + ")
      ARTFormulaFixed = as.formula(strARTFormulaFixed)
      mFixed = art(ARTFormulaFixed, data)
      return(list("Fixed" = mFixed, "Random" = mRandom))
    }
  }
  
  # contrast(emmeans(artlm(m.art, "X1:X2"), ~ X1:X2), method="pairwise")
  # include metaLogRow because change "Omnibus_ART_Warning" to 1 if there the contrast produces a warning
  testOmnibusART = function(data, fullFactorialARTModel, factorSubset, metaLogRow){
    # e.g. factorSubset = ("X1, "X2") -> X1:X2"
    contrastFormulaStr = paste(factorSubset, collapse=":")
    # e.g. factorSubset = ("X1, "X2") -> ~ X1:X2
    contrastFormula = as.formula(paste("~ ", contrastFormulaStr, sep=""))
    
    # do contrasts TODO comment out
    tempResult = contrast(emmeans(artlm(fullFactorialARTModel, contrastFormulaStr), 
                     contrastFormula), method="pairwise", adjust="none")
    
    tc = tryCatch(
      list("result" = contrast(emmeans(artlm(fullFactorialARTModel, contrastFormulaStr), contrastFormula), method="pairwise", adjust="none")),
      warning = function(w) return(list("result" = contrast(emmeans(artlm(fullFactorialARTModel, contrastFormulaStr), contrastFormula), method="pairwise", adjust="none"),
                                         "warning" = w)))
    
    omnibusARTResult = tc[["result"]]
    warning = tc[["warning"]]
    
    if(!is.null(warning)){
      metaLogRow[["Omnibus_ART_Warning"]] = 1
    }
    
    return(list("result" = omnibusARTResult, "metaLogRow" = metaLogRow))
    
  }
  
  # contrast using aligned data from fullFactorialARTModel instead of aligned and ranked.
  # using art.con function
  testart.conAligned = function(data, fullFactorialARTModel, factorSubset){
    contrastFormulaStr = paste(factorSubset, collapse=":")
    art.con(fullFactorialARTModel, contrastFormulaStr, response = "aligned", adjust="none")
  }
  
  # contrasts using aligned data from fullFActorialARTModel instead of aligned and ranked
  # using artlm.con function
  testArtlmconAligned = function(data,fullFactorialARTModel, factorSubset){
    # e.g. factorSubset = ("X1, "X2") -> X1:X2"
    contrastStr = paste(factorSubset, collapse=":")
    
    # e.g. factorSubset = c("X1, "X2") -> "~ X1X2"
    contrastFormulaStr = paste("~ ", paste(factorSubset, collapse=""), sep="")
    # e.g. factorSubset = c("X1, "X2") -> ~ X1X2
    contrastFormula = as.formula(contrastFormulaStr)
    
    # contrast
    contrast(emmeans(artlm.con(fullFactorialARTModel, contrastStr, response = "aligned"), contrastFormula), method="pairwise", adjust="none")
  }
  
  # takes meta log row as argument because will log if there is a warning
  # Returns: list("art.conContrast" = art.con(fullFactorialARTModel, contrastFormulaStr, adjust="none"), "metaLogRow" = metaLogRow)
  testart.con = function(data, fullFactorialARTModel, factorSubset, metaLogRow){
    contrastFormulaStr = paste(factorSubset, collapse=":")
    
    # art con try catch
    art.conTryCatch = tryCatch(
      list("art.conResult" = art.con(fullFactorialARTModel, contrastFormulaStr, adjust="none")),
      warning = function(w) return(list("art.conResult" = art.con(fullFactorialARTModel, contrastFormulaStr, adjust="none"),
                                        "warning" = w)))
    art.conResult = art.conTryCatch[["art.conResult"]]
    warning = art.conTryCatch[["warning"]]
  
    if(!is.null(warning)){
      metaLogRow[["ART_Con_Warning"]] = 1
    }
    
    return(list("art.conContrast" = art.conResult, "metaLogRow" = metaLogRow))
  }
  
  # syntax: contrast(emmeans(artlm.con(m, "X1:X2"), ~ X1X2), method="pairwise")
  testArtlm.con = function(data, fullFactorialARTModel, factorSubset){
    # e.g. factorSubset = ("X1, "X2") -> X1:X2"
    contrastStr = paste(factorSubset, collapse=":")
    
    # e.g. factorSubset = c("X1, "X2") -> "~ X1X2"
    contrastFormulaStr = paste("~ ", paste(factorSubset, collapse=""), sep="")
    # e.g. factorSubset = c("X1, "X2") -> ~ X1X2
    contrastFormula = as.formula(contrastFormulaStr)
    
    # contrast
    contrast(emmeans(artlm.con(fullFactorialARTModel, contrastStr), contrastFormula), method="pairwise", adjust="none")
  }
  
  # modelType is optional param with default fixed. For within subjects, option to use "Random"
  # for between subjects, will always be fixed even if random passed in.
  testParam = function(data, factorNames, factorSubset, between, modelType = "Fixed"){
    # e.g. factorNames = c("X1", "X2", "X3") -> "Y ~ X1*X2*X3"
    modelFormulaStr = paste("Y ~ ", paste(factorNames, collapse="*"), sep="")
    
    # between-Ss
    if(between){
      model = lm(as.formula(modelFormulaStr), data = data)
    }
    # within-Ss
    else {
      if(modelType == "Fixed"){
        # fixed effects
        # e.g. factorNames = c("X1", "X2", "X3") -> "Y ~ X1*X2*X3 + Error(S)"
        modelFormulaStrFixed = paste(modelFormulaStr, "Error(S)", sep="+")
        model = aov(as.formula(modelFormulaStrFixed), data = data)
      }
      else if(modelType == "Random"){
        # random effects
        # e.g. factorNames = c("X1", "X2", "X3") -> "Y ~ X1*X2*X3 + (1|S)"
        modelFormulaStrRandom = paste(modelFormulaStr, "(1|S)", sep="+")
        model = lmer(as.formula(modelFormulaStrRandom), data = data)
      }
      # error if model type isn't fixed or random
      else{
        stop()
      }
    } # end within-Ss
    
    # e.g. factorSubsets = c("X1", "X2") -> ~ X1:X2
    contrastFormulaStr = paste("~ ", paste(factorSubset, collapse=":"), sep="")
    contrastFormula = as.formula(contrastFormulaStr)
    
    contrast(emmeans(model, contrastFormula), method="pairwise", adjust="none")
  }
  
  # for non-param test, concatenate col the same way would for art.con
  # normally would do wilcox.test(df[df$X1 == "a" & df$X2 == "a",]$Y, df[df$X1 == "a" & df$X2 == "b",]$Y)
  # concatenated X1&X2 and running pairwise wilcox on new col is the same thing
  testNonParam = function(data, factorSubset, between){
    
    # e.g. factorSubset = c("X1", "X2") -> "X1X2"
    newColName = paste(factorSubset, collapse="")
    # e.g. factorSubset = c("X1", "X2") create new col X1X2 and drop cols X1 and X2
    concatDf = unite(data, !!newColName, factorSubset, sep = ",", remove = TRUE)
    # make new col a factor
    concatDf[[newColName]] = factor(concatDf[[newColName]])
    
    # run non param tests
    attach(concatDf)
    # mann whitney U test
    if(between){
      contrast = pairwise.wilcox.test(concatDf$Y, concatDf[[newColName]], distribution="exact", p.adjust="none")
    }
    # wilcoxon signed rank test
    else{
      contrast = pairwise.wilcox.test(concatDf$Y, concatDf[[newColName]], paired=TRUE, p.adjust="none")
    }
    detach()
    
    contrast
  }
  
  # makes data frame with one row per trial
  # only logs param, nonparam, and art.con - getFactorSubsetLogDfLarge logs all
  # also only logs random model results for within subjects. Large version does fixed and random
  getFactorSubsetLogDfSmall = function(contrastResult, dataFrameNum, numResponsesPerCondition, between, condStrNumMap){
    
    # col names for factorSubsetLogDf
    if(between){
      colNames = c(
        "Data_Set", "Num_Responses_per_Condition", "BW", 
        "Trial", "C1", "Cstr1", "C2", "Cstr2",
        "DF_ART_Omni", "DF_Param", "DF_ART_Con",
        "T_ART_Omni","T_Param", "T_ART_Con",
        "P_ART_Omni", "P_Param", "P_ART_Con",
        "P_Nonparam"
      )
    } # end between
    
    # within
    else {
      colNames = c(
        "Data_Set", "Num_Responses_per_Condition", "BW", 
        "Trial", "C1", "Cstr1", "C2", "Cstr2",
        "DF_ART_Omni_Random", "DF_Param_Random", "DF_ART_Con_Random",
        "T_ART_Omni_Random", "T_Param_Random", "T_ART_Con_Random",
        "P_ART_Omni_Random", "P_Param_Random", "P_ART_Con_Random",
        "P_Nonparam"
      )
    } # end within
  
    # make data from to log results for factor subset with correct number of cells, all values NA
    factorSubsetLogDf = data.frame(matrix(NA, nrow = nrow(contrastResult@grid), ncol = length(colNames)))
    # add column names to factorSubsetLogDf
    colnames(factorSubsetLogDf) = colNames
  
    # add Data_Set, Num_Responses_per_Condition, and BW to all rows since all the same
    factorSubsetLogDf[["Data_Set"]] = rep(dataFrameNum, nrow(factorSubsetLogDf))
    factorSubsetLogDf[["Num_Responses_per_Condition"]] = rep(numResponsesPerCondition, nrow(factorSubsetLogDf))
    factorSubsetLogDf[["BW"]] = rep(if(between) "between" else "within", nrow(factorSubsetLogDf))
    
    # add C1, Cstr1, C2, Cstr2 to each row
    contrastList = contrastResult@grid[[1]]
    for(i in 1:length(contrastList)){
      currContrast = contrastList[i]
      currContrastStr = toString(currContrast)
      # vector with 2 elems. each is a string rep of one level in current contrast
      currContrastStrVec = strsplit(currContrastStr, ' - ')[[1]]
      # named vector with keys as names and values as values, in increasing sorted order by value
      # e.g. currContrastStrVec = c("a,c", "b,d") -> c("a,c" = 2, "b,d" = 3)
      currContrastVec = sort(values(condStrNumMap, keys = currContrastStrVec))
      
      # add stuff for this trial to factorSubsetLogDf
      # Note: str_replace because names in log separated by spaces not comma
      factorSubsetLogDf[i, c("Trial", "C1", "Cstr1", "C2", "Cstr2")] = c(trialNum, 
                      currContrastVec[[1]], 
                      str_replace_all(names(currContrastVec)[[1]], ",", " "), 
                      currContrastVec[[2]], 
                      str_replace_all(names(currContrastVec)[[2]], ",", " ")
                      )
      
      # increment trialNum
      trialNum <<- trialNum + 1
    } # end iterate through contrasts
    factorSubsetLogDf
  }
  
  # makes data frame with one row per trial
  # logs all - getFactorSubsetLogDfSmall only logs nonparam, param, and art.con
  getFactorSubsetLogDfLarge = function(contrastResult, dataFrameNum, numResponsesPerCondition, between, condStrNumMap){
    
    # col names for factorSubsetLogDf
    if(between){
      colNames = c(
        "Data_Set", "Num_Responses_per_Condition", "BW", 
        "Trial", "C1", "Cstr1", "C2", "Cstr2",
        "DF_ART_Omni", "DF_Param", "DF_ART_Con_Aligned", "DF_ARTlm_Con_Aligned", "DF_ART_Con", "DF_ARTlm_Con",
        "T_ART_Omni", "T_Param", "T_ART_Con_Aligned", "T_ARTlm_Con_Aligned", "T_ART_Con", "T_ARTlm_Con",
        "P_ART_Omni", "P_Param", "P_ART_Con_Aligned", "P_ARTlm_Con_Aligned", "P_ART_Con", "P_ARTlm_Con",
        "P_Nonparam"
      )
    } # end between
    
    # within
    else {
      colNames = c(
        "Data_Set", "Num_Responses_per_Condition", "BW", 
        "Trial", "C1", "Cstr1", "C2", "Cstr2",
        "DF_ART_Omni_Fixed", "DF_Param_Fixed", "DF_ART_Con_Aligned_Fixed", "DF_ARTlm_Con_Aligned_Fixed", "DF_ART_Con_Fixed", "DF_ARTlm_Con_Fixed",
        "T_ART_Omni_Fixed", "T_Param_Fixed", "T_ART_Con_Aligned_Fixed", "T_ARTlm_Con_Aligned_Fixed", "T_ART_Con_Fixed", "T_ARTlm_Con_Fixed",
        "P_ART_Omni_Fixed", "P_Param_Fixed", "P_ART_Con_Aligned_Fixed", "P_ARTlm_Con_Aligned_Fixed","P_ART_Con_Fixed", "P_ARTlm_Con_Fixed",
        "DF_ART_Omni_Random", "DF_Param_Random", "DF_ART_Con_Aligned_Random", "DF_ARTlm_Con_Aligned_Random", "DF_ART_Con_Random", "DF_ARTlm_Con_Random",
        "T_ART_Omni_Random", "T_Param_Random", "T_ART_Con_Aligned_Random", "T_ARTlm_Con_Aligned_Random", "T_ART_Con_Random", "T_ARTlm_Con_Random",
        "P_ART_Omni_Random", "P_Param_Random", "P_ART_Con_Aligned_Random", "P_ARTlm_Con_Aligned_Random", "P_ART_Con_Random", "P_ARTlm_Con_Random",
        "P_Nonparam"
      )
    } # end within
    
    # make data from to log results for factor subset with correct number of cells, all values NA
    factorSubsetLogDf = data.frame(matrix(NA, nrow = nrow(contrastResult@grid), ncol = length(colNames)))
    # add column names to factorSubsetLogDf
    colnames(factorSubsetLogDf) = colNames
    
    # add Data_Set, Num_Responses_per_Condition, and BW to all rows since all the same
    factorSubsetLogDf[["Data_Set"]] = rep(dataFrameNum, nrow(factorSubsetLogDf))
    factorSubsetLogDf[["Num_Responses_per_Condition"]] = rep(numResponsesPerCondition, nrow(factorSubsetLogDf))
    factorSubsetLogDf[["BW"]] = rep(if(between) "between" else "within", nrow(factorSubsetLogDf))
    
    # add C1, Cstr1, C2, Cstr2 to each row
    contrastList = contrastResult@grid[[1]]
    for(i in 1:length(contrastList)){
      currContrast = contrastList[i]
      currContrastStr = toString(currContrast)
      # vector with 2 elems. each is a string rep of one level in current contrast
      currContrastStrVec = strsplit(currContrastStr, ' - ')[[1]]
      # named vector with keys as names and values as values, in increasing sorted order by value
      # e.g. currContrastStrVec = c("a,c", "b,d") -> c("a,c" = 2, "b,d" = 3)
      currContrastVec = sort(values(condStrNumMap, keys = currContrastStrVec))
      
      # add stuff for this trial to factorSubsetLogDf
      # Note: str_replace because names in log separated by spaces not comma
      factorSubsetLogDf[i, c("Trial", "C1", "Cstr1", "C2", "Cstr2")] = c(trialNum, 
                                                                         currContrastVec[[1]], 
                                                                         str_replace_all(names(currContrastVec)[[1]], ",", " "), 
                                                                         currContrastVec[[2]], 
                                                                         str_replace_all(names(currContrastVec)[[2]], ",", " ")
      )
      
      # increment trialNum
      trialNum <<- trialNum + 1
    } # end iterate through contrasts
    factorSubsetLogDf
  }
  
  # contrastAndTestName list("contrast" = contrast, "testName" = testName)
  #   contrast is object of type emmGrid with results of one test
  #   testName is the name of the test (needed for column names in log)
  #   logs Df, T, and P
  #   e.g."Param", "Param_Fixed",...
  #   Note: this is for all contrasts except nonParam
  #   
  # adds result of contrast to factorSubsetLogDf
  logResults = function(contrastAndTestName, factorSubsetLogDf, condStrNumMap){
    contrast = contrastAndTestName[["contrast"]]
    
    # iterate through each row of contrast data frame
    contrastDf = as.data.frame(contrast)
    # list of levels in all contrasts ie list("a,c - a,d", "b,c - a,d")
    contrastList = as.character(contrastDf[["contrast"]])
    for(i in 1:length(contrastList)){
      
      # row of contrastDf for this contrast
      contrastRow = contrastDf[i,]
      # get DF (ie degrees of freedom)
      df = contrastRow$df
      # get T statistic
      t.ratio = contrastRow$t.ratio
      # get p value
      p.value = contrastRow$p.value
      
      # find each condition in contrast so can get corresponding condition numbers
      contrastStr = contrastList[i]
      # vector with 2 elems. each is a string representing one level in the contrast
      contrastStrSplit = strsplit(contrastStr, ' - ')[[1]]
  
      # named vector with keys as names and values as values, in ascending sorted order by value
      # e.g. contrastStrSplit = c("a,c", "b,d") -> c("a,c" = 2, "b,d" = 3)
      contrastVec = sort(values(condStrNumMap, keys = contrastStrSplit))
      
      # Need to find row in factorSubsetLogDf with C1 = contrastVec[[1]], and C2 = contrastVec[[2]]
      correctC1 = contrastVec[[1]]
      correctC2 = contrastVec[[2]]
      # logged string names are separated by space not comma
      correctCstr1 = str_replace_all(names(contrastVec)[[1]], ",", " ")
      correctCstr2 = str_replace_all(names(contrastVec)[[2]], ",", " ")
      
      # get row where C1 and C2 match
      matchRow = factorSubsetLogDf %>% filter (C1 == correctC1 & C2 == correctC2)
      # assert that 1. there is only one row in matchRow 2. Cstr1 and Cstr2 also match
      if(nrow(matchRow) > 1 || matchRow[["Cstr1"]] != correctCstr1 || matchRow[["Cstr2"]] != correctCstr2){
        stop("Either more than one row matchs C1 and C2, or Cstr1 and Cstr2 and C1 and C2 don't all match the same row.")
      }
      
      testName = contrastAndTestName[["testName"]]
      # add values to data frame
      for(colList in list(
                    list("colName" = "DF", "valueName" = df), 
                    list("colName" = "T", "valueName" = t.ratio), 
                    list("colName" = "P", "valueName" = p.value)
                    )){
        # make name e.g. testName = "Param_Fixed", logColName = "DF" -> "DF_Param_Fixed"
        fullLogColName = paste(colList[["colName"]], testName, sep="_")
        valueToLog = colList[["valueName"]]
  
        # log
        # add value to row with correct C1 and C2
        factorSubsetLogDf = factorSubsetLogDf %>% dplyr::mutate(!! fullLogColName := ifelse(C1 == correctC1 & C2 == correctC2, valueToLog, factorSubsetLogDf[[!!fullLogColName]]))
        
      } # end add df, t, and p to log for one contrast
    } # end iterate through all contrasts in this test
    factorSubsetLogDf # return updated data frame
  } # end function
  
  # log results from non parametric contrasts
  # contrast results are formatted differently so have to do them separately 
  logNonParamResults = function(nonParamContrast, factorSubsetLogDf, condStrNumMap){
    # get contrast p values as data frame
    nonParamContrastDf = as.data.frame(nonParamContrast$p.value)
    # add column with rownames so can reshape
    nonParamContrastDf$cond1 = rownames(nonParamContrastDf)
    # reshape so each row is cond1, cond2, p.value
    contrastDf = gather(nonParamContrastDf, cond2, p.value, -cond1, na.rm=TRUE)
    
    ### FROM LOG RESULTS
    for(i in 1:nrow(contrastDf)){
      
      # row of contrastDf for this contrast
      # row is cond1, cond2, p.value
      contrastRow = contrastDf[i,]
      # get p value
      p.value = contrastRow$p.value
      
      contrastStrSplit = c(contrastRow$cond1, contrastRow$cond2)
      # named vector with keys as names and values as values, in ascending sorted order by value
      # e.g. contrastStrVec = c("a,c", "b,d") -> c("a,c" = 2, "b,d" = 3)
      contrastVec = sort(values(condStrNumMap, keys = contrastStrSplit))
      
      # Need to find row in factorSubsetLogDf with C1 = contrastVec[[1]], and C2 = contrastVec[[2]]
      correctC1 = contrastVec[[1]]
      correctC2 = contrastVec[[2]]
      
      fullLogColName = "P_Nonparam"
      valueToLog = p.value
      # log
      factorSubsetLogDf = factorSubsetLogDf %>% dplyr::mutate(!! fullLogColName := ifelse(C1 == correctC1 & C2 == correctC2, valueToLog, factorSubsetLogDf[[!!fullLogColName]]))
    } # end iterate through all contrasts in this test
    ### END FROM LOG RESULTS
    factorSubsetLogDf # return updated data frame
  }
  
  # Only runs tests for param, non-param and art.con.
  # testAndLogOneFactorSubsetLarge runs all tests
  # tests one factorSubset (ie combo of factors).
  # logs results
  # factorNames is all factors in data frame
  # fullFactorialARTModelList
  # if between-Ss list("Fixed" = art(Y ~ X1*X2*X3, data))
  # if within-Ss list("Fixed" = art(Y ~ X1*X2*X3 + (1|S), data), "Random" = art(Y ~ X1*X2*X3 + Error(S)))
  testAndLogOneFactorSubsetSmall = function(data, factorNames, fullFactorialARTModelList, factorSubset, dataFrameNum, 
                                       numResponsesPerCondition, between, condStrNumMap, metaLogRow){
    if(between){
      # full factorial ART model used by most contrasts
      fullFactorialARTModel = fullFactorialARTModelList[["Fixed"]]
      
      # param
      paramContrast = testParam(data, factorNames, factorSubset, between)
      
      # non param
      nonParamContrast = testNonParam(data, factorSubset, between)
      
      # omnibus ART
      omnibusart.contrastResult = testOmnibusART(data, fullFactorialARTModel, factorSubset, metaLogRow)
      omnibusart.contrast = omnibusart.contrastResult[["result"]]
      metaLogRow = omnibusart.contrastResult[["metaLogRow"]]
      
      # art.con
      art.conContrastResult = testart.con(data, fullFactorialARTModel, factorSubset, metaLogRow)
      art.conContrast = art.conContrastResult[["art.conContrast"]]
      metaLogRow = art.conContrastResult[["metaLogRow"]] # meta log row because tracks whether there was a warning
      
      # make initial data frame for logs with one row per trial
      # first argument can be any contrast except non-param
      factorSubsetLogDf = getFactorSubsetLogDfSmall(paramContrast, dataFrameNum, numResponsesPerCondition, between, condStrNumMap)
      # does not include nonParam since results are formatted differently
      contrastAndTestNameList = list(list("contrast" = omnibusart.contrast, "testName" = "ART_Omni"),
                                     list("contrast" = art.conContrast, "testName" = "ART_Con"),
                                     list("contrast" = paramContrast, "testName" = "Param"))
      for(contrastAndTestName in contrastAndTestNameList){
        factorSubsetLogDf = logResults(contrastAndTestName, factorSubsetLogDf, condStrNumMap)
      }
      
      # log nonParam results
      factorSubsetLogDf = logNonParamResults(nonParamContrast, factorSubsetLogDf, condStrNumMap)
      
    } # end between
    else{
      
      # random effects model
      fullFactorialARTModelRandom = fullFactorialARTModelList[["Random"]]
      
      # param
      paramContrastRandom = testParam(data, factorNames, factorSubset, between, modelType = "Random")
      
      # omnibus ART
      omnibusart.contrastRandomResult = testOmnibusART(data, fullFactorialARTModelRandom, factorSubset, metaLogRow)
      omnibusart.contrastRandom = omnibusart.contrastRandomResult[["result"]]
      metaLogRow = omnibusart.contrastRandomResult[["metaLogRow"]]
      
      # art.con
      art.conContrastResultRandom = testart.con(data, fullFactorialARTModelRandom, factorSubset, metaLogRow)
      art.conContrastRandom = art.conContrastResultRandom[["art.conContrast"]]
      metaLogRow = art.conContrastResultRandom[["metaLogRow"]] # adds warning flag to metalogrow
      
      # non param
      nonParamContrast = testNonParam(data, factorSubset, between)
      
      # Notes: within-ss logs are totally separate from between-ss logs
      #        fixed and random results go in same row so still get factorsubsetlogdf from omnibusart.contrastFixed
      # make initial data frame for logs with one row per trial
      factorSubsetLogDf = getFactorSubsetLogDfSmall(paramContrastRandom, dataFrameNum, numResponsesPerCondition, between, condStrNumMap)
      
      # does not include nonParam since results are formatted differently
      contrastAndTestNameList = list(list("contrast" = omnibusart.contrastRandom, "testName" = "ART_Omni_Random"),
                                     list("contrast" = art.conContrastRandom, "testName" = "ART_Con_Random"),
                                     list("contrast" = paramContrastRandom, "testName" = "Param_Random"))
      for(contrastAndTestName in contrastAndTestNameList){
        factorSubsetLogDf = logResults(contrastAndTestName, factorSubsetLogDf, condStrNumMap)
      }
      
      # log nonParam results
      factorSubsetLogDf = logNonParamResults(nonParamContrast, factorSubsetLogDf, condStrNumMap)
    } # end within
    
    # return
    list("factorSubsetLogDf" = factorSubsetLogDf, "metaLogRow" = metaLogRow)
  }
  
  # runs all types of tests
  # testandlogonefactorsubsetsmall only runs param, nonparam, and art.con
  # tests one factorSubset (ie combo of factors).
  # logs results
  # factorNames is all factors in data frame
  # fullFactorialARTModelList
    # if between-Ss list("Fixed" = art(Y ~ X1*X2*X3, data))
    # if within-Ss list("Fixed" = art(Y ~ X1*X2*X3 + (1|S), data), "Random" = art(Y ~ X1*X2*X3 + Error(S)))
  testAndLogOneFactorSubsetLarge = function(data, factorNames, fullFactorialARTModelList, factorSubset, dataFrameNum, 
                                       numResponsesPerCondition, between, condStrNumMap, metaLogRow){
    if(between){
      # full factorial ART model used by most contrasts
      fullFactorialARTModel = fullFactorialARTModelList[["Fixed"]]
      
      # leave this so when can verify that in the data we use for our analysis, only lmer causes errors.
      # run all contrasts
      art.conAlignedContrast = testart.conAligned(data, fullFactorialARTModel, factorSubset)
      artlmconAlignedContrast = testArtlmconAligned(data,fullFactorialARTModel, factorSubset)
      artlmconContrast = testArtlm.con(data, fullFactorialARTModel, factorSubset)
      paramContrast = testParam(data, factorNames, factorSubset, between)
      nonParamContrast = testNonParam(data, factorSubset, between)
      
      # testOmnibusART returns list
      omnibusart.contrastResult = testOmnibusART(data, fullFactorialARTModel, factorSubset, metaLogRow)
      omnibusart.contrast = omnibusart.contrastResult[["result"]]
      metaLogRow = omnibusart.contrastResult[["metaLogRow"]]
      
      # testart.con returns list
      art.conContrastResult = testart.con(data, fullFactorialARTModel, factorSubset, metaLogRow)
      art.conContrast = art.conContrastResult[["art.conContrast"]]
      metaLogRow = art.conContrastResult[["metaLogRow"]]
      
      # make initial data frame for logs with one row per trial
      factorSubsetLogDf = getFactorSubsetLogDfLarge(omnibusart.contrast, dataFrameNum, numResponsesPerCondition, between, condStrNumMap)
      # does not include nonParam since results are formatted differently
      contrastAndTestNameList = list(list("contrast" = omnibusart.contrast, "testName" = "ART_Omni"),
                          list("contrast" = art.conContrast, "testName" = "ART_Con"),
                          list("contrast" = art.conAlignedContrast, "testName" = "ART_Con_Aligned"),
                          list("contrast" = artlmconAlignedContrast, "testName" = "ARTlm_Con_Aligned"),
                          list("contrast" = artlmconContrast, "testName" = "ARTlm_Con"),
                          list("contrast" = paramContrast, "testName" = "Param"))
      for(contrastAndTestName in contrastAndTestNameList){
        factorSubsetLogDf = logResults(contrastAndTestName, factorSubsetLogDf, condStrNumMap)
      }
      
      # log nonParam results
      factorSubsetLogDf = logNonParamResults(nonParamContrast, factorSubsetLogDf, condStrNumMap)
      
    } # end between
    else{
      # fixed effects model
      fullFactorialARTModelFixed = fullFactorialARTModelList[["Fixed"]]
      
      art.conAlignedContrastFixed = testart.conAligned(data, fullFactorialARTModelFixed, factorSubset)
      artlmconAlignedContrastFixed = testArtlmconAligned(data,fullFactorialARTModelFixed, factorSubset)
      artlmconContrastFixed = testArtlm.con(data, fullFactorialARTModelFixed, factorSubset)
      paramContrastFixed = testParam(data, factorNames, factorSubset, between, modelType = "Fixed")
      
      # testOmnibusART returns list
      omnibusart.contrastFixedResult = testOmnibusART(data, fullFactorialARTModelFixed, factorSubset, metaLogRow)
      omnibusart.contrastFixed = omnibusart.contrastFixedResult[["result"]]
      metaLogRow = omnibusart.contrastFixedResult[["metaLogRow"]]
  
      # testart.con returns list
      art.conContrastResultFixed = testart.con(data, fullFactorialARTModelFixed, factorSubset, metaLogRow)
      art.conContrastFixed = art.conContrastResultFixed[["art.conContrast"]]
      metaLogRow = art.conContrastResultFixed[["metaLogRow"]]
      
      # random effects model
      fullFactorialARTModelRandom = fullFactorialARTModelList[["Random"]]
    
      art.conAlignedContrastRandom = testart.conAligned(data, fullFactorialARTModelRandom, factorSubset)
      artlmconAlignedContrastRandom = testArtlmconAligned(data,fullFactorialARTModelRandom, factorSubset)
      artlmconContrastRandom = testArtlm.con(data, fullFactorialARTModelRandom, factorSubset)
      paramContrastRandom = testParam(data, factorNames, factorSubset, between, modelType = "Random")
      
      # testOmnibusART returns list
      omnibusart.contrastRandomResult = testOmnibusART(data, fullFactorialARTModelRandom, factorSubset, metaLogRow)
      omnibusart.contrastRandom = omnibusart.contrastRandomResult[["result"]]
      metaLogRow = omnibusart.contrastRandomResult[["metaLogRow"]]
      
      # testart.con returns list
      art.conContrastResultRandom = testart.con(data, fullFactorialARTModelRandom, factorSubset, metaLogRow)
      art.conContrastRandom = art.conContrastResultRandom[["art.conContrast"]]
      metaLogRow = art.conContrastResultRandom[["metaLogRow"]]
  
      # not dependent on fixed or random model
      nonParamContrast = testNonParam(data, factorSubset, between)
      
      # Notes: within-ss logs are totally separate from between-ss logs
      #        fixed and random results go in same row so still get factorsubsetlogdf from omnibusart.contrastFixed
      # make initial data frame for logs with one row per trial
      factorSubsetLogDf = getFactorSubsetLogDfLarge(omnibusart.contrastFixed, dataFrameNum, numResponsesPerCondition, between, condStrNumMap)
      
      # does not include nonParam since results are formatted differently
      contrastAndTestNameList = list(list("contrast" = omnibusart.contrastFixed, "testName" = "ART_Omni_Fixed"),
                                     list("contrast" = art.conContrastFixed, "testName" = "ART_Con_Fixed"),
                                     list("contrast" = art.conAlignedContrastFixed, "testName" = "ART_Con_Aligned_Fixed"),
                                     list("contrast" = artlmconAlignedContrastFixed, "testName" = "ARTlm_Con_Aligned_Fixed"),
                                     list("contrast" = artlmconContrastFixed, "testName" = "ARTlm_Con_Fixed"),
                                     list("contrast" = paramContrastFixed, "testName" = "Param_Fixed"),
                                     list("contrast" = omnibusart.contrastRandom, "testName" = "ART_Omni_Random"),
                                     list("contrast" = art.conContrastRandom, "testName" = "ART_Con_Random"),
                                     list("contrast" = art.conAlignedContrastRandom, "testName" = "ART_Con_Aligned_Random"),
                                     list("contrast" = artlmconAlignedContrastRandom, "testName" = "ARTlm_Con_Aligned_Random"),
                                     list("contrast" = artlmconContrastRandom, "testName" = "ARTlm_Con_Random"),
                                     list("contrast" = paramContrastRandom, "testName" = "Param_Random"))
      for(contrastAndTestName in contrastAndTestNameList){
        factorSubsetLogDf = logResults(contrastAndTestName, factorSubsetLogDf, condStrNumMap)
      }
      
      # log nonParam results
      factorSubsetLogDf = logNonParamResults(nonParamContrast, factorSubsetLogDf, condStrNumMap)
    } # end within
    
    # return
    list("factorSubsetLogDf" = factorSubsetLogDf, "metaLogRow" = metaLogRow)
  }
  
  # metaLogRow input has all entries, but Mi and SDi for all conditions are NA
  # completes metaLogRow with Mi and Sdi for all conditions
  # logs metaLogRow
  # resultsLogFilename is passed in, but just replace word "results" with "meta"
  # Note: factorSubsets is not sorted
  finishAndLogMetaLogRow = function(data, factorSubsets, metaLogRow, condStrNumMap, resultsLogFilename, between){
    
    # iterate through each small factor subset
    for(factorSubset in factorSubsets){
      # remove other factors from data
      tempData = data %>% select(c(factorSubset, "Y"))
      # tested in excel, this does what i want it to.
      # get mean and sd of Y based on each combo of levels in factorSubset
      summary = ddply(tempData, factorSubset, summarise, mean=mean(Y), sd=sd(Y))
      
      # iterate through each combo of levels in factor subset and log mean and sd
      for(j in 1:nrow(summary)){
        
        currRow = summary[j,]
        
        # get value of each level
        currLevels = as.character(unlist((currRow[factorSubset])))
        # e.g currLevels = c("X1" = "a", "X2" = "b") -> "a,b"
        strCondition = paste(currLevels, collapse = ",")
        # lookup numeric code for condition
        numCondition = condStrNumMap[[strCondition]]
        
        # change Mi and SDi to correct mean for this condition
        currMean = currRow$mean
        metaLogRow[[paste("M", numCondition, sep="")]] = currMean
        
        currSd = currRow$sd
        metaLogRow[[paste("SD", numCondition, sep="")]] = currSd
      } # end iterate through conditions in one factorSubset
    } # end iterate through factorSubsets
    
    # if within subjects, want SDsOf and SDsProp at end
    if(!between){
      # metaLogRow is a list
      #metaLogRow = metaLogRow %>% plyr::mutate(tempOffsetSd = offsetSd, tempOffsetProp = offsetProp, offsetSd = NULL, offsetProp = NULL) %>% lplyr::rename(offsetSd = tempOffsetSd, offsetProp = tempOffsetProp)
      tempOffsetSd = metaLogRow[["offsetSd"]]
      tempOffsetProp = metaLogRow$offsetProp
      
      metaLogRow = metaLogRow %>% plyr::mutate(offsetSd = NULL, offsetProp = NULL)
      metaLogRow = c(metaLogRow, offsetSd = tempOffsetSd, offsetProp = tempOffsetProp)
    }
    
    # log metaRow
    metaLogFilename = str_replace_all(resultsLogFilename, "results", "meta")
    write.table(metaLogRow, metaLogFilename, sep = ",", col.names = !file.exists(metaLogFilename), row.names = FALSE, append = file.exists(metaLogFilename))
    
  } # end function
  
  # test this one data set
  # loop through every combination of factors in data and test all contrasts
  # meta Log
    # finish computing Mi and SDi for meta
    # write meta log to file
  # results log
    # add all results to logDf
    # write logDf to logfile when done
  analyzeOneDataSet <- function(data, dataFrameNum, between, resultsLogFilename, numResponsesPerCondition, condStrNumMap, metaLogRow){
    logDf = data.frame()
    # get all combos of all factors
    # first col of data is S, last col is Y, rest are factors
    # e.g. c("S", X1", "X2", "Y")
    colNames = colnames(data)
    # e.g. c("X1", "X2")
    factorNames = colNames[2:(length(colNames)-1)]
    # e.g. list("X1", "X2", "X1 X2")
    # note: this language isn't great. a condition is really (a,c) not (X1, X2).
    factorSubsets = powerset(factorNames)
    
    # list("Fixed" = art(), "Random" = art()). Random only there if Within-Ss
    fullFactorialARTModelList = makeFullFactorialARTModel(data, factorNames, between)
    # 
    for(factorSubset in factorSubsets){
      # data frame with results for all contrasts in factorSubset
      
      # testAndLogOneFactorSubsetSmall resturns list("factorSubsetLogDf" = factorSubsetLogDf, "metaLogRow" = metaLogRow)
      testAndLogOneFactorSubsetResult = testAndLogOneFactorSubsetSmall(data, factorNames, fullFactorialARTModelList, factorSubset, dataFrameNum,
                                                   numResponsesPerCondition, between, condStrNumMap, metaLogRow)
      factorSubsetLogDf = testAndLogOneFactorSubsetResult[["factorSubsetLogDf"]]
      metaLogRow = testAndLogOneFactorSubsetResult[["metaLogRow"]]
  
      # add factorSubsetLogDf to logDf
      logDf = rbind(logDf, factorSubsetLogDf)
      
    }
    
    rownames(logDf) = NULL
    
    # write log df to log file
    write.table(logDf, resultsLogFilename, sep = ",", col.names = !file.exists(resultsLogFilename), row.names = FALSE, append = file.exists(resultsLogFilename))
    
    # only one meta log row per data frame
    finishAndLogMetaLogRow(data, factorSubsets, metaLogRow, condStrNumMap, resultsLogFilename, between)
    
  } # end analyze data frame
  
  ####### END ANALYZE DATA #######
  
  ####### MAKE DATA #######
  
  # each data point is pulled from y ~ p(mu0, s0)
  # mu0 = g(alpha + beta*x +...)
  # latent mean for one condition is alpha + beta*x for that condition.
  # latent mean is either set at 0, or pulled from some normal distr N(mu1, s1).
  # the sd returned is s1.
  # currMeanSeed is var name for latent mean. ms_ss is var name for sd
  # RETURNS list("latentMean" = currMeanSeed, "latentMeanSd" = ms_ss, "trueSd" = sdSeed, "metaLogRow" = metaLogRow)
  # NOTE: double exponential distrs only have one param, but still return a value for trueSd
  # for consistency, even though it doesn't get used.
  getLatentMeanAndSd = function(infoFilename, metaLogRow, conditionNum){
    # set seed values
      currMeanSeed = 0
      # Sd for distr currMeanSeed was pulled from
      ms_ss = 1
      # sd that will be used to pull y values
      sdSeed = 1
      
      # normal, lognormal, exponential, cauchy, 3df t, double exponential with same mean and same sd
      if(testType %in% c(0,10,20,30,40,50)){
        currMeanSeed = 0
        sdSeed = 1
  
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = %s, sdSeed = %s", currMeanSeed, sdSeed),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
      
      # normal, lognormal, cauchy, 3 df t, double exponential with same mean seeds and diff sd seeds
      if(testType %in% c(1, 11, 31, 41, 51)){
        currMeanSeed = 0
        ss_ms = 0
        ss_ss = 1
        sdSeed = abs(rnorm(1, mean = ss_ms, sd=ss_ss))
        
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = %s, sdSeed = abs(rnorm(1, mean = %s, sd=%s))",
                      currMeanSeed, ss_ms, ss_ss),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
  
      # normal, lognormal, cauchy, double exponential, 3dft distr with diff means but same scale
      if(testType %in% c(101, 111, 132, 151, 141)){
        ms_ms = 0
        ms_ss = 1
        
        currMeanSeed = rnorm(1, mean=ms_ms, sd=ms_ss)
        sdSeed = 1
  
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = rnorm(1, mean=%s, sd=%s), sdSeed = %s",
                      ms_ms, ms_ss, sdSeed),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
  
      # normal, lognormal, cauchy, exponetial, double exponential, 3 df t, with diff mean seeds and diff sd seeds
      if(testType %in% c(100, 110, 120, 130, 140, 150)){
        ms_ms = 0
        ms_ss = 1
        ss_ms = 0
        ss_ss = 1
        currMeanSeed = rnorm(1, mean=ms_ms, sd=ms_ss)
        sdSeed = abs(rnorm(1, mean = ss_ms, sd=ss_ss))
  
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = rnorm(1, mean=%s, sd=%s), sdSeed = abs(rnorm(1, mean = %s, sd=%s))",
                      ms_ms, ms_ss, ss_ms, ss_ss),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
  
      # cauchy with diff mean seeds and sd seeds
      # mean seed and sd seed pulled from normal with mean 10, sd 2
      if(testType %in% c(131)){
        ms_ms = 10
        ms_ss = 5
        ss_ms = 0
        ss_ss = 1
        currMeanSeed = rnorm(1, mean=ms_ms, sd=ms_ss)
        sdSeed = abs(rnorm(1, mean = ss_ms, sd = ss_ss))
  
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = rnorm(1, mean=%s, sd=%s), sdSeed = abs(rnorm(1, mean = %s, sd=%s))",
                      ms_ms, ms_ss, ss_ms, ss_ss),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
  
      # 3df t dist with different noncentrality parameters and same scale, just bigger means
      if(testType == 142){
        ms_ms = 10
        ms_ss = 2
        currMeanSeed = rnorm(1, mean=ms_ms, sd=ms_ss)
        sdSeed = 1
  
        if(!infoFileWritten){
          cat(sprintf("currMeanSeed = rnorm(1, mean=%s, sd=%s)", ms_ms, ms_ss),
              file=infoFilename,append=TRUE,sep="\n")
        }
      }
  
  
      # log seed values
      # add mean seed to logs
      currMeanSeedName = paste("Ms", conditionNum, sep="")
      metaLogRow[[currMeanSeedName]] = currMeanSeed
      # add sd seed to logs
      currSdSeedName = paste("SDs", conditionNum, sep="")
      metaLogRow[[currSdSeedName]] = sdSeed
      
      # add NA entries for Mi and SDi, will add actual values later
      # NOTE: need to do this here so order of metaLogRow is correct
      metaLogRow[[paste("M", conditionNum, sep="")]] = NA
      metaLogRow[[paste("SD", conditionNum, sep="")]] = NA
      
      
      return(list("latentMean" = currMeanSeed, "latentMeanSd" = ms_ss, "trueSd" = sdSeed, "metaLogRow" = metaLogRow))
  }
  
  # add random offsets for random intercepts to latent mean column in data
  # Returns: list("data" = data, "metaLogRow" = metaLogRow)
  addRandomIntercepts = function(data, metaLogRow, latentMeanSd){
    # if random intercepts: for each participant get one value from N(0, p*latentMeanSd)
    # where p is some proportion
    # add this value to all latent means for that participant
    
    # get proportion
    pVec = c(0.1, 0.5, 0.9)
    # use data frame number to choose proportion
    dataFrameNum = metaLogRow[["Data_Set"]]
    p = pVec[dataFrameNum %% length(pVec) + 1]
    
    # get s
    s = p*latentMeanSd
    
    # get offset for each participant
    numSubjects = max(data[["S"]])
    interceptOffsets = rnorm(numSubjects, mean = 0, sd = s)
    
    # log s and p
    metaLogRow[["offsetProp"]] = p
    metaLogRow[["offsetSd"]] = s
    
    # add offset for each participant to all their latent means
    data["latentMean"] = adply(data, 1, function(x) x[["latentMean"]] + interceptOffsets[[x[["S"]]]], .expand=FALSE, .id=NULL)
    
    # return
    return(list("data" = data, "metaLogRow" = metaLogRow))
  }
  
  # applies link function to latent mean. puts result in column named mu
  applyLinkFunction = function(data){
    #testTypes = c(0,10,20,30,100,101,110,111,120,130,131,132,141,142,150,151)
    
    # normalDistrs = c(0,100,101) # linear
    # logNormalDistrs = c(10,110,111) # linear (rlnorm does it for us)
    # exponentialDistrs = c(20,120) # log link (ie use exp) and 1/mu = rate
    # cauchyDistrs = c(30,130,131,132) # linear
    # threeDfTDistrs = c(140,141,142) # linear
    # doubleExponentialDistrs = c(150,151) # linear
    
    # mu = exp(latentMean). then put 1/mu in mu column because exponential takes in rate
    if(testType %in% exponentialDistrs){
      # log transform latent mean
      data = data %>% mutate(mu = 1.0/exp(latentMean), latentMean = NULL)
    }
    
    # linear link
    else{
      data = data %>% mutate(mu = latentMean, latentMean = NULL)
    }
    data
  }
  
  
  # data has rows trueSd and mu
  # make new row Y = distr(mu, trueSd)
  # where distr depends on testType
  # Returns: data with Y  values in each row
  getNewY = function(data){
    # defined globally
    # normalDistrs = c(0,100,101) DONE
    # logNormalDistrs = c(10,110,111) DONE
    # exponentialDistrs = c(20,120) DONE
    # cauchyDistrs = c(30,130,131,132) DONE
    # threeDfTDistrs = c(40,140,141,142)
    # doubleExponentialDistrs = c(150,151)
  
    if(testType %in% normalDistrs){
      data[["Y"]] = apply(data, 1, function(x) rnorm(1, mean = x[["mu"]], sd = x[["trueSd"]]))
      if(!infoFileWritten){
        cat("sample distr: rnorm(1, mean = mu, sd = trueSd)",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
    
    if(testType %in% logNormalDistrs){
      data[["Y"]] = apply(data, 1, function(x) rlnorm(1, meanlog = x[["mu"]], sdlog = x[["trueSd"]]))
      if(!infoFileWritten){
        cat("sample distr: rlnorm(1, meanlog = mu, sdlog = trueSd)",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
   
    # data[["mu"]] is really rate for exp distrs.
    # data[["trueSd"]] has default value 1.
    if(testType %in% exponentialDistrs){
      data[["Y"]] = apply(data, 1, function(x) rexp(1, rate = x[["mu"]]))
      if(!infoFileWritten){
        cat("sample distr: rexp(1, rate = mu)\n Note: column name in data frame was mu, but it was the computed rate",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
    
    if(testType %in% cauchyDistrs){
      data[["Y"]] = apply(data, 1, function(x) rcauchy(1, location = x[["mu"]], scale = x[["trueSd"]]))
      if(!infoFileWritten){
        cat("sample distr: rcauchy(1, location = mu, scale = trueSd)",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
    
    if(testType %in% threeDfTDistrs){
      data[["Y"]] = apply(data, 1, function(x) rTF(1, mu = x[["mu"]], sigma = x[["trueSd"]], nu = 3))
      if(!infoFileWritten){
        cat("sample distr: rTF(1, mu = mu, sigma = trueSd, nu = 3)",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
    
    if(testType %in% doubleExponentialDistrs){
      data[["Y"]] = apply(data, 1, function(x) rdexp(1, location = x[["mu"]], scale = x[["trueSd"]]))
      if(!infoFileWritten){
        cat("sample distr: rdexp(1, location = mu, scale = trueSd)",
            file=infoFilename,append=TRUE, sep="\n")
      }
    }
    
    # remove mu and trueSd columns
    data = data %>% mutate(mu = NULL, trueSd = NULL)
  
    return(data)
  }
  
  # make data frame for this one data set
  # RETURNS: data frame for single data set
  #          metadata for log file (metaLogRow)
  makeOneDataSet = function(infoFilename, dataFrameNum, between, testType, conditionsDf, numConditions, numResponsesPerCondition, condStrNumMap){
    data = data.frame()
    
    # temp row has log info that will be used by every row in log
    metaLogRow = list("Data_Set" = dataFrameNum, "ART_Con_Warning" = 0, "Omnibus_ART_Warning" = 0, "Num_Responses_per_Condition" = numResponsesPerCondition, "Random_Seed" = randomSeed) # Starting row of log data frame.
    metaLogRow[["BW"]] = if(between) "B" else "W"
    # add B or W to metaLogRow
    
    # latent mean sd is same for all conditions in one df (actually same for all conditions in one test type)
    # will update with different value from getLatentMeanAndSd if changed
    latentMeanSd = 1
    
  # get y values for each condition and add to data.
    for(i in 1:numConditions){
      # replicate this condition correct number of times
      currConditionRow = as.vector(conditionsDf[i,])
      currConditionDf = do.call("rbind", replicate(numResponsesPerCondition, currConditionRow, simplify = FALSE))
      
      # add subject ids
      if(between){
        startId = numResponsesPerCondition*(i-1) + 1
        endId = startId + numResponsesPerCondition - 1 # inclusive on end
        S = startId:endId
        currConditionDf = cbind(S, currConditionDf)
      } else {
        S = 1:numResponsesPerCondition
        currConditionDf = cbind(S, currConditionDf)
      }
      
      rownames(currConditionDf) = NULL # need to reset row index
      
      # check that condition's number lines up with its number in condStrNumMap
      # NOTE: since condStrNumMap was created from conditionsDf, the value for a condition in condStrNumMap should be the same as
        # the index of the row in conditionsDf that contains that condition. 
        # assert this is true in each iteration.
      currConditionStr = paste((unlist(currConditionRow, use.names = FALSE)), collapse = ",")
      correctConditionNumber = condStrNumMap[[currConditionStr]]
      if(i != correctConditionNumber){
        stop("Condition order in makeOneDataSet does not match up with condition order in condStrNumMap")
      }
    
      
      # add Numi and Stri to to metaLogRow
      metaLogRow[[paste("Num", i, sep="")]] = i
      # log strings names with space since csv
      metaLogRow[[paste("Str", i, sep="")]] = paste((unlist(currConditionRow, use.names = FALSE)), collapse = " ")
      
      #  for condition get latent mean and sd for latent mean
      # e.g. latent mean pulled from y = normal(m, s). return y and s.
  
      latentMeanAndSd = getLatentMeanAndSd(infoFilename, metaLogRow, i)
      latentMean = latentMeanAndSd[["latentMean"]] # argument to link function (maybe adding offsets). used to be currMeanSeed
      trueSd = latentMeanAndSd[["trueSd"]] # will get response ~ p(mu, trueSd) used to be sdSeed.
      metaLogRow = latentMeanAndSd[["metaLogRow"]]
      
      # latent mean sd is same for all conditions. 
      # latentMeanSd var is defined outside the loop. Update its value here.
      if(i == 1){
        latentMeanSd = latentMeanAndSd[["latentMeanSd"]] # latentMean = Normal(m, latentMeanSd)
      }
      
      # trueSd is one value, but it is same for all rows in this condition.
      # temporarily add it to data so can use it when get new y
      trueSdVec = rep(trueSd, numResponsesPerCondition)
      currConditionDf[["trueSd"]] = trueSdVec
      # latent mean is one value, but is same for all rows in this condition.
      # vector c(latentMean, latentMean,...)
      latentMeanVec = rep(latentMean, numResponsesPerCondition)
      currConditionDf[["latentMean"]] = latentMeanVec
      
      infoFileWritten <<- TRUE
      
      # add currConditionDf to data
      data = rbind(data, currConditionDf)
    }
    
    # WITHIN SUBJECTS could have just random offsets, or random offsets and random slopes.
    if(!between){
      # Random offsets
      # if within subjects, each participant needs a random offset to add to all its latent means
      # this gives each participant a different random intercept.
      randomInterceptResult = addRandomIntercepts(data, metaLogRow, latentMeanSd)
      metaLogRow = randomInterceptResult[["metaLogRow"]]
      data = randomInterceptResult[["data"]]
    }
    
    
    # apply link function
    data = applyLinkFunction(data)
    
    # get new y
    data = getNewY(data)
    
    rownames(data) = NULL
    
    # columns of data are lists right now. Need to convert to traditional data frame
    data <- as.data.frame(lapply(data, unlist))
    # set factor and subject columns of data to factors
    data[,1:(ncol(data)-1)] <- lapply(data[,1:(ncol(data)-1)], as.factor)
    
    return(list("data" = data, "metaLogRow" = metaLogRow))
  }
  
  # make data frame with one row per condition
  # to be used for
  # RETURNS: data frame with one row pers condition
  makeConditionsDf <- function(numFactors, numLevelsPerFactor){
    
    levels = list() #c[0] will be a vector with all levels of factor 1 etc
    asciiA = 97
    for(factor in 0:(numFactors-1)){
      currLevels = list()
      # make list of str levels for each factor
      for(level in 0:(numLevelsPerFactor-1)){
        currLevelInt = asciiA + (factor*numLevelsPerFactor) + level
        currLevelStr = intToUtf8(currLevelInt)
        currLevels[[(level + 1)]] = currLevelStr
      } # end loop through levels for one factor
      # add list of str levels for one factor to levels list
      levels[[factor+1]] = currLevels
    }# end loop through factors
    
    # make data frame with each row as a new condition
    conditionsDf = expand.grid(levels)
    
    # fix colnames of conditionsDf
    colNames = c()
    for(i in 1:numFactors){
      factorStrName = paste("X", i, sep="")
      colNames = c(colNames, factorStrName)
    }
    colnames(conditionsDf) = colNames
    return(conditionsDf)
  }
  
  # metaLogRow has Data_Set, Num_Responses_per_Condition, BW, 
  # Msi, SDsi for newY conditions
  # Mi, SDi for newY conditions but values are NA
  
  # add Numi, Stri, for non newY conditions 
  # add entries for Mi and SDi for non newY conditions but values are NA
  appendToMetaLogRow = function(metaLogRow, condStrNumMap, numConditions){
    # sort condStrNumMap
    sortedCondStrNumMap = sort(values(condStrNumMap))
    # get rid of first numConditions elements in sorted map
    prunedCondStrNumMap = tail(sortedCondStrNumMap, (length(sortedCondStrNumMap)-numConditions))
    # iterate through remaining conditions, log Numi, Stri, Mi, Sdi
    for(i in 1:length(prunedCondStrNumMap)){
      # log Numi and Stri
      Numi = prunedCondStrNumMap[[i]]
      Stri = names(prunedCondStrNumMap)[[i]]
      correctLogIndex = numConditions + i # start from one after conditions that were logged earlier
      metaLogRow[[paste("Num", correctLogIndex, sep="")]] = Numi
      metaLogRow[[paste("Str", correctLogIndex, sep="")]] = Stri
      
      # set mi and sdi but make them NA. 
      # then later when analyzing one factor subset, calculate means and sds and log them
      # NOTE: need to make them NA so order will be 
      # Numi, Stri, Mi, SDi, Num(i+1), Str(i+1), M(i+1), SD(i+1), ...
      metaLogRow[[paste("M", correctLogIndex, sep="")]] = NA
      metaLogRow[[paste("SD", correctLogIndex, sep="")]] = NA
    }
    metaLogRow
  }
  
  ####### END MAKE DATA #######
  
  # create and analyze one data frame
  # log results
  makeAndAnalyzeOneDataSet <- function(resultsLogFilename, infoFilename, dataFrameNum, between, testType,
                                      conditionsDf, condStrNumMap){
    numResponsesPerConditionOptions = c(8,16,24,32,40) # will pick one of these options
    # randomly pick num responses per condition from options
    numResponsesPerCondition = sample(numResponsesPerConditionOptions, 1)
    
    numConditions = nrow(conditionsDf) # number of conditions that newY will be pulled for
    makeOneDataSetResult = makeOneDataSet(infoFilename, dataFrameNum, between, testType, conditionsDf, numConditions, numResponsesPerCondition, condStrNumMap)
    data = makeOneDataSetResult[["data"]]
    metaLogRow = makeOneDataSetResult[["metaLogRow"]]
    
    # Now has all info except Mi and SDi.
    metaLogRow = appendToMetaLogRow(metaLogRow, condStrNumMap, numConditions)
    
    analyzeOneDataSet(data, dataFrameNum, between, resultsLogFilename, numResponsesPerCondition, condStrNumMap, metaLogRow) # test data in "data" dataframe
  } # end makeAndAnalyzeOneDataSet)
  
  
  ####### END MAKE AND ANALYZE DATA #######
  
  ####### MAIN #######
  
  # for each subset of factors, each combination of levels with exactly 
  # one level from each factor is assigned a unique number. key is str, value is num
  # e.g. Str2 = "a c"-> condStrNumMap[["a,c"]] = 2
  # Note: numbers are assigned to conditions in descending order by number of factors involved
  #   i.e. conditions with levels from 3 factors have smaller numbers than conditions with levels from 2 factors
  # clear and reset here
  getCondStrNumMap = function(conditionsDf, numFactors, numLevelsPerFactor){
    condStrNumMap = hash()
    
    # get all factor subsets (powerset of factor names)
    # sorted in descending order by number of factors in subset
    factorNames = colnames(conditionsDf)
    # list of vectors. elems of each vector are names of factors in subset
    allFactorSubsetsUnsorted = powerset(factorNames)
    # in decreasing order by number of factors in subset (ie subsets with most factors first)
    # need conditions is largest number of factors to have smallest numbers since need to log seeds
    allFactorSubsets= allFactorSubsetsUnsorted[order(sapply(allFactorSubsetsUnsorted,length),decreasing=T)]
    
    # iterate through all factor subsets
    currId = 1 # id for current condition
    for(factorSubset in allFactorSubsets){
      # get subset of conditionsDf only containing cols from factorSubset
      subsetConditionsDf = conditionsDf[factorSubset]
      # get rid of duplicate rows
      # every condition with factors in factorSubset is represented by exactly one row in subsetConditionsDfUnique 
      subsetConditiondDfUnique = unique(subsetConditionsDf)
      
      # add entry to condStrNumMap for each row in subsetConditiondDfUnique
      # e.g. row = c("a", "b"), key = "a,b". value = next available value
      for(row in 1:nrow(subsetConditiondDfUnique)){
        # get row
        currRow = subsetConditiondDfUnique[row,]
        # make key
        currKey = paste((unlist(currRow, use.names = FALSE)), collapse = ",")
        # add entry
        condStrNumMap[[currKey]] = currId
        # increment id
        currId = currId + 1
      } # end iterate through conditions
    } # end iterate through factor subsets
    condStrNumMap # return
  } # end function
  
  # makes directory for test type and returns path to dir
  getTestTypeDir = function(parentDirPath, testType){
    # make dir
    dirName = paste(parentDirPath, testType, sep="/")
    dir.create(dirName)
    
    # return
    dirName
  }
  
  # makes directory for design and returns path to dir
  getDesignDirPath = function(parentDirPath, numFactors, numLevelsPerFactor){
    # e.g. numFactors = 2, numLevelsPerFactor = 3 -> "3x3"
    designStr = paste(rep(numLevelsPerFactor, numFactors), collapse="x")
    dirName = paste(parentDirPath, designStr, sep="/")
    # make dir
    dir.create(dirName)
    
    # return
    dirName
  }
  
  # creates string name for path to info file and returns it
  getInfoFilename <- function(parentDirPath, testType){
    infoFilename = paste(parentDirPath, "/", "info-", testType, ".txt", sep="") # full file name
    infoFilename
  }
  
  # makes directory for between/within and returns path to dir
  getBetweenDirPath = function(parentDirPath, between){
    betweenStr = if(between) "between" else "within"
    dirName = paste(parentDirPath, "/", betweenStr, sep="")
    
    # create dir
    dir.create(dirName)
    
    # return
    dirName
  }
  
  # Note: this gets called after get info file name (since info file same for between and and within, but log separate)
  getResultsLogFilename = function(parentDirPath, numFactors, numLevelsPerFactor, between, testType){
    # e.g. numFactors = 3, numLevelsPerFactor = 2 -> 2x2x2
    levelsStr = paste(replicate(numFactors, numLevelsPerFactor), collapse="x")
    # between = TRUE -> "between", between = FALSE -> "within"
    betweenStr = if (between) "between" else "within"
    
    # testType is global
    resultsLogFilename = paste(parentDirPath, "/results-log-", levelsStr,  "-type-", testType, "-", betweenStr, ".csv", sep="") # full file name
    resultsLogFilename
  }
  
  #### NEW MAINLOOP
  mainLoop<- function(testType, timeStampStr, testDirPath){
    
    # make directory for logs for this test type. return path to dir.
    testTypeDirPath = getTestTypeDir(testDirPath, testType)
    
    # get info file path.
    infoFilename = getInfoFilename(testTypeDirPath, testType)
    infoFileWritten <<- FALSE
    
    designs = list(list("numFactors" = 2, "numLevelsPerFactor" = 2),
                   list("numFactors" = 3, "numLevelsPerFactor" = 2),
                   list("numFactors" = 2, "numLevelsPerFactor" = 3)
    )
  
    for(design in designs){
      numFactors = design[["numFactors"]]
      numLevelsPerFactor = design[["numLevelsPerFactor"]]
      
      # make directory for logs for this design. return path to dir.
      designDirPath = getDesignDirPath(testTypeDirPath, numFactors, numLevelsPerFactor)
      
      # data frame with 1 row per condition. same for all data sets with current design.
      conditionsDf = makeConditionsDf(numFactors, numLevelsPerFactor)
      
      # e.g. Str2 = "a c"-> condStrNumMap[["a,c"]] = 2
      # Note: use comma in R since that's how contrasts print
      #       ues space when logging since don't want to use comma in csv
      # Note: numbers are assigned to conditions in descending order by number of factors involved
      #   i.e. conditions with levels from 3 factors have smaller numbers than conditions with levels from 2 factors
      condStrNumMap = getCondStrNumMap(conditionsDf, numFactors, numLevelsPerFactor)
      
      for(between in c(TRUE, FALSE)){
        # reset since separate logs for between and within
        dataFrameNum = 1 # local
        trialNum <<- 1 # global
        
        # between subjects logs and within subjects logs in different directories
        # this is the path to that directory
        betweenDirPath = getBetweenDirPath(designDirPath, between)
        
        # log file name
        resultsLogFilename = getResultsLogFilename(betweenDirPath, numFactors, numLevelsPerFactor, between, testType)
        
        for(i in 1:numDataSets){
          if((i-1)%%25 == 0){ # print every 25 data sets
            cat("test type")
            print(testType)
            cat("design")
            print(paste(rep(numLevelsPerFactor, numFactors), collapse='x'))
            cat("between")
            print(between)
            cat("test #")
            print(i)
            cat("\n")
          }
          makeAndAnalyzeOneDataSet(resultsLogFilename, infoFilename, dataFrameNum, between, testType,
                                   conditionsDf, condStrNumMap) # one data data frame
          dataFrameNum  = dataFrameNum + 1
        } # end iterate through all data sets
      } # end between and within
    }
  }
  #### END NEW MAINLOOP
  
  # makes directory for results from this test and return path
  getTestDir = function(timeStampStr){
    # logs dir doesn't already exist, make it
    if(!dir.exists("logs")){
      dir.create("logs")
    }
    
    # make dir
    # it won't already exist since based on system time
    dirName = paste("logs/test-", timeStampStr, sep="")
    dir.create(dirName)
    
    # return
    dirName
  }
  
  # globals
  trialNum = 1
  numDataSets = 1000
  randomSeed = 1
  infoFileWritten = FALSE # new info file for every testType x Design
  warningCt = 0 # number of times we get a warning
  
  # useful for link functions and pulling from distrs
  normalDistrs = c(0,1,100,101)
  logNormalDistrs = c(10,11,110,111)
  exponentialDistrs = c(20,120)
  cauchyDistrs = c(30,31,130,131,132)
  threeDfTDistrs = c(40,41,140,141,142)
  doubleExponentialDistrs = c(50,51,150,151)
  
  # NOTE no 21 or 121 because exponential only has one param
  testTypes = c(normalDistrs,logNormalDistrs,exponentialDistrs,cauchyDistrs,threeDfTDistrs,doubleExponentialDistrs)
  testType = 0 # make explicitly global since running inside function now
  # Second arg doesn't matter. we're not doing ordered contrasts.
  options(contrasts=c("contr.sum", "contr.sum"))
  
  # get testTypes from command line
  # if no args supplied, use default testTypes above (is all test types)
  args=(commandArgs(TRUE))
  if(length(args) > 0){
    testTypes = sapply(parse(text = args), eval)
  }
  
  runTest = function(){
    emm_options(sep=",")
    # results dir path and time for this test
    timeStamp = as.character(Sys.time()) # timestamp for log file name
    timeStampStr = str_replace_all(timeStamp, c(":" = "-", " " = "-"))
    # path to directory for results from this test
    testDirPath = getTestDir(timeStampStr)
    
    for(localTestType in testTypes){
      set.seed(randomSeed) # set random seed at start of each test type
      testType <<- localTestType
      suppressMessages(mainLoop(localTestType, timeStampStr, testDirPath))
    }
  }
  
  runTest()
  
