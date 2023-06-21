fitPOM_fn <- function(endpointName, 
                      predictorNames_vec, 
                      offsetName=NULL, 
                      design_dat, 
                      adverseEndpoint_vec=nlevels(design_dat[[endpointName]]), 
                      analysisLabel=NULL, 
                      subgroupLabel=NULL, 
                      outputPath, 
                      interactionLabel=NULL, 
                      interactionName=NULL, 
                      boot_size=NULL
                      )
  {
  # Function fits a Proportional Odds regression model to a data-set including a 
  # count variable or ordinal factor (response), a set of predictors (covariates) 
  # and a binary exposure ('treat'), and computes tables of customary diagnostic 
  # statistics and model inferences.
  # 
  # 
  # Arguments:
  # 
  # endpointName        : character string denoting the response variable's name 
  #                       (already an ordered factor)
  # predictorNames_vec  : character vector denoting explanatory variables' names 
  # offsetName          : character string denoting the offset variable's name; 
  #                       if NULL, no offset term is included
  # design_dat          : design data-frame
  # adverseEndpoint_vec : integer vector of level indices identifying adverse outcome 
  #                       levels;
  # analysisLabel       : character string denoting the analysis type label
  # subgroupLabel       : character string denoting the subgroup analysis label
  # outputPath          : character string denoting the path to the output 
  #                       destination folder
  # interactionLabel    : character string denoting the "treat" interaction label
  # interactionName     : character string denoting the name of an explanatory  
  #                       variable interacting with "treat"
  # boot_size           : integer scalar denoting the bootstrap sample size; if
  #                       NULL, no bootstrap is carried out
  # 
  # 
  # Output: 
  # 
  # List of named function arguments and ASCII (.csv) tables of diagnostic 
  # statistics and model inferences.
  # 
  # N.B.: despite apparently successful attempts, coded in foreign packages, to 
  # extract p-values from polr() model objects, original MASS authors caution against 
  # their use, being the Wald-type tests they are derived from inappropriate and 
  # based on potentially misleading asymptotic theory.
  
  
  # endpointName <- "re78_class"
  # predictorNames_vec <- c("age", "educ", "race", "nodegree", "married", "re74", "re75")
  # offsetName <- NULL
  # design_dat <- matched.dat
  # adverseEndpoint_vec <- seq.int(nlevels(design_dat[[endpointName]]))[-nlevels(design_dat[[endpointName]])]
  # analysisLabel <- "full"
  # subgroupLabel <- "max"
  # outputPath <- file.path(project.dir,
  #                         "R/DA/Projects/Methods tool-box")
  # interactionLabel <- NULL
  # interactionName <- NULL
  # boot_size <- 2e3
  
  
  ##############
  ## Preamble ##
  ##############
  
  require(boot)  # Load library required for bootstrap resampling
  
  library(MASS)  # Load library enabling "Modern Applied Statistics with S" routines
  
  require(multcomp)  # Load library required for multiple inferences
  
  require(nnet)  # Load library required for Log-Linear modelling
  
  require(sandwich)  # Load library required for model-robust covariance estimation
  
  
  stopifnot(is.character(endpointName), 
            is.character(predictorNames_vec), 
            is.null(offsetName) | is.character(offsetName), 
            is.data.frame(design_dat), 
            all(adverseEndpoint_vec %in% seq.int(nlevels(design_dat[[endpointName]]))), 
            is.null(analysisLabel) | is.character(analysisLabel), 
            is.null(subgroupLabel) | is.character(subgroupLabel), 
            is.null(interactionLabel) | is.character(interactionLabel), 
            is.null(interactionName) | is.character(interactionName), 
            boot_size > 0 & boot_size %% 1 == 0
            )
  
  if(! endpointName %in% names(design_dat)) 
    stop("The outcome variable is not in the analysis data-set.")
  
  if(! "treat" %in% names(design_dat)) 
    stop("The exposure variable is not in the analysis data-set.")
  
  if(any(sort(unique(design_dat$treat)) != 0:1)) 
    stop("The exposure variable is not coded as binary.")
  
  if(is.character(offsetName)) 
    if(offsetName %in% names(design_dat)) 
      stop("The offset variable is not in the analysis data-set.")
  
  if(! is.ordered(design_dat[[endpointName]])) 
    stop("The outcome variable is not an ordered factor.")
  
  if(! any(c(predictorNames_vec, interactionName) %in% names(design_dat))) 
    stop("At least one of the predictors to the GLM models is not in the analysis data-set.")
  
  if(is.character(outputPath)) 
    if(! dir.exists(outputPath)) 
      dir.create(outputPath, recursive=TRUE)
  
  
  ########################
  ## Model formulations ##
  ########################
  
  design_dat <- droplevels(design_dat)  # Drop unused factor levels
  
  
  predictorNames_vec <- predictorNames_vec[sapply(design_dat[predictorNames_vec], 
                                                  FUN=function(vec) length(unique(vec)) > 1
                                                  )]  # Discard singular explanatory variables
  
  
  adj_frm <- as.formula(paste(endpointName, 
                              paste(c(if(is.character(interactionName)) 
                                paste("treat", interactionName, 
                                      sep=" * "
                                      ) else 
                                        "treat", 
                                setdiff(predictorNames_vec, 
                                        y=if(is.character(interactionName)) 
                                          interactionName
                                        )
                                ), 
                                collapse=" + "
                                ), 
                              sep=" ~ "
                              )
                        )  # Set formula for adjusted POM of response with explanatory variables
  
  
  unadj_frm <- update(adj_frm, 
                      new=paste(".", 
                                paste(".", 
                                      paste("-", setdiff(predictorNames_vec, 
                                                         y=if(is.character(interactionName)) 
                                                           interactionName
                                                         ), 
                                            collapse=" "
                                            ), 
                                      sep=" "
                                      ), 
                                sep=" ~ "
                                )
                      )  # Set formula for unadjusted POM of response with explanatory variables and offset term
  
  
  ###################
  ## Model fitting ##
  ###################
  
  pom_adj_fit <- tryCatch(polr(if(is.null(offsetName)) 
    adj_frm else 
      update(adj_frm, new=. ~ . + offset(log(get(offsetName)))), 
    data=design_dat, 
    # control=list(maxit=2e2, trace=FALSE), 
    Hess=TRUE#, model=FALSE, method="logistic"  # Specifying more options than necessary leads to incompatibilities with the 'multcomp' library
    ), 
    error=function(err) NULL
    )  # Fit adjusted POM to sample
  
  if(class(pom_adj_fit) == "polr" && pom_adj_fit$convergence == 0) 
    cat("Adjusted Proportional Odds regression model successfully fitted (AIC =", 
        round(pom_adj_fit$deviance + 2 * pom_adj_fit$edf, 
              digits=2
              ), 
        "\b) to", endpointName, "in the", analysisLabel, 
        "analysis of the", subgroupLabel, "sample\n\n", 
        sep=" "
        ) else   # Print to terminal AIC of fitted POM
          cat("Adjusted Proportional Odds regression model did not converge for", 
              endpointName, "in the", 
              analysisLabel, "analysis of the", 
              subgroupLabel, "sample\n\n", 
              sep=" "
              )  # Print progress message
  
  
  pom_unadj_fit <- tryCatch(polr(if(is.null(offsetName)) 
    unadj_frm else 
      update(unadj_frm, new=. ~ . + offset(log(get(offsetName)))), 
    data=design_dat, 
    # control=list(maxit=2e2, trace=FALSE), 
    Hess=TRUE#, model=FALSE, method="logistic"  # Specifying more options than necessary leads to incompatibilities with the 'multcomp' library
    ), 
    error=function(err) NULL
    )  # Fit unadjusted POM to sample
  
  if(class(pom_unadj_fit) == "polr" && pom_unadj_fit$convergence == 0) 
    cat("Unadjusted Proportional Odds regression model successfully fitted (AIC =", 
        round(pom_unadj_fit$deviance + 2 * pom_unadj_fit$edf, 
              digits=2
              ), 
        "\b) to", endpointName, 
        "in the", analysisLabel, 
        "analysis of the", subgroupLabel, 
        "sample\n\n", 
        sep=" "
        ) else   # Print to terminal AIC of fitted POM
          cat("Unadjusted Proportional Odds regression model did not converge for", 
              endpointName, "in the", analysisLabel, 
              "analysis of the", subgroupLabel, "sample\n\n", 
              sep=" "
              )  # Print progress message
  
  
  pom_fit_ls <- list(adjusted=if(class(pom_adj_fit) == "polr" && pom_adj_fit$convergence == 0)
    pom_adj_fit else 
      pom_unadj_fit, 
    unadjusted=if(class(pom_unadj_fit) == "polr" && pom_unadj_fit$converge == 0) 
      pom_unadj_fit
    )  # Set list by predictor adjustment type of fitted POM
  
  
  llm_adj_fit <- tryCatch(multinom(adj_frm, 
                                   data=design_dat, 
                                   offset=if(is.character(offsetName)) 
                                     log(get(offsetName)) * model.matrix(~ -1 + endpointName, data=design_dat), 
                                   maxit=1e3, Hess=TRUE, trace=FALSE#, model=FALSE, method="logistic"  # Specifying more options than necessary leads to incompatibilities with the 'multcomp' library
                                   ), 
                          error=function(err) NULL
                          )  # Fit adjusted LLM to sample
  
  if("multinom" %in% class(llm_adj_fit) && llm_adj_fit$convergence == 0) 
    cat("Adjusted Log-Linear regression model successfully fitted (AIC =", 
        round(llm_adj_fit$deviance + 2 * llm_adj_fit$edf, 
              digits=2
              ), 
        "\b) to", endpointName, 
        "in the", analysisLabel, 
        "analysis of the", subgroupLabel, 
        "sample\n\n", 
        sep=" "
        ) else   # Print to terminal AIC of fitted LLM
          cat("Adjusted Log-Linear regression model did not converge for", 
              endpointName, "in the", analysisLabel, 
              "analysis of the", subgroupLabel, "sample\n\n", 
              sep=" "
              )  # Print progress message
  
  
  llm_unadj_fit <- tryCatch(multinom(unadj_frm, 
                                     data=design_dat, 
                                     offset=if(is.character(offsetName)) 
                                       log(get(offsetName)) * model.matrix(~ -1 + endpointName, data=design_dat), 
                                     maxit=1e3, Hess=TRUE, trace=FALSE#, model=FALSE, method="logistic"  # Specifying more options than necessary leads to incompatibilities with the 'multcomp' library
                                     ), 
                            error=function(err) NULL
                            )  # Fit unadjusted LLM to sample
  
  if("multinom" %in% class(llm_unadj_fit)  && llm_unadj_fit$convergence == 0) 
    cat("Unadjusted Log-Linear regression model successfully fitted (AIC =", 
        round(llm_unadj_fit$deviance + 2 * llm_unadj_fit$edf, 
              digits=2
              ), 
        "\b) to", endpointName, 
        "in the", analysisLabel, 
        "analysis of the", 
        subgroupLabel, "sample\n\n", 
        sep=" "
        ) else   # Print to terminal AIC of fitted LLM
          cat("Unadjusted Proportional Odds regression model did not converge for", 
              endpointName, "in the", 
              analysisLabel, "analysis of the", 
              subgroupLabel, "sample\n\n", 
              sep=" "
              )  # Print progress message
  
  
  llm_fit_ls <- list(adjusted=if("multinom" %in% class(llm_adj_fit) && llm_adj_fit$convergence == 0) 
    llm_adj_fit else 
      llm_unadj_fit, 
    unadjusted=if("multinom" %in% class(llm_unadj_fit) && llm_unadj_fit$converge == 0) 
      llm_unadj_fit
    )  # Set list by predictor adjustment type of fitted LLM
  
  
  #####################
  ## Model diagnosis ##
  #####################
  
  diag_mod_ls <- sapply(names(pom_fit_ls), 
                        FUN=function(nm) 
                          tryCatch(
                            data.frame("Pearson's Chi-squared"=NA, 
                                       "Residual deviance"=pom_fit_ls[[nm]]$deviance, 
                                       "Residual d.o.f."=pom_fit_ls[[nm]]$df.residual, 
                                       "Pearson's Chi-squared p-value"=NA, 
                                       "Residual deviance p-value"=pchisq(pom_fit_ls[[nm]]$deviance, 
                                                                          df=pom_fit_ls[[nm]]$df.residual, lower.tail=FALSE
                                                                          ), 
                                       AIC=pom_fit_ls[[nm]]$deviance + 2 * pom_fit_ls[[nm]]$edf, 
                                       RMSE=sqrt(mean((prop.table(table(design_dat[[endpointName]])) - 
                                                         colMeans(pom_fit_ls[[nm]]$fitted.values)) ^ 2
                                                      )
                                                 ), 
                                       "POM vs LLM LRT"=-2 * (logLik(pom_fit_ls[[nm]]) - logLik(llm_fit_ls[[nm]])), 
                                       "POM vs LLM LRT p-value"=pchisq(-2 * (logLik(pom_fit_ls[[nm]]) - logLik(llm_fit_ls[[nm]])), 
                                                                              df=llm_fit_ls[[nm]]$edf - pom_fit_ls[[nm]]$edf, 
                                                                              lower.tail=FALSE
                                                                              ), 
                                       "% observed minimum"=mean(design_dat[[endpointName]] == levels(design_dat[[endpointName]])[1]), 
                                       "% predicted minimum"=mean(pom_fit_ls[[nm]]$fitted.values[, levels(design_dat[[endpointName]])[1]]), 
                                       t(setNames(unclass(table(design_dat$treat)), 
                                                  nm=paste(c("Control", "Intervention"), 
                                                           "subgroup size", 
                                                           sep=" "
                                                           )
                                                  )
                                         ),   # Compute and rename sample size by exposure.  N.B.: transposing required for building single-row table
                                       check.names=FALSE
                                       ), 
                            error=function(err) 
                              data.frame("Pearson's Chi-squared"=NA, 
                                         "Residual deviance"=NA, 
                                         "Residual d.o.f."=NA, 
                                         "Pearson's Chi-squared p-value"=NA, 
                                         "Residual deviance p-value"=NA, 
                                         AIC=NA, 
                                         RMSE=NA, 
                                         "POM vs LLM LRT"=NA, 
                                         "POM vs LLM LRT p-value"=NA, 
                                         "% observed minimum"=NA, 
                                         "% predicted minimum"=NA, 
                                         t(setNames(unclass(table(design_dat$treat)), 
                                                    nm=paste(c("Control", "Intervention"), 
                                                             "subgroup size", 
                                                             sep=" "
                                                             )
                                                    )
                                           ),   # Compute and rename sample size by exposure.  N.B.: transposing required for building single-row table
                                         check.names=FALSE
                                         )
                            ), 
                        simplify=FALSE
                        )  # Derive list by adjustment type of commonplace diagnostic regression statistics for fitted POM
  
  
  ##########################
  ## Bootstrap inferences ##
  ##########################
  
  if(is.numeric(boot_size)){
    risk_boot_fn <- function(data_dat, idx_vec, adj){
      prd.sng_vec <- which(sapply(data_dat[idx_vec, ], 
                                  FUN=function(vec) all(na.omit(vec)[-1] == na.omit(vec)[1]))
                           )  # Derive vector of singular predictor indices from bootstrap analysis data-frame
      
      adj_boot_frm <- update(adj_frm, 
                             new=paste(".", 
                                       paste(".", 
                                             paste(names(prd.sng_vec), collapse=" - "), 
                                             sep=" - "
                                             ), 
                                       sep=" ~ "
                                       )
                             )  # Set formula for adjusted POM of response with non-singular explanatory variables from bootstrap analysis data-frame
      
      
      pom_boot_fit <- tryCatch(polr(switch(match(adj, table=names(pom_fit_ls)), 
                                           if(is.null(offsetName)) 
                                             adj_boot_frm else 
                                               update(adj_boot_frm, new=. ~ . + offset(log(get(offsetName)))), 
                                           if(is.null(offsetName)) 
                                             unadj_frm else 
                                               update(unadj_frm, new=. ~ . + offset(log(get(offsetName)))), 
                                           ), 
                                    data=data_dat[idx_vec, -prd.sng_vec], 
                                    # control=list(maxit=2e2, trace=FALSE), 
                                    Hess=TRUE#, model=FALSE, method="logistic"  # Specifying more options than necessary leads to incompatibilities with the 'multcomp' library
                                    ), 
                               error=function(err) NULL
                               # cat(sub("(^[[:alpha:]]{1})(.+$)", 
                               #         replacement="\\U\\1\\L\\2", 
                               #         x=adj, 
                               #         perl=TRUE
                               #         ), 
                               #     "Proportional Odds regression model did not converge for", 
                               #     endpointName, "in the", analysisLabel, 
                               #     "analysis of the", subgroupLabel, "bootstrap sample\n\n", 
                               #     sep=" "
                               #     )  # Print progress message
                               )  # Fit POM to bootstrap sample
      
      prob_adv_prd_boot_vec <- sapply(sort(unique(data_dat[idx_vec, "treat"])), 
                                   FUN=function(trt) 
                                     if(class(pom_boot_fit) == "polr" && pom_boot_fit$convergence == 0)
                                       {
                                       prob_adv_prd_boot_arr <- predict(pom_boot_fit, 
                                                                        newdata=data.frame(treat=trt, 
                                                                                           data_dat[idx_vec, 
                                                                                                    -match(c("treat", names(prd.sng_vec)), 
                                                                                                           table=names(data_dat)
                                                                                                           )]
                                                                                           ), 
                                                                        type="probs")[, adverseEndpoint_vec, drop=FALSE]  # Derive 2d-array by Draw, Outcome of bootstrap POM predictions of adverse outcome
                                       
                                       return(setNames(mean(rowSums(prob_adv_prd_boot_arr)), nm=trt)
                                              )  # Return as output named vector of bootstrap POM predicted probabilities of adverse outcome
                                       } else 
                                         return(setNames(NA, nm=trt))  # Return as output named NA
                                   )  # Derive vector of predicted bootstrap POM probabilities of adverse outcome
      
      
      rd_boot <- prob_adv_prd_boot_vec[["0"]] - prob_adv_prd_boot_vec[["1"]]  # Derive risk difference (RD) bootstrap estimate
      
      rr_boot <- prob_adv_prd_boot_vec[["1"]] / prob_adv_prd_boot_vec[["0"]]  # Derive relative risk (RR) bootstrap estimate
      
      rrr_boot <- 1e2 * (1 - rr_boot)  # Derive RR reduction (RRR, %) bootstrap estimate
      
      nnt_boot <- ceiling(1 / rd_boot)  # Derive ceiling bootstrap estimate of Number Needed to Treat (NNT)
      
      
      return(c(RD=rd_boot, RR=rr_boot, RRR=rrr_boot, NNT=nnt_boot))  # Return as output named vector of bootstrap risk estimates
      }
    
    
    risk_boot_ls <- sapply(names(pom_fit_ls), 
                           FUN=function(nm){
                             risk_boot_out <- boot(design_dat, 
                                                   statistic=function(data_dat, idx_vec) 
                                                     risk_boot_fn(data_dat, idx_vec, adj=nm), 
                                                   R=boot_size, sim="ordinary", stype="i", 
                                                   strata=design_dat$treat
                                                   )  # Compute stratified bootstrap samples from POM
                             
                             risk_boot_ls <- list(Estimate=risk_boot_out$t0,   # Set bootstrap point estimates for risk parameters
                                                  "95% confidence interval"=sapply(seq_along(risk_boot_out$t0), 
                                                                                   FUN=function(idx_vec) 
                                                                                     boot.ci(risk_boot_out, 
                                                                                             conf=.95, type="perc", index=idx_vec
                                                                                             )$percent[4:5]
                                                                                   ),   # Derive 95% bootstrap confidence intervals for risk parameters
                                                  "P-value"=mean(abs(risk_boot_out$t[, 1] - mean(risk_boot_out$t[, 1], na.rm=TRUE)) > 
                                                                   abs(risk_boot_out$t0[1]), 
                                                                 na.rm=TRUE
                                                                 )  # Compute RD p-value
                                                  )
                             
                             dimnames(risk_boot_ls$"95% confidence interval") <- list(Bounds=paste0(c(2.5, 97.5), " %"), 
                                                                                      Parameter=names(risk_boot_out$t0)
                                                                                      )  # Set dimension names for 2d-array of 95% confidence interval for risk parameters
                             
                             
                             return(risk_boot_ls)  # Return as output list of risk parameters bootstrap inferences
                             }, 
                           simplify=FALSE
                           )  # List by regression adjustment lists of risk parameters bootstrap inferences
    }
  
  
  #######################
  ## Output tabulation ##
  #######################
  
  mod_out_ls <- sapply(names(pom_fit_ls), 
                       FUN=function(nm) 
                         tryCatch(cbind(
                           data.frame(Analysis=analysisLabel, 
                                      Subgroup=unname(subgroupLabel), 
                                      Model="pom", 
                                      Type=nm, 
                                      Endpoint=endpointName, 
                                      Predictors=rep(paste(names(pom_fit_ls[[nm]]$coefficients), 
                                                           paste0("(estimate = ", 
                                                                  round(exp(pom_fit_ls[[nm]]$coefficients), 
                                                                        digits=2
                                                                        ), 
                                                                  ", p-value = NA)"
                                                                  ), 
                                                           collapse=";\n"
                                                           ), 
                                                     times=length(grep("treat", 
                                                                       x=names(pom_fit_ls[[nm]]$coefficients), 
                                                                       value=TRUE
                                                                       )
                                                                  ) + 
                                                       ifelse(is.null(interactionName), 0, 1)
                                                     ), 
                                      Coefficient=c(grep("treat", 
                                                         x=names(pom_fit_ls[[nm]]$coefficients), 
                                                         value=TRUE
                                                         ), 
                                                    if(! is.null(interactionName)) 
                                                      paste(grep("treat", 
                                                                 x=names(pom_fit_ls[[nm]]$coefficients), 
                                                                 value=TRUE
                                                                 ), 
                                                            collapse=" + "
                                                            )
                                                    ), 
                                      "Odds ratio"=exp(c(pom_fit_ls[[nm]]$coefficients[grep("treat", 
                                                                                            x=names(c(pom_fit_ls[[nm]]$coefficients, 
                                                                                                      pom_fit_ls[[nm]]$zeta)
                                                                                                    ), 
                                                                                            value=TRUE
                                                                                            )], 
                                                         if(! is.null(interactionName)) 
                                                           coef(glht(pom_fit_ls[[nm]], 
                                                                     linfct=t(ifelse(grepl("treat", 
                                                                                           x=names(pom_fit_ls[[nm]]$coefficients)
                                                                                           ), 
                                                                                     1, 0
                                                                                     )
                                                                              ), 
                                                                     vcov.=vcovCL(pom_fit_ls[[nm]], 
                                                                                  cluster= ~ subclass
                                                                                  )[names(pom_fit_ls[[nm]]$coefficients), 
                                                                                    names(pom_fit_ls[[nm]]$coefficients)]
                                                                     )
                                                                )
                                                         )
                                                       ), 
                                      rbind(exp(rbind(confint(pom_fit_ls[[nm]], 
                                                              parm=grep("treat", 
                                                                        x=names(pom_fit_ls[[nm]]$coefficients), 
                                                                        value=TRUE
                                                                        ), 
                                                              level=.95, 
                                                              vcov.=vcovCL, cluster= ~ subclass
                                                              ), 
                                                      if(! is.null(interactionName)) 
                                                        confint(glht(pom_fit_ls[[nm]], 
                                                                     linfct=t(ifelse(grepl("treat", 
                                                                                           x=names(pom_fit_ls[[nm]]$coefficients)
                                                                                           ), 
                                                                                     1, 0
                                                                                     )
                                                                              ), 
                                                                     vcov.=vcovCL(pom_fit_ls[[nm]], 
                                                                                  cluster= ~ subclass
                                                                                  )[names(pom_fit_ls[[nm]]$coefficients), 
                                                                                    names(pom_fit_ls[[nm]]$coefficients)]
                                                                     ), 
                                                                level=.95
                                                                )$confint[, c("lwr", "upr")]
                                                      )
                                                )
                                            ),   # N.B.: row-binding required for building single-row table
                                      "P-value"=NA,   # N.B.: p-values derived from foreign packages are based on inappropriate Wald-type tests and misleading asymptotic theory
                                      check.names=FALSE
                                      ), 
                           diag_mod_ls[[nm]][rep.int(1, times=length(grep("treat", 
                                                                          x=names(pom_fit_ls[[nm]]$coefficients), 
                                                                          value=TRUE
                                                                          )
                                                                     ) + 
                                                       ifelse(is.null(interactionName), 0, 1)
                                                     ), ]
                           ),   # Set data-frame of key inferences from Poisson regression models
                           error=function(err) 
                             cbind(
                               data.frame(Analysis=analysisLabel, 
                                          Subgroup=unname(subgroupLabel), 
                                          Model="pom", 
                                          Type=nm, 
                                          Endpoint=endpointName, 
                                          Predictors=paste(predictorNames_vec, 
                                                           collapse=";\n"
                                                           ), 
                                          Coefficient="treat", 
                                          "Odds ratio"=NA, 
                                          "2.5 %"=NA, 
                                          "97.5 %"=NA, 
                                          "P-value"=NA, 
                                          check.names=FALSE
                                          ), 
                               diag_mod_ls[[nm]]
                               )
                           ), 
                       simplify=FALSE
                       )  # Derive list by adjustment type of POM parameter inferences
  
  if(is.numeric(boot_size)) 
    mod_out_ls <- sapply(names(pom_fit_ls), 
                         FUN=function(nm) 
                           tryCatch(data.frame(
                             mod_out_ls[[nm]], 
                             RD=risk_boot_ls[[nm]]$Estimate["RD"], 
                             t(risk_boot_ls[[nm]]$"95% confidence interval"[, "RD"]), 
                             "P-value"=risk_boot_ls[[nm]]$"P-value", 
                             RR=risk_boot_ls[[nm]]$Estimate["RR"], 
                             t(risk_boot_ls[[nm]]$"95% confidence interval"[, "RR"]), 
                             RRR=risk_boot_ls[[nm]]$Estimate["RRR"], 
                             t(risk_boot_ls[[nm]]$"95% confidence interval"[, "RRR"]), 
                             NNT=risk_boot_ls[[nm]]$Estimate["NNT"], 
                             t(risk_boot_ls[[nm]]$"95% confidence interval"[, "NNT"]), 
                             check.names=FALSE
                             ), 
                             error=function(err) 
                               data.frame(
                                 mod_out_ls[[nm]], 
                                 RD=NA, 
                                 "2.5 %"=NA, 
                                 "97.5 %"=NA, 
                                 "P-value"=NA, 
                                 RR=NA, 
                                 "2.5 %"=NA, 
                                 "97.5 %"=NA, 
                                 RRR=NA, 
                                 "2.5 %"=NA, 
                                 "97.5 %"=NA, 
                                 NNT=NA, 
                                 "2.5 %"=NA, 
                                 "97.5 %"=NA, 
                                 check.names=FALSE
                                 )  # Augment data-frame of POM parameters inferences with bootstrap inferences for risk parameters
                             ), 
                         simplify=FALSE
                         )  # Derive list by adjustment type of inferences for POM and risk parameters
  
  
  mod_out_dat <- do.call(rbind, args=mod_out_ls)  # Merge inferences from POMs into single data-frame
  
  
  #####################################
  # Export fitted regression models  ##
  # by analysis data-set and outcome ##
  #####################################
  
  write.table(mod_out_dat, 
              file=file.path(outputPath, 
                             paste0(paste(sub("(^.+)(_class$)", 
                                              replacement="\\1", 
                                              x=endpointName
                                              ), 
                                          "pom", subgroupLabel, 
                                          ifelse(is.character(interactionName), 
                                                 paste("int", 
                                                       paste(interactionLabel, collapse="."), 
                                                       sep="."
                                                       ), 
                                                 "noint"
                                                 ), 
                                          sep="_"
                                          ), 
                                    "_", analysisLabel, ".csv"
                                    )
                             ), 
              quote=TRUE, sep=",", row.names=FALSE, col.names=TRUE
              )  # Export  in .csv format inferences data-frame by outcome, model and subgroup
  
  return(list(endpointName=endpointName, 
              predictorNames=predictorNames_vec, 
              offsetName=offsetName, 
              adverseEndpoint=levels(design_dat[[endpointName]])[adverseEndpoint_vec], 
              analysisLabel=analysisLabel, 
              subgroupLabel=subgroupLabel, 
              interactionNames=interactionName
              )
         )  # Return as output a list of models' terms
  }