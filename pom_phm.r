##############
## Preamble ##
##############

rm(list=ls())  # Clear work-space


library(MatchIt)  # Load library enabling Optimal Matching routines

library(optmatch)  # Load library with core Optimal Matching engine


old.par <- par(no.readonly=TRUE)  # Set default graphical parameters


project.dir <- "C:/Users/StefanoConti/OneDrive - NHS England/Documents/Work/NSHE"  # Set project folder


source(file.path(project.dir, 
                 "R/DA/Projects/Ongoing/Methods tool-box/pom_code.r")
       )  # Source R script with POM regression code


##########################
## Build counterfactual ##
##########################

original.dat <- lalonde  # Set original (unmatched) data-frame: 
                         # a sub-sample from the intervention group 
                         # in the National Supported Work Demonstration 
                         # (NSW) and the comparison sample from the 
                         # Population Survey of Income Dynamics (PSID).



original.dat <- within(original.dat, expr=
                         re78_class <- cut(re78, 
                                           breaks=c(min(re78), 25e2 * seq.int(4), max(re78)), 
                                           labels=c("- $2,499", "$2,500 - $4,999", "$5,000 - $7,499", "$7,500 - $9,999", "$10,000 -"), 
                                           include.lowest=TRUE, right=FALSE, ordered_result=TRUE
                                           )  # Derive coarsened version of 're78' as ordered factor
                       )

ps.frm <- treat ~ age + educ + race + nodegree + married + re74 + re75  # Set propensity score model for general matching example


om.out <- matchit(ps.frm, data=original.dat, method="optimal", 
                  distance="glm", link="logit", estimand="ATT", 
                  ratio=1
                  )  # Derive Optimal Matching output for general matching example


matched.dat <- match.data(om.out, 
                          data=original.dat, distance="prop.score"
                          )  # Derive matched data-frame for general matching example


####################
## Validating the ##
## PO assumption  ##
####################

lno.ls <- with(matched.dat, expr=
                 sapply(c("treat", "race", "married", "nodegree"), 
                        FUN=function(prd) 
                          tapply(unclass(re78_class) - 1, 
                                 INDEX=list(get(prd)), 
                                 FUN=function(y) sapply(seq.int(nlevels(re78_class) - 1) - 1, 
                                                        FUN=function(z) qlogis(mean(y <= z))
                                                        )
                                 )
                        ), 
               simplify=FALSE
               )  # Derive list by predictors levels of log-odds of cumulative outcome


lno.arr <- do.call(cbind, 
                   args=unlist(lno.ls, recursive=FALSE)
                   )  # Format list of cumulative outcome log-odds as 2d-array by Response, Predictor

dimnames(lno.arr) <- list(Response=sub("(^.+)(- .+$)", 
                                       replacement="\\2", 
                                       x=levels(matched.dat$re78_class)[-nlevels(matched.dat$re78_class)]
                                       ), 
                          Predictor=sub("\\.", replacement=" = ", x=colnames(lno.arr))
                          )  # Rename margins of cumulative outcome log-odds 2d-array


par(mar=c(5, 7, 4, 2) + .1)  # Move left plot margin inwards

plot(lno.arr, col(lno.arr), 
     type="p", axes=FALSE, 
     col=seq_along(dimnames(lno.arr)$Response), pch=15, 
     main="Cumulative Log-Odds\nof 1978 Income", 
     xlab="Log-odds", ylab="", 
     xlim=range(pretty(range(lno.arr), n=5))
     )  # Plot log-odds of cumulative outcome

axis(1, at=pretty(range(lno.arr), n=5))  # Overlay x-axis labels onto plot

axis(2, at=seq_along(dimnames(lno.arr)$Predictor), 
     labels=dimnames(lno.arr)$Predictor, las=1
     )  # Overlay x-axis labels onto plot

abline(h=seq_along(dimnames(lno.arr)$Predictor), 
       col=8, lty=2
       )  # Overlay horizontal lines across variables for readability

legend("topleft", 
       legend=dimnames(lno.arr)$Response, 
       fill=seq_along(dimnames(lno.arr)$Response), 
       bty="o"
       )  # Overlay legend with "response" keys

dev.print(png, 
          file=file.path(project.dir, 
                         "R/DA/Projects/Methods tool-box/lno.png"), 
          width=1.5 * 480, height=480
          )  # Export in .png format log-odds of cumulative response plot

# unlink(file.path(project.dir, 
#                  "R/DA/Projects/Methods tool-box/lno.png")
#        )  # Remove .png file of log-odds of cumulative response plot

par(old.par)  # Restore default graphical parameters


########################
## Outcome regression ##
########################

fitPOM_fn("re78_class", 
          predictorNames_vec=c("age", "educ", "race", "nodegree", "married", "re74", "re75"), 
          offsetName=NULL, 
          design_dat=matched.dat, 
          adverseEndpoint_vec <- seq.int(nlevels(matched.dat$re78_class))[-nlevels(matched.dat$re78_class)], 
          analysisLabel="full", 
          subgroupLabel="max", 
          outputPath=file.path(project.dir, "R/DA/Projects/Methods tool-box"), 
          interactionLabel=NULL, 
          interactionName=NULL, 
          boot_size=5e3
          )  # Fit POM regression to response from general matched data-frame