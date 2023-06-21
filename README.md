## NHSEI Methods: Proportional Odds Models

### About this method

The method consists in applying a Proportional Odds Model (POM) ordinal regression strategy to a comparative (intervention vs control) study set-up, with the aim of estimating not only an intervention effect, but also an array of diagnostic measures that are commonplace in randomised controlled trial analyses.

In the provided illustrative code (script ‘pom_phm.r’), a POM is fitted to a fictitious rectangular data-set to represent the relationship between on one hand a set of predictors, including a binary intervention assignment variable, and on the other a count variable coarsened into an ordered categorical variable (factor).  The code produces as output a .csv spreadsheet subsuming inferences (consisting of point and interval estimates) of the intervention effect, of the Absolute Risk Difference, of the Relative Risk, of the Relative Risk Reduction and of the Number Needed to Treat, as well as an array of and model diagnostic statistics. The script ‘pom_code.r’), which encodes the model fitting engine called upon by the ‘pom_phm.r’ script, is not supposed to be altered by the user unless model capability expansion is desired.

The proposed code should be considered for use by an analyst with statistical background in regression in the presence of data from a comparative study contrasting cases (e.g. patients, NHS trusts, ICSs…), on which a variety of characteristics have been recorded, assigned to an intervention or a control group.  The aim of the study is to ascertain the effect of the intervention of interest when this is measured via either a count metric (e.g. number of emergency admissions into hospital) or an ordinal factor (e.g. a coarsened version of a count metric). 

For more information about the method, including how is being used, when and by who, please refer to the Methods toolbox documentation on our [FutureNHS workspace](https://future.nhs.uk/DataMeth/grouphome).

### Requirements

The method requires sourcing in the R platform the ‘pom_code.r’ script, and following the steps outlined in the illustrative ‘pom_phm.r’ example script.  File and folder paths in the ‘pom_code.r’ script will need editing by the analyst to reflect their file system; additionally, the R libraries ‘boot’, ‘MASS’, ‘multcomp’, ‘nnet’ and ‘sandwich’ are required dependencies.  Additional libraries loaded in the preamble of the ‘pom_phm.r’ script (i.e. ‘MatchIt’ and ‘optmatch’) are only needed to carry out the analysis outlined in the script as an illustration.  In particular it should be noted that the data input to the code should include a binary intervention assignment variable named ’treat’.

### Summary of the code

* [pom_phm.r]  This script proceeds to derive a matched sample (‘matched.dat’) from applying the Optimal Matching algorithm to a freely available data-set (‘original.dat’) including a binary intervention assignment variable (‘treat’).  It subsequently produces from the matched data-set a plot of selected categorical predictor levels vs cumulative log-odds of the response variable, which is useful to visually check the validity of the proportional odds assumption underpinning a POM.  It subsequently proceeds to fitting the POM to the data-set via the core ‘fitPOM_fn()’ function, which is coded in the ‘pom_code.r’ script.
  
* [pom_code.r]  This script defines a single function (namely ‘fitPOM_fn()’) by checking the congruency of the arguments fed by the analyst as inputs.  It subsequently defines the regression models as described by the analyst among the function’s arguments, and fits a POM and a Log-Linear Model (LLM) in both their adjusted (to specified predictors in the data-set) and unadjusted variants.  A data-frame (in essence a flat table) comprising an array of diagnostic statistics (including residuals, deviance, Chi-Square goodness-of-fit test, LMM vs POM comparison, AIC…) is then built and merged with another data-frame comprising statistical inferences (consisting of point estimates and 95% confidence intervals) around the intervention effect, as well as bootstrap estimates of risk parameters (namely the Absolute Risk Difference, the Relative Risk, the Relative Risk Reduction and the Number Needed to Treat).  Lastly the code dumps the combined data-frame onto the hard drive at a user-specified location.


### Data sources

The convenience ‘lalonde’ data-set bundled with the ‘MatchIt’ R library.

### Authors

* Dr Stefano Conti, Senior Statistician - stefano.conti@nhs.net

