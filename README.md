# PCORI Project Code

This is a Patient Centered Outcomes Research Institute (PCORI)-sponsored project to improve methods for analyzing comparative effectiveness studies

## Antiretroviral Therapy in HIV Patients and Cardiovascular Disease

Advances in antiretroviral therapy (ART) have dramatically reduced mortality from HIV, enabling reclassification of HIV as a chronic condition. Numerous studies suggest that some drugs increase the risk of cardiovascular disease, although findings are inconsistent. Studies differ largely due to methodological choices, including study design, definition of exposure, and approaches to handling missing data. It is crucial to incorporate information on drug exposure and other confounders over time; patients vary their regimens over time and for reasons that may be related to their condition. Thus, in order to not implicate the wrong drug, information on changes in regimen and other factors must be considered.

## Studies with Time-Varying Confounders Can Induce a Missing Data Problem

Including information over time complicates the analysis, however. One such complication is the introduction of missing data. Common methods for handling missing data yield misleading descriptions of relationships. Appropriate methods for handling missing data are computationally burdensome; software does not exist for many situations and the analyst must rely on his/her own programming skills to implement specialized techniques. Multiple imputation (MI) is a reasonably accessible and theoretically sound method for handling missing data. Available in mainstream software, its special application is required due to the unique issues posed by time-varying covariates and outcomes that are only partially observed for those individuals who do not experience a cardiovascular event during the observation period. We propose an extensive simulation study to evaluate commonly applied methods to this setting, to investigate the performance of standard MI in this context, and to adapt and evaluate MI methods utilized in a longitudinal setting where the outcome is fully observed to this particular setting. Based on our findings we will develop concrete guidelines on how to use MI in the context of partially observed outcomes and time-varying covariates. We will develop user-friendly open-source software in order to optimize the use of recommended methods and to eliminate lack of software as a barrier to employing missing data methods. Finally, we will illustrate methods considered on data from the US veteran population of HIV-infected individuals using the Veterans Health Administration's rich longitudinal Clinical Case Registry (CCR), the analysis of which motivated this proposal. This work has the potential to greatly impact patients living with HIV. Currently, there is no consensus on which ART agents increase cardiovascular risk.

## Goals

To address the implications of methodological choices for handling missing data when conducting comparative effectiveness research in the longitudinal setting.

## Aims

### Aim 1: Develop methods to simulate correlated time-varying covariates and right-censored outcomes

### Aim 2: Characterize variability in common missing data methods applied in this context

### Aim 3: Develop new methods for handling missing data for these studies

### Aim 4: Provide open-source software for implementing recommended methods

## About this Code

This code is a work in progress and is provided as-is. More comprehensive documentation should be forthcoming. All code is written in R and parallelization files are optimized for used on a SLURM-based cluster computing system.

Questions should be directed to [Kristopher Kapphahn](https://med.stanford.edu/profiles/kristopher-kapphahn).


## About the Subfolders

-<b>GENCOV</b>: R Code to generate correlated covariates (binary drug indicators, continuous proportion of time on drug variables, normal variables, and categorical variables), with or without time effects within subjects

Covariate generation scripts are split into several separate code files that need to be run in sequence (see ```simulate.R``` for the workflow). The first, ```init_variables.R```, assigns global variables that control simulation parameters such as covariate means, covariate standard deviations, and correlation structures. It loads several matrices of parameters that can be manually edited as csv files by the user: ```parameters_time_vary.csv``` and ```pcor_time_vary.csv``` (a population correlation matrix). Files ```load_functions.R``` and ```jointly_generate_binary_normal_modified_v2.R``` define necessary functions, including a modified version of a function in package ```BinNor```. Finally, ```simulate.R``` is a parallelized wrapper script (requiring modification for non-SLURM cluster computing systems) allowing user-inputted arguments (representing the sample size, number of observations per subject, and number of simulations), computes statistical performance measures over multiple simulations, and creates performance plots.

The file ```corrBound.R``` defines accessory functions that compute lower and upper bounds on correlations between variables. While these functions are not called internally in the simulation process, they may be a useful pre-processing tool to guide choices of variable parameters.

The file ```generateSbatch_mm.r``` is for use with a SLURM-based cluster computing system and automatically generates sbatch files to run the simulation code. A few examples are in the directory ```sbatch_files```. File ```push_to_sherlock.R``` is used to push files to the cluster. 

-<b>GENSURV</b>: Code to generate survival times from an extended Cox model

Program "generateByScenario" calls the function mysim which is located on filemysimRangeScenarioSB.R. This file contains all the code necessary to run a simulation where various parameters are tested. This include different assumptions on the covariates used, several censoring methods, among others. Once data is generated, models are fitted and statistics computed using program "fitSurvivalModelsSHER_COM.R".

A program that looks at run times depending on the various parameters tested can be found here.

-<b>IMPMISS</b>: Code to generate auxiliary variables and use them to impose missingness. 
