Please cite as:
Dorigatti et. al., SARS-CoV-2 antibody dynamics, within-household transmission and the impact of contact tracing from community-wide serological testing in the Italian municipality of Vo’. Available as a pre-print (not peer-reviewed) at SSRN: http://dx.doi.org/10.2139/ssrn.3792892. 

# Repository structure
This repository contains: 

Data 
- the data collected in Vo' in two sequential serological and nasopharyngeal swab surveys conducted in May and November 2020; for completeness, the dataset also includes the data collected in the surveys conducted in February and March 2020 (see Lavezzo et al., Suppression of a SARS-CoV-2 outbreak in the Italian municipality of Vo', 2020, Nature);
- the contact tracing data collected independently from this study by the local public health authorities in Vo' at the start of the epidemic; 
	
Code 
- to estimate SARS-CoV-2 seroprevalence; 
- to assess the association between antibody response and demographics, health indicators and disease severity; 
- to estimate SARS-CoV-2 attack rates by household size;  
- to fit a within-household SARS-CoV-2 transmission model to the observed data; 
- to estimate the sensitivity of contact tracing; 
- to estimate the impact of contact tracing on the epidemic dynamics and final size; 
- to assess the association between SARS-CoV-2 infection and comorbidities/medication history. 

## Data
The data folder contains the file `Vo_serology_data.xlsx`. This is a line list with the results of the serological assays performed in May and November 2020, along with the contact data, dates and results of swab testing, symptoms and hospitalisation data collected in two previous surveys conducted in the same population in February and  March 2020 (Lavezzo et al., 2020, Nature). The file `Vo_serology_data.csv` contains the same data included in spreadsheet 2 of `Vo_serology_data.xlsx`. 

The file `Contact_tracing_Vo.xlsx` contains the reported named contacts and dates of isolation for the traced contacts. 

The folder 'household_model' contains the files with the default model parameters, model variants and MCMC settings for the original model developed by Fraser et al., 2011, AJE, and the extended and 2-groups models. 

## Scripts
The scripts folder contains the main scripts for each separate analysis
- `Descriptive_analysis_script_1.R` is the script producing the estimates reported in Figure 1 and the attack rates by household size used for the within-household transmission model; 
- `Estimate_true_seroprevalence_script_2.R` is the main script for estimating the true eroprevalence and the assay performances using a multinomial likelihood calibrated to the data using Markov Chain Monte Carlo (MCMC);
- `Seroprevalence_script_3.R` is the main script to estimate and compare the seroprevalence using the ground truth definitions and the estimated assay-specific performances; 
- `Association_antibody_titres_script_4.R` is the main script to run the association analysis between the antibody titres and demographic, health and disease severity indicators; 
- `Association_antibody_decay_script_5.R` is the main script to run the association analysis between the antibody decay rates and demographic, health and disease severity indicators; 
- `Fit_original_household_model_script_6.R` is the main script to fit the within-household SARS-CoV-2 transmission model developed by Fraser et al., 2011, AJE to the within-household attack rates observed in Vo'. We use MCMC for parameter inference. This script requires script 1 to be run beforehand. 
- `Fit_extended_household_model_script_7.R` is the main script to fit an extension of the model developed by Fraser et al., 2011, AJE to the within-household attack rates observed in Vo'. We use MCMC for parameter inference. This script requires script 1 to be run beforehand. 
- `Fit_2_groups_household_model_script_8.R` is the main script to fit a 2-group model to the within-household attack rates observed in Vo'. We use MCMC for parameter inference. This script requires script 1 to be run beforehand. 
- `Plot_DIC_script_9.R` is the main script to compare the performance of the original, extended and 2-groups within-household transmission models. This script requires scripts 6, 7 and 8 to be run beforehand. 
- `Plot_SITP_overdisp_fit_original_model_script_10.R` is the main script to output the model fit, the SITP and the overdispersionin the offspring distribution of the best model. This script requires script 6 to be run beforehand. 
- `Impact_of_contact_tracing_11.R` is the main script to compute the impact of contact tracing on the spread of COVID-19 in Vo' in spring 2020. Under the hood, this script runs a SEIR model that has previously been published (Lavezzo et. al., 2020, Suppression of COVID-19 outbreak in the municipality of Vo, Italy, medRxiv, doi: 10.1101/2020.04.17.20053157). This script is independent of all other scripts.
- `Complementary_SI_info_12.Rmd` computes complementary info for Figures S1, and Tables S5 and S6.
- `Neutralization_decay_script_13.R` is the main script that calculates the halflife of neutralising antibodies observed between May and November, considering only subjects who were positive in May and not increasing between May and November.


## R 
The R folder contains the functions used by the main scripts. 

- `Descriptive_analysis_script_1.R` uses script functions_descriptive_analysis.R
- `Estimate_true_seroprevalence_script_2.R` uses script functions_multinomial_likelihood.R
- `Seroprevalence_script_3.R` uses script functions_seroprevalence.R
- `Association_antibody_titres_script_4.R` uses script functions_association_antibody_titres.R
- `Association_antibody_decay_script_5.R` uses script functions_association_antibody_slopes.R
- `Fit_original_household_model_script_6.R` uses scripts functions_fit_model.R and functions_plot_fitted_model.R
- `Fit_extended_household_model_script_7.R` uses scripts functions_fit_extended_model.R and functions_plot_fitted_model.R
- `Fit_2_groups_household_model_script_8.R` uses scripts functions_fit_2_groups_model.R and functions_plot_fitted_model.R
- `Plot_DIC_script_9.R` uses scripts functions_plot_deviance_DIC.R and functions_plot_fitted_model.R
- `Plot_SITP_overdisp_fit_original_model_script_10.R` uses scripts functions_fit_model.R and functions_plot_fitted_model.R
- `Impact_of_contact_tracing_11.R` depends on functions defined in the files functions_model_11.R, functions_analysis_11.R, functions_figures_11.R

## Contact tracing analysis 
The contact_tracing_analysis folder contains the data and Python code used to estimate the performance of contact tracing. Instructions on how to run the code are given in the file Command_Line_Statistic_forge. 

