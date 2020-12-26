# Monte Carlo Simulation Study of EM-algorithm for missing data

Following Caudill, Ayuso and Guillén (2005) we conduct an investigation of the performance of the EM-algorithm for missing data, by the means of a Monte Carlo Simulation study. The example we have in mind, is imperfect detection rates in audits, and its consequences and the ability of the EM-algorithm to address this. 

Our discussion paper is available from: https://www.nhh.no/en/research-centres/nocet/research/discussion-papers/

param_config_table.xlsx shows the 11 configurations we present here. 

scripts contains two files; v2_functions.R which defines our most important functions, which v2_mc_em_nonparam_boot.R sources, and then runs all 11 Monte Carlo Simulations.

The simulation results can be found in individual data files in the folder simulation_results, while an HTML overview can be found in the folder reports_and_data_summariers

Reference: 
Caudill, S. B., Ayuso, M., & Guillén, M. (2005). Fraud detection using a multinomial logit model with missing information. Journal of Risk and Insurance, 72(4), 539-550.

Feel free to contact us with questions and suggestions: 

Jonas Andersson
Jonas.Andersson@nhh.no
https://www.nhh.no/en/employees/faculty/jonas-andersson/

Andreas Olden
Andreas.Olden@nhh.no
https://www.nhh.no/en/former-employees/andreas-olden/
