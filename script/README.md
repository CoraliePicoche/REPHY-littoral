This folder contains all the scripts that were used for the analyses in the manuscript "Strong self-regulation and widespread facilitative interactions in phytoplankton communities" (Picoche & Barraquand 2020). They correspond to results in the main text, the supporting information and the responses to the reviewers which were not necessarily shown in the final paper. Some of the scripts are in the exploratory folder, corresponding to previous analyses which were not used in the final paper. 

The main analyses are as follows. 

#### Main text
*  `check_continuity.r` computes the length of the gaps in the datasets of the different REPHY sites in order to choose the time series which have the less gaps (detailed in the Methods of the main text).
*  `exploitation_FLORTOT.r` turns the raw datafile (directly extracted from the database of IFREMER, Quadrige) and turns into easier-to-use .txt files
*  `hydro_times_series_par_classif.r` does the same for the hydrological variable  
*  `MARSS_clean.r` is a function that automatizes and standardizes the MARSS analysis we perform on the datasets 
*  `MAR_single.r and MAR_single_for_Arcachon.r` use the function `MARSS_clean.r` to actually analyse the datasets. Arcachon is treated separately because files were already treated in a slightly different way in a previous paper (Barraquand et al. 2018, Oikos)
*  `matrix_MAR_clean.r` is a function that takes the MARSS object into an easier-to-treat matrix 
*  `complexity_stability_MainFig.r` compares the resilience (measured by the dominant eigenvalue) to different network metrics and produces Fig. 2 in the main text
*  `biotic_interaction_matrices_MainFig_v4.r` produces Fig. 1 in the main text and Fig. S8 in the supporting information (respectively pennate-centric and unconstrained interaction scenario), that is the interaction strengths for all sites in different regions, and shows the percentage of positive interactions among them.  
*  `generality_vulnerability_MainFig.r` computes the impact and vulnerability of each species, compares them to their self-regulation and produces Fig. 3 in the main text 

#### Supporting Information
*  `compute_BIC.r` computes the BIC of the MAR models associated to different interaction scenarios based on the phylogeny of all groups of species. It produces Fig. S3 in the suporting information
*  `map_geom.py` produces a map of the different sites used in the analyses (Fig. S1)
*  `abiotic_effects.r` produces a figure showing the effects of temperature and salinity on the different groups of species in the MAR models (Fig. S4)  
*  `interaction_moy_var_SIFig.r` compares the mean and variance of the intra- and intertaxa interaction strengths (Fig. S6) 
*  `new_hessian.r` computes the Fisher Information matrix of the interaction parameters to look for artefacts due to statistical constraints (discussed in section *MAR(1) models* of the supporting information)
*  `trying_mutual_commensalist.r` compare the types of interactions (mutualist, commensalist, competitive) obtained in the interaction matrices and produces Fig. S5.
*  `quantif_unconstrained.r` computes number of indicators on the unconstrained matrices and compares them to the more constrained pennate-centric interaction matrices used in the main paper. It produces the results for section *Comparison with a full interaction matrix* in the supporting information
*  `regulation_vs_abundance.r` shows the self-regulation/intraspecific interaction strength against the abundances of each groups of species, producing Fig. S7 in the supporting information 
*  `time_series_SI.r` shows the time series of the 5 most abundant groups of species at each site (Fig. S2)
*  `compare_logB.r` compares log(B) and B-I where B is the interaction matrix obtained from MAR analyses. This is detailed in section *Connection to continuous-time models* in the supporting information
*  `ratio_abundances_conversion_MARBH.r` computes the ratios of abundances at equilibrium to compare orders of magnitude of the interaction strengths between MAR and Beverton-Holt models (discussed in the section *Connection to Lotka-Volterra competition dynamics* in the supporting information'


#### Responses to referees only
*  `biotic_interaction_matrices_with_segment.r` produces a figure similar to Fig. 1 in the main text, but replacing the bars by point estimates + standard deviations. This was only used to answer the referees. 
*  `spectre_eigenvalues.r` compares the spectra of the eigenvalues of the interaction matrices. This was only used in the response to the referees.


The folder `MARinteraction_review` contains all data and the script for Fig. 4 in the main text and section *MAR references and analysis* in the supporting information
