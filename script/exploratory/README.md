This folder contains exploratory analyses which were not used in the final paper. 

* `abundance_proportion.r` compares the proportion of diatoms and dinoflagellates, and computes the proportion of studied genera against the total population
* `comparaison_REPHY_SRN.r` compares the data from the REPHY monitoring and the SRN monitoring
* `compare_evt_interaction.r` shows the mean effect of the environment (using both temperature and salinity at the same time) against the mean strength of inter- and intra-group interaction strength in each site
* `compare_quadratic.r` computes AIC and BIC on MAR models using a quadratic effect of the environmental variables to try and characterize a niche effect. Quadratic effects can be viewed thanks to `view_niche_quadratic.r` 
* `comparison_matrix_B.r` produces a figure with different metrics (weighted connectance, linkage density, mean values) when considering only the negative or only the positive values of the community marix B
* `compute_common_taxa.r` computes the mean number of taxa each site has with each other
* `correlation_covariates.r` computes the cross- and auto-correlation in the two environmental variables (salinity and temperature) at each site
* `distrib_interaction_coefficients.r` checks the distribution of the interaction coefficients and test their gaussian distribution
* `exploration_repetition.r` shows the variation in abundances of the species in all sites of one region (basically, to show that they are synchronous)
* `Gantt_REPHY.r` presents the number of points for each variable and each site of the REPHY monitoring (i.e., not only the ones we kept)
* `get_taxonomy_par_site.r` presents the richness and proportion at each site for genus, family, class, family
* `interaction_moy_var_per_module.r` compares mean and variane of inter and intra for each plankton cluster (pennate, centric, dino, others) 
* `nutrients.r` computes the concentration of nutrients and plots them for sites in which the data is available
* `search_for_FLORPAR.r` checks partial flora instead of total flora for phytoplankton. 
* `test_coeff_vs_seasonality.r` compares the strength of the competition to the strength of seasonality, based on an idea by Usinowicz 2017 in Nature
* `times_series_reconstruite.r` compares the random and spectrum-based reconstructions of the time series
