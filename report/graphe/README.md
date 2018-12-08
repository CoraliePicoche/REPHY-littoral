This folder contains all graphs intended to explore our data set, or to visualize results from different analyses (MAR estimates, for instance).

The *localisation* folder contains maps locating the different REPHY monitoring points as well as MeteoFrance sites.

The *exploration_raw* folder contains multiple ways of looking at the data:
- all files ending with "var_time.pdf" show all variables which were monitored, for each sampling point. They do not represent the time series themselves, but only the presence of a variable, in order to know which covariate we can use for each site.
- "boxplot_gap_all.pdf" and "gap_per_season.pdf" shows the average duration between two sampling dates, especially when it is over 2 weeks, which is REPHY prescribed frequency. The point was to remove sites that were not monitored regularly enough.
- "Sites_REPHY.pdf" and "Gantt_REPHY.pdf" focus on plankton data: they show the same information as above, but only for the "FLORTOT" variable: you can have access to the length and number of gaps for each point at a glance
- the file "comparaison_REPHY_SRN" was used to compare data between the two programs for our northern points, as Hernandez has focused on these points. However, none of them seems appropriate for a MAR analysis (they are too short and gappy).
- the same was plotted for the "FLORPAR" and "FLORIND" variables: these code correspond to a focus on toxic and/or abundant species. We wanted to check that synchrony analyses would be better if performed on this variable instead of FLORTOT. Does not seem to be the case.
- still focusing on planktonic data, the files ending with "per_species.pdf" show a first attempt (and failure) at using exactly the same groups of species as in our first paper on Arcachon, with the same reconstruction method.
- after such failure, in order to choose the right species, the file "repartition_taxonomy" and "nb_groups" showed the number of groups of species we could have according to the taxonomic level we wanted to keep, and the number of missing values we were ready to accept, especially if we wanted to have homogeneous groups of variables. We then chose to focus on the genus level, with Hernandez et al. (2015) classification. For other classification, see the folder *time_series_classif_old*

The *time_series* folder contains visualization of the data. 
- the *Hernandez* folder only shows raw data, with the proportion of missing value and the average abundance of each group, according to Hernandez et al. (2015) classification. A color code indicates for each point how dominant the species is at one date.
- the *hydro* folder focuses on covariates, and their correlation and cross-correlation
- the *plankton* folder also compares covariation of plankton groups. The "comparaison_SSA_random" files compare a reconstruction method based on SSA (with different SSA method), and a random replacement of missing value as done by Hampton et al. 2013. We finally chose to use Hampton's method, to be consistent with our first article on Arcachon.

The *MAR_estimates* folder contains all visualization of interaction matrices and comparison of AICc and BIC for MAR estimates. They are separated per location (BZ=South Brittany, MO=Marennes Oleron, AR=Arcachon, SU=South of France, the Mediterranen sites), and per groups of species selected. The interaction models are specified in the names of the files (null=diagonal matrix, unconstrained means a full interaction matrix, pencen means that pennate and centric cannot interact with each other, diatdin means that diatoms and dinoflagellates cannot interact together but pennate and centric diatoms can, inter means that there are interactions inter-phyla but no intraphyla: a diatom can interact with a dinoflagellate, but not with another diatom).
There were at first some attempts at building models using NEI either as a covariate or as a species: all results were put in the *with_NEI* folder. There were also models built using only common species from BZ, MO and SU, not taking into account AR: results are in COMMON/without_AR. Both types of models can be ignored now. 
For each location, there is a difference between "subsite_specific" and "site_specific". In the first folder, groups of species can differ between sites of the same location (in Marennes-Oléron, Auger and L'Eperon can have different interacting species), while the "site_specific" folder contains only results for common species between different subsites (in Marennes Oléron, we only consider species that are well-resolved in Auger, L'Eperon and Cornard).
Finally, in the *COMMON* folder, files ending with "subsampling" show interaction coefficients for species that are common to our 4 locations (which means only 4 species), but which were estimated with a model taking into account all species common at one place (which means 14 species at BZ, instead of the resulting 4, for instance). It was just a test while real estimates were ran with only 4 species taken into account.