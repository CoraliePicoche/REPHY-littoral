This folder contains all data, scripts, graphs and notes on the plankton times series along the French coastline:
- data/raw contains raw files extracted directly from Quadrige (see notes, in French)
- data/ contains a self-made file for all locations (lieu corres), the intermediate files ([0-9]* flortot only.csv and [A-Z]* flortot only.csv correspond to the first classification, made with the same scripts used in Arcachon, which means that groups are the same as those given by D. Maurer. Files beginning with a number correspond to a single point which can change location and therefore bear a different id. They are gathered in files beginning with a letter. Files beginning by subclass or corres correspond to the sites we have decided to work on, with different ways of groupins plankton species (in different groups of genera, or in subclasses - pennate or centric)
- script/ contains all scripts used on data.
- for now, graphe/ mostly contains time series. Files ending with "var time" show the variables that were collected, and how long, for each site (which can help see if nutrients were collected for a certain site, for example). Files ending with "per species" show the time series and the handling of missing values we used in Arcachon (not satisfactory). Subfolders contain different time series corresponding to different groups of genera (according to Tania Hernandez, to the ones we used in Arcachon, or to Subclass).
- Rapport/ is pretty self-explanatory for now

There are miscellaneous files in the folder:
-corres.csv are only files with codes and groups together
-db_taxonomy_REPHY and its corresponding file taxonomy_pencen.RData contain all matches between REPHY and WORMS

