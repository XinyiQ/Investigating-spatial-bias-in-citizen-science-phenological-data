# Investigating-spatial-bias-in-citizen-science-phenological-data
Xinyi_Qiu MSc Thesis 

The folder "Data Filtering" contains data for processing voluntary phenological observations and codes for processing.
- Data sheet "ancillary_individual_plant_data.csv", "ancillary_site_data.csv" and "individual_phenometrics_data.csv" are input data.
- Data sheet "one_individual_site.csv" is the output data from processing.
-"Data preparation for observations.py" is the python script for processing data.

The folder "Modelling" contains data and code for the construction of the model.
- In folder "DATA", there are 17 spatial covariates used to construct the model and are in the form of raster images ".tif". "Obs_bound_Proj" is the study area. "Obspoint_Proj" is the observations used to construct the model.
- "Model.R" is the R script for constructing the model.
