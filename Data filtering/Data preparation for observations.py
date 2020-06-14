#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np


obs = pd.read_csv('individual_phenometrics_data.csv')
site = pd.read_csv('ancillary_site_data.csv')

# Screen out the required column names
obscols = ['ObservedBy_Person_ID', 'Partner_Group', 'Site_ID','Latitude',
          'Longitude','Individual_ID','Multiple_Observers']
sitecols = ['Site_ID', 'Site_Type']
obs=obs[obscols]
site=site[sitecols]

# Merge the table
obssite = pd.merge(obs, site, on='Site_ID',how='inner')

# Filter one observation per Individual_ID
ois = obssite.groupby(['Individual_ID','Site_ID']).head(1)
ois.to_csv('one_individual_site.csv',index=False)





