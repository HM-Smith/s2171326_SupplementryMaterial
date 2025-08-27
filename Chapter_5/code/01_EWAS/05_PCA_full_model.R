python

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import datatable as dt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import statsmodels.api as sm
import colorsys
import math


# Import data
################################################################

df = dt.fread("methylation_data.txt", sep = " ")
df = df.to_pandas()
df = df.set_index("IID")
df = df.drop("FID", axis = 1)
names = df.index.values.tolist()

cpgs = [i.rstrip() for i in open("bmi.txt", "r").readlines()]


# Filter to CpGs of interest
############################

df1 = df[cpgs]


# Standardize
################################################################

df1 = StandardScaler().fit_transform(df1)


# PCA
################################################################

pca = PCA(n_components=18413)
pcs = pca.fit_transform(df1)
pca_df = pd.DataFrame(data = pcs, columns = ["PC%s" % i for i in range(1,18414)], index = names)
expvar_df = pd.DataFrame(data = pca.explained_variance_ratio_, index = ["PC%s" % i for i in range(1,18414)], columns = ["ExpVar"])


# Export results
################################################################

pca_df.to_csv("bmi_PCA_df.tsv", sep = "\t", index_label = "ID")
expvar_df.to_csv("bmi_PCA_explainedvariance_full_model.tsv", sep = "\t", index_label = "PC")





