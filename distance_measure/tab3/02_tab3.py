import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import pearsonr, spearmanr, ttest_ind
from statsmodels.regression.linear_model import OLS

def bootstrap(df, measure):
   corrs = []
   for _ in range(1000):
      df_smpl = df.sample(n = df.shape[0], replace = True)
      corrs.append(spearmanr(df_smpl["infection"], df_smpl[measure])[0])
   return corrs

def format_coefficient(df, measure, baseline_bootstrap):
   corr = spearmanr(df["infection"], df[measure])
   rsquared = "%1.4f" % corr[0]
   bootstrap_model = bootstrap(df, measure)
   ttest = ttest_ind(baseline_bootstrap, bootstrap_model, equal_var = False)
   if ttest[0] < 0 and ttest[1] < 0.01:
      rsquared = "\\textbf{%s}" % rsquared
   if corr[0] < 0:
      rsquared = "\\textcolor{red}{%s}" % rsquared
   if corr[1] < 0.001:
      rsquared = "%s***" % rsquared
   elif corr[1] < 0.01:
      rsquared = "%s**" % rsquared
   elif corr[1] < 0.1:
      rsquared = "%s*" % rsquared
   return rsquared


netlabel = {"er": "ER", "ba": "BA", "pc": "PC", "lfr": "LFR"}
measures = ("ge", "mmc", "annihil", "emd", "spl single", "spl avg", "spl complete", "gft", "mapp", "count")

print("\\begin{table}")
print("\\centering")
print("\\begin{tabular}{l|rrrrrrrrr}")
print("Topology & %s\\\\" % ' & '.join(measures).upper())
print("\\hline")

for network in netlabel:
   df = pd.read_csv("vm_table_%s.csv" % network, sep = "\t")
   baseline_bootstrap = bootstrap(df, "count")
   line = []
   line.append("%s" % netlabel[network])
   for measure in measures:
      line.append(format_coefficient(df, measure, baseline_bootstrap))
   print("%s\\\\" % (" & ".join(line)))

print("\\end{tabular}")
print("\\caption{Correlation coefficients of each measure with the infection parameter in the cascade model, for all types of random networks. Bold indicates the measures with a correlation coefficient significantly higher than the Count baseline (t-test based on bootstrapping). Red indicates measures negatively correlated with the infection parameter. * $p < 0.1$, ** $p < 0.01$, *** $p < 0.001$.}\n\\label{tab:vm}")
print("\\end{table}")

