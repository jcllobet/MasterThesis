import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from collections import defaultdict

df = pd.read_csv("world_bank_gdp_pc.csv").set_index("Country Code")

country_decade_growth = []
for country in df.index:
   country_row = df.loc[country]
   country_decade_growth.append((country, 1, np.abs(np.log(country_row["1990"]) - np.log(country_row["1980"]))))
   country_decade_growth.append((country, 2, np.abs(np.log(country_row["2000"]) - np.log(country_row["1990"]))))
   country_decade_growth.append((country, 3, np.abs(np.log(country_row["2010"]) - np.log(country_row["2000"]))))
   country_decade_growth.append((country, 4, np.abs(np.log(country_row["2018"]) - np.log(country_row["2010"]))))

df = pd.DataFrame(data = country_decade_growth, columns = ("country", "decade", "growth"))

distance_digit_corr = defaultdict(lambda : defaultdict(lambda : (np.nan, 1)))
for digit in (1, 2, 3, 4):
   df2 = pd.read_csv("country_distances_%d.csv" % digit, sep = "\t")
   df2 = pd.pivot_table(data = df2, index = ["country", "decade"], columns = "measure", values = "distance").reset_index()
   df2 = df2.merge(df, on = ["country", "decade"]).dropna()
   for distance in ("annihil", "gft", "lapl", "mmc", "otp", "spla", "splc", "spls", "count", "eci"):
      distance_digit_corr[distance][digit] = spearmanr(df2["growth"], df2[distance])

print("\\begin{table}")
print("\\centering")
print("\\begin{tabular}{l|llll}")
print("Measure & 1 Digit & 2 Digit & 3 Digit & 4 Digit\\\\")
print("\\hline")
for distance in ("lapl", "mmc", "annihil", "otp", "spls", "spla", "splc", "gft", "count", "eci"):
   digits = {}
   for digit in (1, 2, 3, 4):
      digits[digit] = ("%1.4f$^{***}$" % distance_digit_corr[distance][digit][0] if distance_digit_corr[distance][digit][1] < 0.01 else
                      "%1.4f$^{**}$" % distance_digit_corr[distance][digit][0] if distance_digit_corr[distance][digit][1] < 0.05 else
                      "%1.4f$^{*}$" % distance_digit_corr[distance][digit][0] if distance_digit_corr[distance][digit][1] < 0.1 else
                      "%1.4f" % distance_digit_corr[distance][digit][0])
   print("%s & %s & %s & %s & %s\\\\" % (distance, digits[1], digits[2], digits[3], digits[4]))
print("\\end{tabular}")
print("\\caption{}")
print("\\label{tab:ps-corr}")
print("\\end{table}")
