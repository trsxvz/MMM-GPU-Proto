import sys

import numpy as np
import pandas as pd

# relative difference accepted, in %
tolerance = 4

# --- Lecture du fichier .res ---
data_res = pd.read_csv("cermet.res", sep=r"\s+", header=None)
data_res.columns = [f"col{i+1}" for i in range(data_res.shape[1])]
data_res['col1'] = data_res['col1'].astype(float)

# --- Lecture du CSV ---
data_csv = pd.read_csv("results/cermet_cast3m.csv", sep=";", header=None)
data_csv[0] = data_csv[0].astype(str).str.strip()  # enlever espaces
data_csv = data_csv[data_csv[0] != '']            # supprimer lignes vides
data_csv[0] = data_csv[0].astype(float)

# --- Fusion sur la première colonne ---
# the times are matched once rounded: in single precision, the program
# writes 0.40000000596046 for 0.4
data_res['time_key'] = data_res['col1'].round(6)
data_csv['time_key'] = data_csv[0].round(6)
merged = pd.merge(
    data_res,
    data_csv,
    on='time_key',
    how='inner',
    suffixes=('_res', '_csv')
).drop(columns='time_key')

# --- Comparaison ---
col_res = 'col7'      # 7ème colonne du .res
col_csv = merged.columns[13]           # 7ème colonne du CSV (index 6 après merge)

# --- Conversion de la colonne CSV en float (nettoyage espaces/virgules si nécessaire) ---
merged[col_csv] = merged[col_csv].astype(str).str.replace(',', '.').str.strip()  # remplacer virgules éventuelles
merged[col_csv] = merged[col_csv].astype(float)

# Calcul de l'écart relatif en %
merged['RelDiff_%'] = (merged[col_csv] - merged[col_res]) / merged[col_res] * 100
#merged['RelDiff_%'] = (merged.iloc[:,col_csv] - merged[col_res]) / merged[col_res] * 100

# Définition du statut
merged['Status'] = np.where(merged['RelDiff_%'].abs() <= tolerance, 'OK', 'MISMATCH')

# Sélection des colonnes pour sortie
output = merged[[merged.columns[0], col_res, col_csv, 'RelDiff_%', 'Status']]
output = merged[[merged.columns[0], col_res, col_csv, 'RelDiff_%', 'Status']]

# Affichage et sauvegarde
output.columns = ['Time', 'MFEM/MGIS', 'CAST3M', 'RelDiff_%', 'Status']
print(output)
output.to_csv("comparison_results.csv", index=False)

# Check if the 'Status' column contains any "MISMATCH": the exit status is
# not null on a mismatch, or if no time is common to both files
if output.empty:
    print("No time common to the results and to the reference!")
    sys.exit(1)
if (output['Status'] == 'MISMATCH').any():
    print("There is at least one MISMATCH in the data!")
    sys.exit(1)
print("Check PASS.")
