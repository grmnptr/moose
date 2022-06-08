import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14


failure_probability = 349.34465596472

original_db = pd.read_csv("concatinated_realistic.csv", names=["T_avg"], usecols=[0])
nonoptimized_db = pd.read_csv("original_aout_acc.csv", skiprows=2, names=["T_avg"])
optimized_db = pd.read_csv("Main_GP_4D_aout_acc_5001.csv", skiprows=2, names=["T_avg"])

nonoptimized_db_count = pd.read_csv("original_cout.csv", skiprows=1, names=["flag"], usecols=[1])
optimized_db_count = pd.read_csv("Main_GP_4D_cout.csv", skiprows=1, names=["flag"], usecols=[1])

combined_nonoptimized = pd.concat([nonoptimized_db, nonoptimized_db_count], axis=1)
combined_optimized = pd.concat([optimized_db, optimized_db_count], axis=1)

filtered_nonoptimized = combined_nonoptimized[combined_nonoptimized.flag == False]
filtered_optimized = combined_optimized[combined_optimized.flag == False]

failed_nonoptimized = filtered_nonoptimized[filtered_nonoptimized.T_avg > failure_probability].count()
failed_optimized = filtered_optimized[filtered_optimized.T_avg > failure_probability].count()

fp_nonoptimized = failed_nonoptimized.T_avg/len(filtered_nonoptimized)
fp_optimized = failed_optimized.T_avg/len(filtered_optimized)
print("HF evaluations: ", 5000 - len(filtered_nonoptimized), "Failure probability for nonoptimized: ", fp_nonoptimized)
print("HF evaluations: ", 5000 - len(filtered_optimized), "Failure probability for optimized: ", fp_optimized)

print("COV reference: ", np.sqrt((0.96)*0.04**2+0.04*(0.96)**2)/np.sqrt(len(original_db))/0.04)
print("COV for nonoptimized: ", np.sqrt((1-fp_nonoptimized)*fp_nonoptimized**2+fp_nonoptimized*(1-fp_nonoptimized)**2)/np.sqrt(len(filtered_nonoptimized))/fp_nonoptimized)
print("COV for optimized: ", np.sqrt((1-fp_optimized)*fp_optimized**2+fp_optimized*(1-fp_optimized)**2)/np.sqrt(len(filtered_optimized))/fp_optimized)

dbs = [original_db, filtered_nonoptimized, filtered_optimized]
labels = ["High-Fidelity Model", "Non-optimized Active Learner", "Optimized Active Learner"]

fig, ax = plt.subplots()

for dbi in range(len(dbs)):
    sns.histplot(dbs[dbi].T_avg, ax=ax, element="step", fill=False,
                 bins=range(240, 370, 5), stat="probability", label=labels[dbi])

plt.plot([failure_probability, failure_probability], [0,1], '--k', label="Failure Threshold (349.345 K)")
ax.set_ylim([0,0.07])
ax.set_xlabel("Average Temperature (K)")

plt.legend()
plt.savefig("histograms.pdf")
