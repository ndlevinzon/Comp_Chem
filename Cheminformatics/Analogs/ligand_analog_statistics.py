import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import ttest_ind
from itertools import combinations

# Read the CSV file
csv_file = 'C:/Users/ndlev/PycharmProjects/shoichet/alpha2a_ampc_summary.csv'
df = pd.read_csv(csv_file)

# Group the data by TARGET
grouped = df.groupby('TARGET')

# Create a single figure for all box plots
plt.figure(figsize=(10, 6))

# Create box plots for each group's PROP_IMRPOV values
box_plot_data = [group_data['PROP_IMPROV'] for target, group_data in grouped]
box_plot = plt.boxplot(box_plot_data, labels=[target for target, group_data in grouped])
plt.xlabel('TARGET')
plt.ylabel('PROP_IMPROV')
plt.title('Box Plot of PROP_IMRPOV for Different Groups \n')
plt.xticks(rotation=45)
plt.tight_layout()

# Perform pairwise t-tests and add colored asterisks for significant results
alpha = 0.05
group_names = [target for target, group_data in grouped]
num_groups = len(group_names)

significant_pairs = []

for i, j in combinations(range(num_groups), 2):
    group1_data = box_plot_data[i]
    group2_data = box_plot_data[j]

    _, p_value = ttest_ind(group1_data, group2_data)

    if p_value < alpha:
        significant_pairs.append((group_names[i], group_names[j], p_value))

# Add legend with asterisk and p-value for significant comparisons
legend_text = '\n'.join([f'{group1} vs {group2}, p={p:.4f}' for group1, group2, p in significant_pairs])
plt.legend([legend_text], loc='upper left')

plt.show()
