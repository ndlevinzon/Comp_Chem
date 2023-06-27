# ND Levinzon, UCSF 2023
#  __     ________ _   _ _______
#  \ \   / /  ____| \ | |__   __|/\
#   \ \_/ /| |__  |  \| |  | |  /  \
#    \   / |  __| | . ` |  | | / /\ \
#     | |  | |____| |\  |  | |/ ____ \
#     |_|  |______|_| \_|  |_/_/    \_\

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, kstest
from sqlalchemy import create_engine

def main():
    # Connect to the SQLite database
    database_file = 'ligands/histogram_data_filtered.db'
    engine = create_engine(f'sqlite:///{database_file}')

    # Retrieve the updated data from the database
    query = "SELECT * FROM histogram_data"
    df = pd.read_sql_query(query, con=engine)

    # Filter out values not within the valid range
    filtered_data = df['Prop_10X'][(df['Prop_10X'] > 0) & (df['Prop_10X'] < np.inf)]

    # Calculate z-scores for each data point
    z_scores = (filtered_data - np.mean(filtered_data)) / np.std(filtered_data)

    # Define the threshold for outliers (e.g., z-score > 3)
    outlier_threshold = 3

    # Identify the outliers based on the z-scores
    outliers = filtered_data[z_scores > outlier_threshold]

    # Remove outliers from the filtered_data dataset
    filtered_data = filtered_data[z_scores <= outlier_threshold]

    # Name of Target
    subtitle = 'ChEMBL Top 25 Single Proteins with the Most Active Compounds'

    # Set up the bar chart
    fig, ax = plt.subplots()

    # Define the bin edges and width
    bin_edges = np.arange(0, 0.55, 0.05)
    bin_width = bin_edges[1] - bin_edges[0]

    # Calculate the histogram
    hist, _ = np.histogram(filtered_data, bins=bin_edges)

    # Normalize the histogram
    hist_normalized = hist / np.sum(hist)

    # Plot the bar chart
    ax.bar(bin_edges[:-1], hist_normalized, width=bin_width, align='edge', alpha=0.5, edgecolor='black')

    # Fit the normal distribution to the filtered data
    mu, std = norm.fit(filtered_data)

    # Generate the x-values for the fitted normal distribution
    x = np.linspace(0, 1, 100)

    # Evaluate the fitted normal distribution at the x-values
    fitted_pdf = norm.pdf(x, loc=mu, scale=std)

    # Scale the fitted distribution by the maximum value of the histogram
    scaling_factor = np.max(hist_normalized) / np.max(fitted_pdf)
    fitted_pdf *= scaling_factor

    # Clip the fitted distribution to ensure values are within [0, 1]
    fitted_pdf = np.clip(fitted_pdf, 0, 1)

    # Plot the fitted normal distribution
    ax.plot(x, fitted_pdf, color='red',
            label=f'Fitted Normal Distribution\nEquation: Normal({mu:.3g}, {std:.3g})')

    # Calculate the Kolmogorov-Smirnov statistic and p-value
    kstest_stat, kstest_pvalue = kstest(filtered_data, 'norm', args=(mu, std))

    ax.set_xlabel('Success Rate \n (10-fold Increase in Potency from Parent with Same B-M Scaffold)')
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.suptitle(subtitle)
    plt.title('Distribution of Success Rates for Analogs within ChEMBLE for Single Protein Targets')

    # Print outlier information
    if len(outliers) > 0:
        outlier_info = f"Outliers ({len(outliers)}): {', '.join(outliers.round(3).astype(str))}"
        print(outlier_info)

    ax.legend(loc='best', title=f'Goodness-of-Fit\nKS Statistic: {kstest_stat:.3g}\np-value: {kstest_pvalue:.3g}')
    plt.show()

    print(df)

if __name__ == '__main__':
    main()
