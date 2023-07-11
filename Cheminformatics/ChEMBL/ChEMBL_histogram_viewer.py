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
from scipy.stats import binom, kstest, gaussian_kde
from sqlalchemy import create_engine
from statsmodels.distributions.empirical_distribution import ECDF

def main():
    # List of database files
    database_files = [
        '/mnt/nfs/home/nlevinzon/chembl/total/EC50_ChemblDB.db',
        '/mnt/nfs/home/nlevinzon/chembl/total/IC50_ChemblDB.db',
	'/mnt/nfs/home/nlevinzon/chembl/total/Potency_ChemblDB.db',
        '/mnt/nfs/home/nlevinzon/chembl/total/Ki_ChemblDB.db'
        # Add more database files as needed
    ]

    # Connect to the databases and retrieve the data
    data_frames = []
    for database_file in database_files:
        engine = create_engine(f'sqlite:///{database_file}')
        query = "SELECT * FROM histogram_data"
        df = pd.read_sql_query(query, con=engine)
        data_frames.append(df)

    # Concatenate the data frames into a single data frame
    df = pd.concat(data_frames)

    # Filter out values not within the valid range
    filtered_data = df['Prop_10X'][(df['Prop_10X'] > 0) & (df['Prop_10X'] < np.inf)]

    # Name of Target
    subtitle = 'ChEMBL: All Active Compounds for Single Protein Targets'

    # Set up the bar chart
    fig, ax = plt.subplots()

    # Define the bin edges and width
    bin_edges = np.arange(0, 0.55, 0.025)
    bin_width = bin_edges[1] - bin_edges[0]

    # Calculate the histogram
    hist, _ = np.histogram(filtered_data, bins=bin_edges)

    # Normalize the histogram
    hist_normalized = hist / np.sum(hist)

    # Plot the bar chart
    ax.bar(bin_edges[:-1], hist_normalized, width=bin_width, align='edge', alpha=0.5, edgecolor='black')

    # Perform KDE on the filtered data
    kde = gaussian_kde(filtered_data)

    # Generate the x-values for the binomial distribution
    x = np.linspace(0, 1, 100)

    # Evaluate the KDE at the x-values
    kde_values = kde(x)

    # Normalize the KDE values to have an area of 1
    kde_normalized = kde_values / np.sum(kde_values)

    # Plot the fitted binomial distribution
    ax.plot(x, kde_normalized * 2, color='red', label='Fitted Binomial Distribution')

    # Calculate the empirical CDF of the filtered data
    ecdf = ECDF(filtered_data)

    # Calculate the Kolmogorov-Smirnov statistic and p-value
    kstest_stat, kstest_pvalue = kstest(filtered_data, ecdf)

    ax.set_xlabel('Success Rate')
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.suptitle(subtitle)
    plt.title('Distribution of Success Rates for Analogs (All Standard Types)')

    n = len(filtered_data)
    p = np.mean(filtered_data)
    binomial_params = f'n = {n}, p = {p:.3f}'
    legend_text = f'Goodness-of-Fit\nKS Statistic: {kstest_stat:.3g}\np-value: {kstest_pvalue:.3g}\n{binomial_params}'

    ax.legend(loc='best', title=legend_text)

    # Save the histogram as a high-resolution PNG file
    plt.savefig('/mnt/nfs/home/nlevinzon/chembl/total/total_histogram.png', dpi=300)

    # Show the plot in fullscreen mode
    plt.figure(figsize=(20, 12))
    mng = plt.get_current_fig_manager()
    mng.window.state('zoomed')
    plt.show()
    plt.show()

    print(df)

if __name__ == '__main__':
    main()
