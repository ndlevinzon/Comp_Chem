import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sqlalchemy import create_engine


def main():
    # Connect to the SQLite database
    database_file = 'ligands/histogram_data_filtered.db'
    engine = create_engine(f'sqlite:///{database_file}')

    # Retrieve the updated data from the database
    query = "SELECT * FROM histogram_data"
    df = pd.read_sql_query(query, con=engine)

    subtitle = 'ChEMBL Top 25 Single Proteins with the Most Active Compounds'   # Name of Target
    fig, ax = plt.subplots()                                                    # Set up the bar chart

    # Define the bin edges and width
    bin_edges = np.arange(0, 1.05, 0.05)
    bin_width = bin_edges[1] - bin_edges[0]

    hist, _ = np.histogram(df['Prop_10X'], bins=bin_edges)                      # Calculate the histogram
    hist_normalized = hist / np.sum(hist)                                       # Normalize the histogram  
    
    # Plot the bar chart
    ax.bar(bin_edges[:-1], hist_normalized, width=bin_width, align='edge', alpha=0.5, edgecolor='black')
    ax.set_xlabel('Success Rate \n (10-fold Increase in Potency from Parent with Same B-M Scaffold)')
    ax.set_ylabel('Normalized Frequency')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.suptitle(subtitle)
    plt.title('Distribution of Success Rates for Analogs')
    plt.show()

    print(df)


if __name__ == '__main__':
    main()
