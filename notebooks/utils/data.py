import pandas as pd

from notebooks.consts import PROCESSED_OLIGO_CSV_GZ

INDEXED_OLIGO_DATA = 'aso_inhibitions_with_canonical_gene_indexed.csv'

def read_oligo_indexed_data():
    data = pd.read_csv(PROCESSED_OLIGO_CSV_GZ)

    # Check if the column exists
    if 'index_oligo' not in data.columns:
        # Create the column starting at 1
        data['index_oligo'] = range(1, len(data) + 1)

        # Save the file
        data.to_csv('aso_inhibitions_with_canonical_gene_indexed.csv', index=False)
        print("Added 'index_oligo' and saved file.")
    else:
        print("'index_oligo' already exists. No changes made.")

    return data

if __name__ == "__main__":
    read_oligo_indexed_data()
