import pandas as pd

from notebooks.consts import ORIGINAL_OLIGO_CSV_WITH_CANONICAL, OLIGO_CSV_INDEXED


def read_oligo_indexed_data():
    data = pd.read_csv(ORIGINAL_OLIGO_CSV_WITH_CANONICAL)

    # Check if the column exists
    if "index_oligo" not in data.columns:
        # Create the column starting at 1
        data["index_oligo"] = range(1, len(data) + 1)

        # Save the file
        data.to_csv(OLIGO_CSV_INDEXED, index=False)
        print("Added 'index_oligo' and saved file.")
    else:
        print("'index_oligo' already exists. No changes made.")

    return data


if __name__ == "__main__":
    read_oligo_indexed_data()
