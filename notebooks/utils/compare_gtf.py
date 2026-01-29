import pickle

path_1 = "/home/oni/ASOdesign/scripts/data_genertion/cell_line_expression/"
pkl_good = path_1 + '_gtf_annotations_v34.pkl'
pkl_bad = path_1 + '_gtf_annotations.pkl'

def compare_gtf(pkl1, pkl2, max_diffs=20):
    with open(pkl1, 'rb') as f:
         dict1 = pickle.load(f)
    with open(pkl2, 'rb') as f:
         dict2 = pickle.load(f)

    keys1, keys2 = set(dict1.keys()), set(dict2.keys())

    print(f"File 1 has {len(keys1)} transcripts")
    print(f"File 2 has {len(keys2)} transcripts")

    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1

    if only_in_1:
        print(f"\nTranscripts only in {pkl1}: {len(only_in_1)}")
        print(list(only_in_1)[:max_diffs])
    if only_in_2:
        print(f"\nTranscripts only in {pkl2}: {len(only_in_2)}")
        print(list(only_in_2)[:max_diffs])

    shared = keys1 & keys2
    print(f"\nShared transcripts: {len(shared)}")

    # Check structural differences
    diffs = []
    for k in shared:
        if dict1[k] != dict2[k]:
            diffs.append(k)
            if len(diffs) <= max_diffs:
                print(f"\nDifference in transcript {k}:")
                print(f" - {pkl1}: {dict1[k]}")
                print(f" - {pkl2}: {dict2[k]}")

    print(f"\nFound {len(diffs)} differing transcripts (showing {min(len(diffs), max_diffs)})")

    return only_in_1, only_in_2, diffs


# Example usage:
only_in_good, only_in_bad, diffs = compare_gtf(pkl_good, pkl_bad)