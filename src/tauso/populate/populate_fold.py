import uuid

import numpy as np
import pandas as pd
from pandarallel import pandarallel

from ..data.consts import CANONICAL_GENE, SENSE_START, SENSE_LENGTH
from ..features.fold.vienna_fold import calculate_avg_mfe_over_sense_region
from ..features.rna_access.access_calculator import (
    AccessCalculator,
    get_cache,
    get_sense_with_flanks,
)
from ..features.rna_access.sense_accessibility import compute_sense_accessibility_value

PRE_MRNA_SEQUENCE = "pre_mrna_sequence"


def validate_cols_in_df(df, cols):
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")


def populate_mfe_features(df, gene_to_data, n_jobs=1):  # 2. Add n_jobs param
    if n_jobs > 1 or n_jobs == -1:
        pandarallel.initialize(nb_workers=n_jobs, verbose=0)

    required_cols = [CANONICAL_GENE, SENSE_START, SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    settings = [
        (120, 45, 7),
        (120, 45, 15),
        (120, 45, 10),
        (120, 45, 4),
        (120, 35, 15),
        (120, 35, 10),
        (120, 35, 7),
        (120, 35, 4),
        (120, 55, 15),
        (120, 55, 10),
        (120, 55, 7),
        (120, 55, 4),
        (120, 65, 15),
        (120, 65, 10),
        (120, 65, 7),
        (120, 65, 4),
    ]
    feature_names = []
    for setting in settings:
        flank_size, window_size, step = setting
        print(f"Starting MFE population with Flank={flank_size}, Window={window_size}, Step={step}")

        feature_name = f"mfe_win{window_size}_flank{flank_size}_step{step}"
        feature_names.append(feature_name)

        # TODO: extract outside
        def _process_row(row, flank_size=flank_size, window_size=window_size, step=step):
            gene_name = row[CANONICAL_GENE]
            global_start = row[SENSE_START]
            sense_len = row[SENSE_LENGTH]

            if gene_name not in gene_to_data or global_start == -1:
                return np.nan

            locus_info = gene_to_data[gene_name]
            full_mrna = str(locus_info.full_mrna)

            cut_start = max(0, global_start - flank_size)
            cut_end = min(len(full_mrna), global_start + sense_len + flank_size)
            sub_sequence = full_mrna[cut_start:cut_end]

            if len(sub_sequence) < window_size:
                return np.nan

            # Return the score directly
            return calculate_avg_mfe_over_sense_region(
                sequence=sub_sequence,
                sense_start_in_flank=(global_start - cut_start),
                sense_length=sense_len,
                window_size=window_size,
                step=step,
            )

        if n_jobs > 1:
            df[feature_name] = df.parallel_apply(_process_row, axis=1)
        else:
            df[feature_name] = df.apply(_process_row, axis=1)

    return df, feature_names


SENSE_AVG_ACCESSIBILITY = "sense_avg_accessibility"

# --- CONSTANTS ---
FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZE = 13
SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
ACCESS_WIN_SIZE = 80


def populate_sense_accessibility(aso_dataframe, gene_to_data):
    access_cache = get_cache(SEED_SIZES, access_size=ACCESS_SIZE)

    valid_rows = [(idx, row) for idx, row in aso_dataframe.iterrows() if row[SENSE_START] != -1]

    for idx, row in valid_rows:
        sense_start = row[SENSE_START]
        sense_length = row[SENSE_LENGTH]
        full_mrna_seq = gene_to_data[row[CANONICAL_GENE]].full_mrna

        flanked_sense = get_sense_with_flanks(full_mrna_seq, sense_start, sense_length, flank_size=FLANK_SIZE)

        if not flanked_sense or len(flanked_sense) < ACCESS_SIZE:
            aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = 0
            continue

        avg_sense_access = compute_sense_accessibility_value(
            sense_start,
            sense_length,
            flank=flanked_sense,
            flank_size=FLANK_SIZE,
            access_win_size=ACCESS_WIN_SIZE,
            seed_sizes=SEED_SIZES,
            access_size=ACCESS_SIZE,
            cache=access_cache,
        )
        aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = avg_sense_access


def populate_sense_accessibility_batch(
    aso_dataframe,
    gene_to_data,
    batch_size=1000,
    flank_size=FLANK_SIZE,
    access_size=ACCESS_SIZE,
    seed_sizes=SEED_SIZES,
    access_win_size=ACCESS_WIN_SIZE,
    n_jobs=1,
):
    seeds_list = seed_sizes if isinstance(seed_sizes, (list, tuple)) else [seed_sizes]
    seeds_str = "-".join(map(str, seeds_list))
    feature_name = f"access_{flank_size}flank_{access_size}access_{seeds_str}seed_sizes"

    # Setup Output
    df_out = aso_dataframe.copy()

    # Create a purely positional, hidden tracker that ignores the real index
    df_out["_temp_id"] = range(len(df_out))

    valid_mask = df_out[SENSE_START] != -1
    if not valid_mask.any():
        df_out[feature_name] = np.nan
        df_out.drop(columns=["_temp_id"], inplace=True)
        return df_out, feature_name

    # Prepare minimal data for workers using our hidden tracker
    valid_df = df_out.loc[valid_mask, ["_temp_id", SENSE_START, SENSE_LENGTH, CANONICAL_GENE]].copy()
    access_cache = get_cache(seed_sizes, access_size=access_size)

    # 2. Worker Function
    def _process_batch(batch_rows):
        rna_seqs = []
        sense_len_map = {}
        sense_start_map = {}
        batch_uuid = str(uuid.uuid4())

        # Track rows using the hidden positional ID
        for _, row in batch_rows.iterrows():
            current_id = str(row["_temp_id"])

            full_mrna_seq = str(gene_to_data[row[CANONICAL_GENE]].full_mrna)
            current_sense_start = row[SENSE_START]
            flank_start = max(0, current_sense_start - flank_size)
            relative_start = current_sense_start - flank_start

            flanked_sense = get_sense_with_flanks(
                full_mrna_seq,
                current_sense_start,
                row[SENSE_LENGTH],
                flank_size=flank_size,
            )

            if not flanked_sense or len(flanked_sense) < access_size:
                continue

            flanked_sense = flanked_sense.upper().replace("T", "U")
            rna_seqs.append((current_id, flanked_sense))
            sense_len_map[current_id] = row[SENSE_LENGTH]
            sense_start_map[current_id] = relative_start

        if not rna_seqs:
            return pd.DataFrame(columns=["_temp_id", feature_name])

        df = AccessCalculator.calc_new(
            rna_seqs,
            access_size=access_size,
            min_gc=0,
            max_gc=100,
            gc_ranges=1,
            uuid_str=batch_uuid,
            access_win_size=access_win_size,
            access_seed_sizes=seed_sizes,
            cache=access_cache,
        )

        df["pos_in_flank"] = df.groupby("rna_id", as_index=False).cumcount()
        df["sense_length"] = df["rna_id"].map(sense_len_map)
        df["rel_start"] = df["rna_id"].map(sense_start_map)

        mask = (df["pos_in_flank"] >= df["rel_start"]) & (df["pos_in_flank"] < df["rel_start"] + df["sense_length"])

        # Aggregate strictly as a DataFrame
        batch_result = df[mask].groupby("rna_id", as_index=False)["avg_access"].mean()

        # Rename tracker back and ensure it's an integer
        batch_result.rename(columns={"rna_id": "_temp_id", "avg_access": feature_name}, inplace=True)
        batch_result["_temp_id"] = batch_result["_temp_id"].astype(int)

        return batch_result

    # 3. Execution
    valid_df["batch_group"] = np.arange(len(valid_df)) // batch_size

    if n_jobs > 1:
        from pandarallel import pandarallel

        pandarallel.initialize(nb_workers=n_jobs, verbose=0)
        results_df = valid_df.groupby("batch_group").parallel_apply(_process_batch)
    else:
        results_df = valid_df.groupby("batch_group").apply(_process_batch)

    # Clean up structure returned by apply
    results_df = results_df.reset_index(drop=True)

    # 4. Safe Assignment (Dictionary Mapping)
    # Convert results to a fast dictionary {temp_id: calculated_value}
    result_map = dict(zip(results_df["_temp_id"], results_df[feature_name]))

    # Map values directly to the output column. This 100% preserves your original index.
    df_out[feature_name] = df_out["_temp_id"].map(result_map)

    # Cleanup hidden tracker
    df_out.drop(columns=["_temp_id"], inplace=True)

    return df_out, feature_name
