from pathlib import Path

from Bio import SeqIO

from tauso.util import get_antisense

from .hits import risearch_hits_dataframe

_DATA = Path(__file__).parent / "data"


def get_gfp_seq_and_context():
    gfp_context_path = _DATA / "GFP_context.txt"
    gfp_first_exp_path = _DATA / "GFP_first_exp.fasta"

    gfp_obj = next(SeqIO.parse(str(gfp_first_exp_path), "fasta"))
    gfp_seq = str(gfp_obj.seq.upper())

    with open(str(gfp_context_path), "r") as f:
        gfp_context = f.read().upper()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    return (gfp_seq, gfp_context)


def get_gfp_first_exp(gap=100):
    # TODO: gap should be always 100 in this function
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gfp_ext = gfp_context[gfp_start - gap : gfp_start + len(gfp_seq) + gap]

    return gfp_ext


def get_gfp_second_exp():
    right_gap = 50
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    gfp_ext = gfp_context[gfp_start : gfp_start + len(gfp_seq) + right_gap]

    return gfp_ext


def run_risearch(sample_seq, name_to_seq, *, transpose):
    """The hits for one query against one target, as plain rows."""
    frame = risearch_hits_dataframe(
        [("query", sample_seq)],
        name_to_seq,
        minimum_score=900,
        neighborhood=30,
        transpose=transpose,
    )
    return [
        {
            "query": str(row.query),
            "query_start": int(row.query_start),
            "query_end": int(row.query_end),
            "target": str(row.target),
            "target_start": int(row.target_start),
            "target_end": int(row.target_end),
            "score": int(row.score),
            "energy": float(row.energy),
        }
        for row in frame.itertuples()
    ]


def mfe_scores(hits):
    """The energies of each target, in the order its hits came out."""
    by_target = {}
    for hit in hits:
        by_target.setdefault(hit["target"], []).append(hit["energy"])
    return list(by_target.values())


def test_risearch_gfp_modified_transpose(data_regression):
    # --- 2. Data Preparation ---
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data["metadata"] = {
        "gfp_ontarget_slice": gfp_seq[695:714],
        "sample_seq": sample_seq,
        "sample_antisense": get_antisense(sample_seq),
    }

    # --- 3. Test "Good" Sample ---
    hits_good = run_risearch(sample_seq, name_to_seq, transpose=True)
    mfe_scores_good = mfe_scores(hits_good)

    # Add results to our data snapshot
    regression_data["good_sample_run"] = {
        "hits": hits_good,
        "mfe_scores": mfe_scores_good,
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data["bad_sample_runs"] = []

    for bad_sample in bad_samples:
        hits_bad = run_risearch(bad_sample, name_to_seq, transpose=True)
        mfe_scores_bad = mfe_scores(hits_bad)

        # Append each iteration's result to the list
        regression_data["bad_sample_runs"].append(
            {
                "input_sequence": bad_sample,
                "hits": hits_bad,
                "mfe_scores": mfe_scores_bad,
            }
        )

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)


def test_risearch_gfp_modified_original(data_regression):
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data["metadata"] = {
        "gfp_ontarget_slice": gfp_seq[695:714],
        "sample_seq": sample_seq,
        "sample_antisense": get_antisense(sample_seq),
    }

    # --- 3. Test "Good" Sample ---
    hits_good = run_risearch(sample_seq, name_to_seq, transpose=False)
    mfe_scores_good = mfe_scores(hits_good)

    # Add results to our data snapshot
    regression_data["good_sample_run"] = {
        "hits": hits_good,
        "mfe_scores": mfe_scores_good,
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data["bad_sample_runs"] = []

    for bad_sample in bad_samples:
        hits_bad = run_risearch(bad_sample, name_to_seq, transpose=False)
        mfe_scores_bad = mfe_scores(hits_bad)

        # Append each iteration's result to the list
        regression_data["bad_sample_runs"].append(
            {
                "input_sequence": bad_sample,
                "hits": hits_bad,
                "mfe_scores": mfe_scores_bad,
            }
        )

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)
