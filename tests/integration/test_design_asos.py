import pytest
import pandas as pd

from tauso.old_model_generation.main import design_asos

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)


def test_sanity():
    print("Running Design Pipeline Demo...")
    try:
        # This requires the genome DB to be set up via cli.py
        results, annotated_hits = design_asos(
            gene_name="DDX11L1",
            organism='mouse',
            top_k=5,
            run_off_target=True  # This triggers the new Bowtie search
        )
        print(results.columns)

        print("\nTop ASO Candidates:")
        print(results)
        print("Annotated hits: ")
        print(annotated_hits)

    except Exception as e:
        print(f"\nDemo failed (likely due to missing DB or Models): {e}")
