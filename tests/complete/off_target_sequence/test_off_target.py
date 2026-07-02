import pandas as pd
import pytest

from notebooks.data.OligoAI.curate_gene_labels import alignment_stats
from tauso.off_target.search import annotate_hits, run_bowtie_search
from tauso.util import get_antisense


def generate_200_test_sequences():
    # Real transcript segments for high-confidence mapping
    seeds = {
        "KRAS": "AUGACUGAAUAUAAACUUGUGGUAGUUGGAGCUGGUGGCGUAGGCAAGAGUGCCUUGACGAUACAGCUAAUUCAGAAUCAUUUUGUGUAC",
        "SMN": "AUUGGUGGAGUGCUCCUCCUCCCCACCUCCCCUAUAUGUGUAGUCCCAGUUUUUUUUUUUUUUUUCUUUUAAAAGUAAAAGUCAUUUAG",
        "APOL1": "GAAUUAACUCCUACAGCAAUGGAGGGAGCUGCUCAUCUGAGAGCAGUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUC",
        "MYO6": "AUGGAGGACGAGAAGCUGAAGUUCGAGUUCGACGAGGAGUUCGAGCUCGAGUUCGACGAGGAGUUCGAGCUCGAGUUCGACGAGGAGUUC",
        "HTT": "AUGGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUA",
    }

    test_data = []

    # Generate 40 tiled sequences for each of the 5 genes (200 total)
    for gene, transcript in seeds.items():
        for i in range(40):
            # Extract a 20bp window
            seq = transcript[i : i + 20]
            if len(seq) == 20:
                test_data.append({"target_gene": gene, "rna_sequence": seq})

    # Add 10 "Scrambled" negatives at the end for controls
    negatives = [
        "GCACUACCAGAGCUAACUCA",  # eGFP
        "UGCAUACUUAACGCGUAUCG",  # Firefly Luciferase
        "CCGGAAUUUCGACGUGAAUA",  # LacZ
        "UACGUCGAUCAGUCGACUGC",  # Synthetic repeat
        "CGAAUACGGAUUACGCGAUG",  # Synthetic repeat
        "AUCGACGCUAGCGUAAUCGG",  # Synthetic repeat
    ]
    for neg in negatives:
        test_data.append({"target_gene": "UNKNOWN", "rna_sequence": neg})

    return test_data


TEST_DATA = generate_200_test_sequences()


def test_alignment_stats_real_genome(tmp_path):
    """A true integration test hitting the real GRCh38 genome and Bowtie binary (no mocks).

    `alignment_stats` is the diagnostic that aligns each ASO to the genome and, per patent
    target_gene, reports the gene its ASOs most often hit. 40 KRAS-derived windows must all
    align to KRAS; scrambled controls must align to nothing.
    """
    input_csv = tmp_path / "test_oligos_input.csv"
    pd.DataFrame(TEST_DATA).to_csv(input_csv, index=False)

    stats = alignment_stats(input_csv=str(input_csv), genome="GRCh38", seq_col="rna_sequence", threads=4)

    assert not stats.empty
    expected_columns = [
        "original_target_gene",
        "total_asos",
        "most_popular_alignment",
        "original_hit_count",
        "popular_hit_count",
    ]
    for col in expected_columns:
        assert col in stats.columns, f"Missing expected column in stats: {col}"

    # We fed 40 KRAS / 40 SMN / 40 APOL1 / 40 MYO6 windows and 6 UNKNOWN controls.
    kras_stats = stats[stats["original_target_gene"] == "KRAS"].iloc[0]
    smn_stats = stats[stats["original_target_gene"] == "SMN"].iloc[0]
    assert kras_stats["total_asos"] == 40
    assert smn_stats["total_asos"] == 40

    # Bowtie resolves the gene: every KRAS window aligns to KRAS.
    assert kras_stats["most_popular_alignment"] == "KRAS"
    assert kras_stats["original_hit_count"] == 40

    # Negative control: scrambled sequences hit nothing.
    unknown_stats = stats[stats["original_target_gene"] == "UNKNOWN"].iloc[0]
    assert unknown_stats["total_asos"] == 6
    assert pd.isna(unknown_stats["most_popular_alignment"])
    assert unknown_stats["popular_hit_count"] == 0


_CLINICAL_ASOS = [
    ("TCTTGGTTACATGAAATCCC", "TTR"),  # Inotersen
    ("GCCTCAGTCTGCTTCGCACC", "APOB"),  # Mipomersen
]


@pytest.mark.parametrize("aso, gene", _CLINICAL_ASOS)
def test_annotate_hits_antisense_maps_to_target_gene(aso, gene):
    hits, _ = run_bowtie_search(aso, max_mismatches=0)
    on_target = annotate_hits(hits).query("gene_name == @gene")
    assert not on_target.empty
    assert (on_target["mismatches"] == 0).any()
    assert on_target["region_type"].isin(["exon", "CDS", "intron", "gene"]).all()


@pytest.mark.parametrize("aso, gene", _CLINICAL_ASOS)
def test_annotate_hits_sense_strand_is_not_counted_as_the_gene(aso, gene):
    sense = get_antisense(aso)
    hits, _ = run_bowtie_search(sense, max_mismatches=0)
    assert hits
    assert gene not in set(annotate_hits(hits)["gene_name"])
