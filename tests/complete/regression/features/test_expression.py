import os

import pandas as pd
import pytest

from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.data.data import get_data_dir
from tauso.dependencies.depmap import load_cell_line_gene_expression
from tauso.populate.populate_context import (
    _SPECIAL_GENES,
    _SPECIAL_TRANSCRIPTS,
    populate_special_gene_expression,
    populate_special_transcript_expression,
    populate_target_dominant_transcript,
    populate_target_expression,
)
from tauso.timer import Timer


def _missing_data(reason):
    """Skip locally, fail in CI.

    These fixtures are the only thing covering the transcript features, so a CI run that
    silently skips them is indistinguishable from one that passes. GitHub Actions sets CI.
    """
    if os.environ.get("CI"):
        pytest.fail(f"{reason} (CI must build this; see the transcript expression step)")
    pytest.skip(reason)


@pytest.fixture(scope="session")
def expression_transcriptomes(base_data, target_genes):
    """Expression data for the target genes plus whatever genes _SPECIAL_GENES names."""
    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")
    depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))
    valid_genes = list(set(target_genes) | set(_SPECIAL_GENES.values()))
    with Timer("Load Expression Transcriptomes"):
        return load_cell_line_gene_expression(depmap_ids, valid_genes, expression_dir=expression_dir)


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_target_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_target_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])


@pytest.mark.skipif(not _SPECIAL_GENES, reason="no special genes configured")
@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_special_gene_expression_regression(mini_sampled_data, expression_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_special_gene_expression(data, expression_transcriptomes)
    dataframe_regression.check(result[feats])


@pytest.fixture(scope="session")
def transcript_transcriptomes(base_data):
    """Transcript-level expression for the transcripts the special-transcript features name."""
    expression_dir = os.path.join(get_data_dir(), "processed_transcript_expression")
    if not os.path.isdir(expression_dir):
        _missing_data("processed_transcript_expression not built; run build-cohort-transcript-expression")
    depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))
    wanted = {name for names in _SPECIAL_TRANSCRIPTS.values() for name in names}
    frames = {}
    with Timer("Load Transcript Transcriptomes"):
        for ach_id in depmap_ids:
            path = os.path.join(expression_dir, f"{ach_id}_transcript_expression.csv")
            if not os.path.exists(path):
                continue
            frame = pd.read_csv(path)
            kept = frame[frame["TranscriptName"].isin(wanted)]
            if not kept.empty:
                frames[ach_id] = kept
    if not frames:
        _missing_data("no transcript expression for the cohort's cell lines")
    return frames


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_special_transcript_expression_regression(mini_sampled_data, transcript_transcriptomes, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_special_transcript_expression(data, transcript_transcriptomes)
    dataframe_regression.check(result[feats])


@pytest.fixture(scope="session")
def target_gene_transcripts(base_data, target_genes):
    """Transcript-level expression for the cohort's target genes."""
    expression_dir = os.path.join(get_data_dir(), "processed_transcript_expression")
    if not os.path.isdir(expression_dir):
        _missing_data("processed_transcript_expression not built; run build-cohort-transcript-expression")
    depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))
    wanted = set(target_genes)
    frames = {}
    with Timer("Load Target Gene Transcripts"):
        for ach_id in depmap_ids:
            path = os.path.join(expression_dir, f"{ach_id}_transcript_expression.csv")
            if not os.path.exists(path):
                continue
            frame = pd.read_csv(path)
            kept = frame[frame["Gene"].isin(wanted)]
            if not kept.empty:
                frames[ach_id] = kept
    if not frames:
        _missing_data("no transcript expression for the cohort's target genes")
    return frames


@pytest.mark.parametrize("mini_sampled_data", [1000], indirect=True)
def test_target_dominant_transcript_regression(mini_sampled_data, target_gene_transcripts, dataframe_regression):
    data = mini_sampled_data.copy()
    result, feats = populate_target_dominant_transcript(data, target_gene_transcripts)
    dataframe_regression.check(result[feats])
