import pandas as pd

from tauso.data.consts import (
    CANONICAL_GENE,
    CELL_LINE,
    CELL_LINE_DEPMAP,
    CELL_LINE_DEPMAP_PROXY,
    CELL_LINE_ORGANISM,
    CELL_LINE_TO_DEPMAP,
    CELL_LINE_TO_DEPMAP_PROXY_DICT,
    CHEMICAL_PATTERN,
    MODIFICATION,
    PS_PATTERN,
    SEQUENCE,
    VOLUME,
)
from tauso.data.data import get_paths
from tauso.features.names import *
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import (
    GeneCoordinateMapper,
    build_gene_sequence_registry,
)
from tauso.populate.calculators.cache import AssetCache
from tauso.populate.calculators.calculator import Calculator
from tauso.util import get_antisense


def get_initial_data(target_mrna, aso_sizes, canonical_name):
    candidates = []
    sense_lengths = []

    for aso_size in aso_sizes:
        for i in range(0, len(target_mrna) - (aso_size - 1)):
            target = target_mrna[i : i + aso_size]
            candidates.append(get_antisense(str(target)))
            sense_lengths.append(aso_size)

    df = pd.DataFrame(
        {
            SEQUENCE: candidates,
            SENSE_LENGTH: sense_lengths,
            CANONICAL_GENE: canonical_name,
        }
    )
    return df


def get_organism_name(genome):
    if genome == "GRCh38":
        return "human"
    elif genome == "GRCm39":
        return "mouse"
    else:
        raise ValueError(f"Unknown genome type: {genome}")


class Transfection:
    ELECTROPORATION = "Electroporation"
    LIPOFECTION = "Lipofection"
    GYMNOSIS = "Gymnosis"
    OTHER = "Other"


def populate_transfection_features(data, transfection_method):
    df = data.copy()

    # Cast the boolean result to an int (1 or 0)
    df["Gymnosis"] = int(transfection_method == Transfection.GYMNOSIS)
    df["Electroporation"] = int(transfection_method == Transfection.ELECTROPORATION)
    df["Lipofection"] = int(transfection_method == Transfection.LIPOFECTION)
    df["Other"] = int(transfection_method == Transfection.OTHER)

    return df, ["Gymnosis", "Electroporation", "Lipofection", "Other"]


def generate_stub_data(
    target_gene: str,
    gene_sequence: str,
    first_n: int = None,
    transfection=Transfection.GYMNOSIS,
    genome: str = "GRCh38"
):
    data = get_initial_data(gene_sequence, aso_sizes=[20], canonical_name=target_gene)

    if first_n is not None:
        data = data[300 : 300 + first_n]  # TODO: change at some point

    data[MODIFICATION] = "MOE/5-methylcytosines/deoxy"
    data[CHEMICAL_PATTERN] = "MMMMMddddddddddMMMMM"

    organism_name = get_organism_name(genome)

    data[CELL_LINE] = "T24"  # TODO: handle empty case
    data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(CELL_LINE_TO_DEPMAP_PROXY_DICT)
    data[CELL_LINE_DEPMAP] = data[CELL_LINE_DEPMAP_PROXY].map(CELL_LINE_TO_DEPMAP)
    data[PS_PATTERN] = 19 * "*"

    data[CELL_LINE_ORGANISM] = organism_name
    data["transfection_method"] = transfection

    data[VOLUME] = 100  # A reasonable stub
    data["cells_per_well"] = 10000

    data.insert(0, "index_generated", range(1, len(data) + 1))
    return data


def generate_aso_features(
    data,
    cache: AssetCache,
    n_jobs=1,
):
    original_columns = set(data.columns)

    print("version is None, not saving features to disk")
    calculator = Calculator(data=data, data_version=None, overwrite=True, cpus=n_jobs, cache=cache)
    calculator.calculate_all()

    final_data = calculator.data

    # Dynamically determine all new features added by the pipeline
    new_features = list(set(final_data.columns) - original_columns)

    return final_data, new_features


if __name__ == "__main__":
    target_gene = "TMSB10"
    genome = "GRCh38"
    gene_to_data = get_locus_to_data_dict(gene_subset=[target_gene], genome=genome)

    paths = get_paths(genome)
    mapper = GeneCoordinateMapper(paths["db"])

    ref_registry = build_gene_sequence_registry(genes=[target_gene], gene_to_data=gene_to_data, mapper=mapper)

    print(
        generate_aso_features(
            target_gene,
            ref_registry,
            mapper,
            gene_to_data,
            genome=genome,
        )
    )
