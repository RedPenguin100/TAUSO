import os
from pathlib import Path

import pandas as pd

from tauso.algorithms.genomic_context_windows import (
    add_external_mrna_and_context_columns,
)
from tauso.common.gtf import filter_gtf_genes
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
from tauso.data.data import get_data_dir, get_paths, load_db
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_expression
from tauso.features.context.mrna_halflife import (
    populate_mrna_halflife_features,
)

# from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
# from tauso.features.hybridization_off_target.off_target_feature import populate_off_target_general
from tauso.features.context.ribo_seq import add_genomic_coordinates
from tauso.features.hybridization_off_target.add_off_target_feat import AggregationMethod
from tauso.features.hybridization_off_target.common import get_general_expression_of_genes
from tauso.features.hybridization_off_target.off_target_feature import (
    populate_off_target_general,
    populate_off_target_specific,
)
from tauso.features.hybridization_off_target.off_target_specific_gene import (
    off_target_specific_seq_pandarallel,
    on_target_total_hybridization,
)
from tauso.features.names import *
from tauso.features.rbp.pwm_helper import build_rbp_expression_matrix
from tauso.features.rbp.rbp_annotations import rbp_role_map_strict
from tauso.features.rbp.RBP_features import create_positional_sequence_columns
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import (
    GeneCoordinateMapper,
    build_gene_sequence_registry,
)
from tauso.populate.calculators.cache import AssetCache
from tauso.populate.calculators.calculator import Calculator
from tauso.populate.populate_codon_usage import populate_cai, populate_enc, populate_tai
from tauso.populate.populate_context import populate_ribo_seq
from tauso.populate.populate_fold import populate_mfe_features, populate_sense_accessibility_batch
from tauso.populate.populate_hybridization import populate_hybridization
from tauso.populate.populate_modification import populate_modifications
from tauso.populate.populate_rbp import (
    populate_complexity_features,
    populate_functional_features,
    populate_rbp_affinity_features,
    populate_rbp_interaction_features,
    populate_rbp_region,
)
from tauso.populate.populate_rnase import populate_rnase_features
from tauso.populate.populate_sequence import populate_sequence_features, populate_sequence_one_hot_encoded
from tauso.populate.populate_sequence_chemistry import populate_sequence_chemistry_features
from tauso.populate.populate_structure import get_populated_df_with_structure_features
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
    target_gene: str, gene_sequence: str, first_n: int = None, transfection=Transfection.GYMNOSIS, genome: str = "GRCh38"
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

    from tauso.populate.calculators.calculator import Calculator

    print("version is None, not saving features to disk")
    calculator = Calculator(data=data, data_version=None, overwrite=True, cpus=n_jobs, cache=cache)
    calculator.calculate_all()

    final_data = calculator.data

    # Dynamically determine all new features added by the pipeline
    new_features = list(set(final_data.columns) - original_columns)

    return final_data, new_features


def generate_aso_features_bak(
    target_gene,
    ref_registry,
    mapper,
    attract_data,
    half_life_provider,
    gene_to_data,
    genome="GRCh38",
    transfection=Transfection.GYMNOSIS,
    n_jobs=1,
    first_n=None,
):
    data = get_initial_data(gene_to_data[target_gene].full_mrna, aso_sizes=[20], canonical_name=target_gene)

    # Fast data
    if not first_n is None:
        data = data[300:310]

    data[MODIFICATION] = "MOE/5-methylcytosines/deoxy"
    data[CHEMICAL_PATTERN] = "MMMMMddddddddddMMMMM"

    organism_name = get_organism_name(genome)
    data[CELL_LINE_ORGANISM] = organism_name
    data[CELL_LINE_DEPMAP] = "ACH-000018"  # TODO: handle empty case
    data[CELL_LINE] = "T24"  # TODO: handle empty case

    data[VOLUME] = 100  # A reasonable stub
    data["cells_per_well"] = 10000

    data["is_all_ps"] = 1  # TODO: deduce
    data["rnase_expression"] = 1  # TODO: deduce
    data["target_expression"] = 1  # TODO: deduce
    data["max_consecutive_PO"] = 1  # TODO: deduce
    data["po_percentage"] = 0  # TODO: deduce
    data["ps_end_score"] = 38

    data.insert(0, "index_generated", range(1, len(data) + 1))

    calculator = Calculator(data=data, data_version="generated", overwrite=True, cpus=8)

    features = []
    data = get_populated_df_with_structure_features(data, [target_gene], gene_to_data, use_mask=False)
    features.extend(
        [
            SENSE_START,
            SENSE_LENGTH,
            SENSE_START_FROM_END,
            SENSE_EXON,
            SENSE_INTRON,
            SENSE_UTR,
        ]
    )

    data, transfection_features = populate_transfection_features(data, transfection)
    features.extend(transfection_features)

    data, hybridization_features = populate_hybridization(data, n_cores=n_jobs)
    features.extend(hybridization_features)

    data, modification_features = populate_modifications(data, n_cores=n_jobs)
    features.extend(modification_features)

    data = add_genomic_coordinates(data, mapper)
    data, ribo_features = populate_ribo_seq(organism_name, data)  # TODO: generalize beyond human
    features.extend(ribo_features)

    FLANK_SIZES_PREMRNA = [20, 30, 40, 50, 60, 70]
    CDS_WINDOWS = [20, 30, 40, 50, 60, 70]

    # Run the optimized context generator
    data = add_external_mrna_and_context_columns(
        df=data,
        mapper=mapper,
        gene_registry=ref_registry,
        flank_sizes_premrna=FLANK_SIZES_PREMRNA,
        flank_sizes_cds=CDS_WINDOWS,
    )

    print("Genomic context added.")

    # Codon usage bias

    data, tai_features = populate_tai(data, CDS_WINDOWS, ref_registry)
    features.extend(tai_features)

    data, cai_features = populate_cai(data, CDS_WINDOWS, ref_registry, n_jobs=n_jobs)
    features.extend(cai_features)

    data, enc_features = populate_enc(data, CDS_WINDOWS, ref_registry, n_jobs=n_jobs)
    features.extend(enc_features)

    data, sequence_features = populate_sequence_features(data, cpus=n_jobs)
    features.extend(sequence_features)

    data, sequence_one_hot_encoded_features = populate_sequence_one_hot_encoded(data, cpus=n_jobs)
    features.extend(sequence_one_hot_encoded_features)

    data, seq_chem_features = populate_sequence_chemistry_features(data, cpus=n_jobs)
    features.extend(seq_chem_features)

    # Off-target
    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")
    depmap_ids = list(set(data[CELL_LINE_DEPMAP]))

    db = load_db()

    valid_genes = filter_gtf_genes(db, filter_mode="non_mt")

    # 2. Run the Optimized Loader
    #    This loads the DB/FASTA once, filters for valid genes, and enriches sequences.
    transcriptomes = load_cell_line_gene_expression(
        depmap_ids,
        valid_genes,
        expression_dir=expression_dir,  # Directory containing expression CSVs
    )

    mean_exp_data = get_general_expression_of_genes(
        Path(get_data_dir()) / "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv", valid_genes
    )
    transcriptomes["general"] = mean_exp_data

    gene_to_data_full = get_locus_to_data_dict(include_introns=True)

    top_n_list = [50, 100, 200]
    cutoff_list = [800, 1000, 1200]

    for n in top_n_list:
        for cutoff in cutoff_list:
            data, off_target_general_features = populate_off_target_general(
                ASO_df=data,
                gene_to_data=gene_to_data_full,
                cell_line2data=transcriptomes,
                top_n_list=[n],
                cutoff_list=[cutoff],
                method=AggregationMethod.ARTM,
                n_cores=n_jobs,
            )
            features.extend(off_target_general_features)

            data, off_target_specific_features = populate_off_target_specific(
                ASO_df=data,
                gene_to_data=gene_to_data,
                cell_line2data=transcriptomes,
                top_n_list=[n],
                cutoff_list=[cutoff],
                method=AggregationMethod.ARTM,
                n_cores=n_jobs,
            )
            features.extend(off_target_specific_features)

    cutoffs = [0, 1200]

    for cutoff in cutoffs:
        data, feature_name = off_target_specific_seq_pandarallel(
            data, "RNASEH1", gene_to_data_full, cutoff=cutoff, n_jobs=n_jobs, verbose=True
        )
        features.append(feature_name)

        data, feature_name = off_target_specific_seq_pandarallel(
            data, "ACTB", gene_to_data_full, cutoff=cutoff, n_jobs=n_jobs, verbose=True
        )
        features.append(feature_name)

        data, feature_name = on_target_total_hybridization(
            data, gene_to_data_full, cutoff=cutoff, n_jobs=n_jobs, verbose=True
        )
        features.append(feature_name)

    # RBP features
    rbp_map, pwm_db = attract_data
    expression_matrix = build_rbp_expression_matrix(df=data, pwm_db=pwm_db)

    for flank_size in [50]:
        window_col = f"flank_sequence_{flank_size}"

        data, individual_features = populate_rbp_interaction_features(
            data,
            rbp_map,
            pwm_db,
            expression_matrix,
            gene_to_data,
            window_col,
            n_jobs=n_jobs,
        )

        data, global_features = populate_complexity_features(
            data, individual_features, suffix=str(flank_size), type="expression"
        )

    features.extend(individual_features)
    features.extend(global_features)

    for flank_size in [50]:
        window_col = f"flank_sequence_{flank_size}"

        data, individual_features = populate_rbp_affinity_features(
            data, rbp_map, pwm_db, gene_to_data, window_col, n_jobs=32
        )

        data, global_features = populate_complexity_features(
            data, individual_features, suffix=str(flank_size), type="generic"
        )

    features.extend(individual_features)
    features.extend(global_features)

    data, functional_features = populate_functional_features(data, rbp_role_map_strict, flank_size=50)
    features.extend(functional_features)

    # Most important RBP: RNASEH-1
    data, rnase_features = populate_rnase_features(data)
    features.extend(rnase_features)

    data = create_positional_sequence_columns(data, "flank_sequence_50", flank_size=50)

    # 3. LOOP REGIONS
    for region in ["left", "core", "right"]:
        # Fix Memory Fragmentation
        data = data.copy()

        # --- CALL THE FUNCTION ---w
        data, new_features = populate_rbp_region(
            df=data,
            region_name=region,
            rbp_map=rbp_map,
            pwm_db=pwm_db,
            expr_matrix=expression_matrix,
            role_map=rbp_role_map_strict,
            gene_to_data=gene_to_data,
        )

        features.extend(new_features)

    data, half_life_features = populate_mrna_halflife_features(data, half_life_provider)
    features.extend(half_life_features)

    ####################
    # Fold features    #

    data, mfe_features = populate_mfe_features(data, gene_to_data, n_jobs=n_jobs)
    features.extend(mfe_features)

    configurations = [
        {"flank": 120, "access": 20, "seeds": [13]},
        {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
        {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
        {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
    ]

    feature_names = []
    for config in configurations:
        c_flank = config["flank"]
        c_access = config["access"]
        c_seeds = config["seeds"]

        print(f"Running: Flank={c_flank}, Access={c_access}, Seeds={c_seeds})...")

        data, access_feature = populate_sense_accessibility_batch(
            data, gene_to_data, batch_size=1000, flank_size=c_flank, access_size=c_access, seed_sizes=c_seeds, n_jobs=16
        )
        features.extend(access_feature)

    return data, features


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
