import os

import pandas as pd

from tauso.algorithms.genomic_context_windows import (
    add_external_mrna_and_context_columns,
)
from tauso.common.gtf import filter_gtf_genes
from tauso.data.consts import (
    CANONICAL_GENE,
    CELL_LINE,
    CELL_LINE_DEPMAP,
    CELL_LINE_ORGANISM,
    CHEMICAL_PATTERN,
    MODIFICATION,
    SEQUENCE,
)
from tauso.data.data import get_data_dir, get_paths, load_db
from tauso.features.context.mrna_halflife import (
    HalfLifeProvider,
    load_halflife_mapping,
    populate_mrna_halflife_features,
)
from tauso.features.context.ribo_seq import add_genomic_coordinates
from tauso.features.names import *
from tauso.features.rbp.load_rbp import load_attract_data
from tauso.features.rbp.pwm_helper import build_rbp_expression_matrix
from tauso.features.rbp.rbp_annotations import rbp_role_map_strict
from tauso.features.rbp.RBP_features import create_positional_sequence_columns
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import (
    GeneCoordinateMapper,
    build_gene_sequence_registry,
)
from tauso.populate.populate_codon_usage import populate_cai, populate_enc, populate_tai
from tauso.populate.populate_context import populate_ribo_seq
from tauso.populate.populate_fold import populate_mfe_features
from tauso.populate.populate_hybridization import populate_hybridization
from tauso.populate.populate_modification import populate_modifications
from tauso.populate.populate_rbp import (
    populate_complexity_features,
    populate_functional_features,
    populate_rbp_affinity_features,
    populate_rbp_interaction_features,
    populate_rbp_region,
)
from tauso.populate.populate_sequence import populate_sequence_features
from tauso.populate.populate_structure import get_populated_df_with_structure_features
from tauso.util import get_antisense


def get_initial_data(target_mrna, aso_sizes, canonical_name):
    candidates = []
    sense_lengths = []

    for aso_size in aso_sizes:
        for i in range(0, len(target_mrna) - (aso_size - 1)):
            target = target_mrna[i: i + aso_size]
            candidates.append(get_antisense(str(target)))
            sense_lengths.append(aso_size)

    df = pd.DataFrame({SEQUENCE: candidates, SENSE_LENGTH: sense_lengths,
                       CANONICAL_GENE: canonical_name})
    return df


def get_organism_name(genome):
    if genome == 'GRCh38':
        return 'human'
    elif genome == 'GRCm39':
        return 'mouse'
    else:
        raise ValueError('Unknown genome type: {}'.format(genome))


class Transfection:
    ELECTROPORATION = 'ELECTROPORATION'
    LIPOFECTION = 'LIPOFECTION'
    GYMNOSIS = 'GYMNOSIS'
    OTHER = 'OTHER'


def populate_transfection_features(data, transfection_method):
    df = data.copy()

    # Cast the boolean result to an int (1 or 0)
    df['Gymnosis'] = int(transfection_method == Transfection.GYMNOSIS)
    df['Electroporation'] = int(transfection_method == Transfection.ELECTROPORATION)
    df['Lipofection'] = int(transfection_method == Transfection.LIPOFECTION)
    df['Other'] = int(transfection_method == Transfection.OTHER)

    return df, ['Gymnosis', 'Electroporation', 'Lipofection', 'Other']


def generate_aso_features(target_gene, ref_registry, mapper, attract_data, half_life_provider, gene_to_data,
                          genome='GRCh38', transfection=Transfection.GYMNOSIS, n_jobs=1):
    data = get_initial_data(gene_to_data[target_gene].full_mrna, aso_sizes=[20], canonical_name=target_gene)

    # Fast data
    data = data[300:310]

    data[MODIFICATION] = 'MOE/5-methylcytosines/deoxy'
    data[CHEMICAL_PATTERN] = 'MMMMMddddddddddMMMMM'

    organism_name = get_organism_name(genome)
    data[CELL_LINE_ORGANISM] = organism_name
    data[CELL_LINE_DEPMAP] = "ACH-000018"  # TODO: handle empty case
    data[CELL_LINE] = "T24"  # TODO: handle empty case

    features = []
    data = get_populated_df_with_structure_features(data, [target_gene], gene_to_data, use_mask=False)
    features.extend([SENSE_START, SENSE_LENGTH, SENSE_START_FROM_END, SENSE_EXON, SENSE_INTRON, SENSE_UTR])

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
        flank_sizes_cds=CDS_WINDOWS
    )

    print("Genomic context added.")

    # Codon usage bias

    data, tai_features = populate_tai(data, CDS_WINDOWS, ref_registry)
    features.extend(tai_features)

    data, cai_features = populate_cai(data, CDS_WINDOWS, ref_registry)
    features.extend(cai_features)

    data, enc_features = populate_enc(data, CDS_WINDOWS, ref_registry)
    features.extend(enc_features)

    data, sequence_features = populate_sequence_features(data, cpus=32)
    features.extend(sequence_features)

    # Off-target
    data_dir = get_data_dir()
    expression_dir = os.path.join(data_dir, "processed_expression")
    depmap_ids = list(set(data[CELL_LINE_DEPMAP]))

    db = load_db()

    valid_genes = filter_gtf_genes(db, filter_mode='non_mt')

    # # 2. Run the Optimized Loader
    # #    This loads the DB/FASTA once, filters for valid genes, and enriches sequences.
    # transcriptomes = load_cell_line_gene_expression(
    #     depmap_ids,
    #     valid_genes,
    #     expression_dir=expression_dir,  # Directory containing expression CSVs
    # )
    #
    # mean_exp_data = get_general_expression_of_genes(
    #     Path(get_data_dir()) / 'OmicsExpressionTPMLogp1HumanAllGenesStranded.csv', valid_genes)
    # transcriptomes['general'] = mean_exp_data

    rbp_map, pwm_db = attract_data
    expression_matrix = build_rbp_expression_matrix(df=data, pwm_db=pwm_db)

    for flank_size in [50]:
        window_col = f'flank_sequence_{flank_size}'

        data, individual_features = populate_rbp_interaction_features(
            data, rbp_map, pwm_db, expression_matrix, gene_to_data,
            window_col, n_jobs=n_jobs
        )

        data, global_features = populate_complexity_features(
            data,
            individual_features,
            suffix=str(flank_size), type='expression'
        )

    features.extend(individual_features)
    features.extend(global_features)

    for flank_size in [50]:
        window_col = f'flank_sequence_{flank_size}'

        data, individual_features = populate_rbp_affinity_features(
            data, rbp_map, pwm_db, gene_to_data,
            window_col, n_jobs=32
        )

        data, global_features = populate_complexity_features(
            data,
            individual_features,
            suffix=str(flank_size), type='generic'
        )

    features.extend(individual_features)
    features.extend(global_features)

    data, functional_features = populate_functional_features(
        data,
        rbp_role_map_strict,
        flank_size=50
    )

    data = create_positional_sequence_columns(data, "flank_sequence_50", flank_size=50)

    # 3. LOOP REGIONS
    for region in ['left', 'core', 'right']:
        # Fix Memory Fragmentation
        data = data.copy()

        # --- CALL THE FUNCTION ---
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

    data, mfe_features = populate_mfe_features(data, gene_to_data, n_jobs=n_jobs)

    return data, features


if __name__ == '__main__':
    target_gene = 'TMSB10'
    genome = 'GRCh38'
    gene_to_data = get_locus_to_data_dict(gene_subset=[target_gene], genome=genome)

    paths = get_paths(genome)
    mapper = GeneCoordinateMapper(paths['db'])

    ref_registry = build_gene_sequence_registry(
        genes=[target_gene],
        gene_to_data=gene_to_data,
        mapper=mapper
    )

    attract_data = load_attract_data()
    mapping = load_halflife_mapping()
    half_life_provider = HalfLifeProvider(mapping)

    print(generate_aso_features(target_gene, ref_registry, mapper, attract_data, half_life_provider, gene_to_data,
                                genome=genome))
