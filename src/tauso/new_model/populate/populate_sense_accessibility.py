from ..consts_dataframe import SENSE_START
from ...features.feature_names import SENSE_LENGTH
from ...features.rna_access.access_calculator import get_cache, get_sense_with_flanks
from ...features.rna_access.sense_accessibility import compute_sense_accessibility_value

SENSE_AVG_ACCESSIBILITY = 'sense_avg_accessibility'


def populate_sense_accessibility(aso_dataframe, locus_info):
    FLANK_SIZE = 120
    ACCESS_SIZE = 13
    SEED_SIZE = 13
    SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
    ACCESS_WIN_SIZE = 80

    access_cache = get_cache(SEED_SIZES, access_size=ACCESS_SIZE)

    for idx, row in aso_dataframe.iterrows():
        sense_start = row[SENSE_START]
        sense_length = row[SENSE_LENGTH]

        flanked_sense = get_sense_with_flanks(
            str(locus_info.full_mrna), sense_start, sense_length,
            flank_size=FLANK_SIZE
        )
        avg_sense_access = compute_sense_accessibility_value(
            sense_start, sense_length, flank=flanked_sense, flank_size=FLANK_SIZE, access_win_size=ACCESS_WIN_SIZE,
            seed_sizes=SEED_SIZES, access_size=ACCESS_SIZE, cache=access_cache
        )
        aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = avg_sense_access
