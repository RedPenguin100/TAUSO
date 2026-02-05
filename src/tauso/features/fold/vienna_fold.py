import ViennaRNA as RNA
import numpy as np

def calculate_avg_mfe_over_sense_region(sequence, sense_start_in_flank, sense_length, window_size=120, step=1):
    sequence = str(sequence).upper().replace('T', 'U')
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i:i + window_size]
        fc = RNA.fold_compound(subseq)
        _, mfe = fc.mfe()
        mfe_per_nt = mfe / window_size

        for j in range(i, i + window_size):
            energy_values[j] += mfe_per_nt
            counts[j] += 1

    avg_energies = np.divide(energy_values, counts, out=np.full_like(energy_values, np.nan), where=counts != 0)

    sense_end_in_flank = sense_start_in_flank + sense_length

    if 0 <= sense_start_in_flank < sequence_length and sense_end_in_flank <= sequence_length:

        return np.nanmean(avg_energies[sense_start_in_flank:sense_end_in_flank])
    else:
        return np.nan
