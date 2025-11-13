from .access_calculator import AccessCalculator


def compute_sense_accessibility_value(sense_start, sense_length, flank, flank_size, access_win_size, seed_sizes, access_size, min_gc=0, max_gc=100,
                                      gc_ranges=1, cache=None):
    try:
        # Skip invalid rows
        if sense_start == -1:
            print("Sense start bad!")
            return 0

        # Calculate accessibility
        df_access = AccessCalculator.calc(
            flank, access_size,
            min_gc, max_gc, gc_ranges,
            access_win_size, seed_sizes, cache=cache
        )

        flank_start = max(0, sense_start - flank_size)
        sense_start_in_flank = sense_start - flank_start
        sense_end_in_flank = sense_start_in_flank + sense_length

        if 0 <= sense_start_in_flank < len(df_access) and sense_end_in_flank <= len(df_access):
            values = df_access['avg_access'].iloc[sense_start_in_flank:sense_end_in_flank].dropna()
            return values.mean() if not values.empty else None
        else:
            # print(f"sense_start_in_flank: {sense_start_in_flank}")
            # print(f"len(df_access): {len(df_access)}")
            # print(f"sense_end_in_flank: {sense_end_in_flank}")
            return 0
    except Exception as e:
        print(f"Error at {sense_start}, {sense_length} | error: {e}")
        return 0