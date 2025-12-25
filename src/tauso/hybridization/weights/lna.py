# Parameters from Owczarzy et al. (2011) - Table 3
# LNA Single Mismatches (Differential Parameters)
# Format: { '5-Top-3/3-Bot-5': {'dH': kcal/mol, 'dS': cal/mol*K} }
# These values must be ADDED to the standard DNA:DNA Mismatch parameters (SantaLucia 2004).

LNA_DNA_WEIGHTS = {
    '+A+A/TT': {'dH': -2.091, 'dS': -4.975},
    '+A+C/TG': {'dH': -2.989, 'dS': -6.563},
    '+A+G/TC': {'dH': -4.993, 'dS': -10.607},
    '+A+T/TA': {'dH': -7.503, 'dS': -20.350},

    '+C+A/GT': {'dH': -5.677, 'dS': -12.798},
    '+C+C/GG': {'dH': -7.399, 'dS': -16.475},
    '+C+G/GC': {'dH': -3.958, 'dS': -8.039},
    '+C+T/GA': {'dH': -7.937, 'dS': -20.218},

    '+G+A/CT': {'dH': -5.759, 'dS': -12.897},
    '+G+C/CG': {'dH': -6.309, 'dS': -16.338},
    '+G+G/CC': {'dH': -5.022, 'dS': -9.773},
    '+G+T/CA': {'dH': -8.961, 'dS': -23.458},

    '+T+A/AT': {'dH': -3.118, 'dS': -4.808},
    '+T+C/AG': {'dH': -0.966, 'dS': 0.665},   # Note: Positive entropy delta
    '+T+G/AC': {'dH': -1.546, 'dS': 0.109},   # Note: Positive entropy delta
    '+T+T/AA': {'dH': -2.519, 'dS': -5.483},

    # --- +A·A Mismatch ---
    '+A+A/AT': {'dH': 4.074, 'dS': 9.091},
    '+A+C/AG': {'dH': 6.033, 'dS': 15.078},
    '+A+G/AC': {'dH': 2.951, 'dS': 7.993},
    '+A+T/AA': {'dH': 2.151, 'dS': 2.886},
    '+A+A/TA': {'dH': 3.671, 'dS': 7.040},
    '+C+A/GA': {'dH': 2.622, 'dS': 5.037},
    '+G+A/CA': {'dH': -0.358, 'dS': -1.776},
    '+T+A/AA': {'dH': 9.274, 'dS': 24.746},

    # --- +C·C Mismatch ---
    '+C+A/CT': {'dH': 10.718, 'dS': 27.450},
    '+C+C/CG': {'dH': 9.127, 'dS': 21.726},
    '+C+G/CC': {'dH': -0.303, 'dS': -4.825},
    '+C+T/CA': {'dH': 5.747, 'dS': 10.483},
    '+A+C/TC': {'dH': 9.465, 'dS': 20.997},
    '+C+C/GC': {'dH': -1.522, 'dS': -7.124},
    '+G+C/CC': {'dH': 5.033, 'dS': 9.503},
    '+T+C/AC': {'dH': 12.314, 'dS': 31.458},

    # --- +G·G Mismatch ---
    '+G+A/GT': {'dH': 5.280, 'dS': 12.813},
    '+G+C/GG': {'dH': 1.661, 'dS': 2.616},
    '+G+G/GC': {'dH': 2.851, 'dS': 7.392},
    '+G+T/GA': {'dH': -0.591, 'dS': -4.911},
    '+A+G/TG': {'dH': 2.820, 'dS': 5.574},
    '+C+G/GG': {'dH': 6.159, 'dS': 15.042},
    '+G+G/CG': {'dH': -5.505, 'dS': -16.121},
    '+T+G/AG': {'dH': 5.725, 'dS': 13.414},

    # --- +T·T Mismatch ---
    '+T+A/TT': {'dH': 3.456, 'dS': 9.151},
    '+T+C/TG': {'dH': 3.813, 'dS': 8.680},
    '+T+G/TC': {'dH': 2.154, 'dS': 6.071},
    '+T+T/TA': {'dH': 0.203, 'dS': -2.849},
    '+A+T/TT': {'dH': 2.993, 'dS': 6.093},
    '+C+T/GT': {'dH': -0.376, 'dS': -1.962},
    '+G+T/CT': {'dH': 1.159, 'dS': 1.778},
    '+T+T/AT': {'dH': 5.849, 'dS': 15.145},

    # --- +A·C Mismatch ---
    '+A+A/CT': {'dH': 6.538, 'dS': 16.649},
    '+A+C/CG': {'dH': 6.641, 'dS': 15.889},
    '+A+G/CC': {'dH': 1.251, 'dS': 2.927},
    '+A+T/CA': {'dH': 3.637, 'dS': 6.295},
    '+A+A/TC': {'dH': 5.822, 'dS': 12.112},
    '+C+A/GC': {'dH': 2.632, 'dS': 5.748},
    '+G+A/CC': {'dH': -0.277, 'dS': -2.365},
    '+T+A/AC': {'dH': 9.890, 'dS': 26.265},

    # --- +C·A Mismatch ---
    '+C+A/AT': {'dH': -1.344, 'dS': -6.973},
    '+C+C/AG': {'dH': 4.239, 'dS': 8.696},
    '+C+G/AC': {'dH': 0.755, 'dS': -0.116},
    '+C+T/AA': {'dH': 4.411, 'dS': 8.483},
    '+A+C/TA': {'dH': 9.153, 'dS': 21.897},
    '+C+C/GA': {'dH': -4.714, 'dS': -15.655},
    '+G+C/CA': {'dH': -2.858, 'dS': -11.329},
    '+T+C/AA': {'dH': 6.481, 'dS': 15.177},

    # --- +A·G Mismatch ---
    '+A+A/GT': {'dH': 10.093, 'dS': 26.574},
    '+A+C/GG': {'dH': -0.053, 'dS': -0.272},
    '+A+G/GC': {'dH': 6.636, 'dS': 18.468},
    '+A+T/GA': {'dH': -0.218, 'dS': -3.666},
    '+A+A/TG': {'dH': 5.937, 'dS': 13.187},
    '+C+A/GG': {'dH': -0.212, 'dS': -1.079},
    '+G+A/CG': {'dH': 0.325, 'dS': 0.539},
    '+T+A/AG': {'dH': 10.407, 'dS': 28.456},

    # --- +G·A Mismatch ---
    '+G+A/AT': {'dH': 5.286, 'dS': 12.798},
    '+G+C/AG': {'dH': 0.669, 'dS': -0.947},
    '+G+G/AC': {'dH': 5.846, 'dS': 16.029},
    '+G+T/AA': {'dH': -0.115, 'dS': -3.913},
    '+A+G/TA': {'dH': 1.109, 'dS': -0.148},
    '+C+G/GA': {'dH': 6.640, 'dS': 16.612},
    '+G+G/CA': {'dH': -4.898, 'dS': -14.756},
    '+T+G/AA': {'dH': 8.834, 'dS': 22.260},

    # --- +C·T Mismatch ---
    '+C+A/TT': {'dH': 8.882, 'dS': 22.121},
    '+C+C/TG': {'dH': 5.284, 'dS': 11.900},
    '+C+G/TC': {'dH': 0.237, 'dS': -2.115},
    '+C+T/TA': {'dH': 2.017, 'dS': 0.827},
    '+A+C/TT': {'dH': 7.708, 'dS': 17.122},
    '+C+C/GT': {'dH': -2.299, 'dS': -8.603},
    '+G+C/CT': {'dH': 0.738, 'dS': -1.956},
    '+T+C/AT': {'dH': 10.273, 'dS': 26.168},

    # --- +T·C Mismatch ---
    '+T+A/CT': {'dH': 1.715, 'dS': 3.953},
    '+T+C/CG': {'dH': 9.651, 'dS': 23.756},
    '+T+G/CC': {'dH': 1.287, 'dS': 2.572},
    '+T+T/CA': {'dH': 5.503, 'dS': 10.829},
    '+A+T/TC': {'dH': 6.567, 'dS': 14.599},
    '+C+T/GC': {'dH': 0.932, 'dS': 0.000},
    '+G+T/CC': {'dH': 2.547, 'dS': 5.757},
    '+T+T/AC': {'dH': 8.111, 'dS': 20.754},

    # --- +G·T Mismatch ---
    '+G+A/TT': {'dH': 2.649, 'dS': 6.802},
    '+G+C/TG': {'dH': -5.143, 'dS': -15.748},
    '+G+G/TC': {'dH': -0.110, 'dS': 1.551},
    '+G+T/TA': {'dH': -5.813, 'dS': -17.641},
    '+A+G/TT': {'dH': 0.670, 'dS': 0.214},
    '+C+G/GT': {'dH': -4.262, 'dS': -12.230},
    '+G+G/CT': {'dH': -6.622, 'dS': -17.610},
    '+T+G/AT': {'dH': 1.797, 'dS': 4.589},

    # --- +T·G Mismatch ---
    '+T+A/GT': {'dH': 2.588, 'dS': 7.261},
    '+T+C/GG': {'dH': -1.598, 'dS': -4.206},
    '+T+G/GC': {'dH': 3.981, 'dS': 11.635},
    '+T+T/GA': {'dH': 3.377, 'dS': 6.507},
    '+A+T/TG': {'dH': 4.836, 'dS': 11.566},
    '+C+T/GG': {'dH': -3.596, 'dS': -9.732},
    '+G+T/CG': {'dH': 2.167, 'dS': 6.467},
    '+T+T/AG': {'dH': 4.940, 'dS': 12.895},

    # Patricia M. McTigue, Raymond J. Peterson, and Jason D. Kahn
    # Biochemistry 2004 43 (18), 5388-5405
    # DOI: 10.1021/bi035976d

    # --- 5' DNA to 3' LNA (Gap to Wing) ---
    # McTigue Notation: 5' AAL 3' (DNA A, LNA A)
    'A+A/TT': {'dH': 0.992, 'dS': 4.065},
    'A+C/TG': {'dH': 2.890, 'dS': 10.576},
    'A+G/TC': {'dH': -1.200, 'dS': -1.826},
    'A+T/TA': {'dH': 1.816, 'dS': 6.863},

    # McTigue Notation: 5' CAL 3' (DNA C, LNA A)
    'C+A/GT': {'dH': 1.358, 'dS': 4.367},
    'C+C/GG': {'dH': 2.063, 'dS': 7.565},
    'C+G/GC': {'dH': -0.276, 'dS': -0.718},
    'C+T/GA': {'dH': -1.671, 'dS': -4.070},

    # McTigue Notation: 5' GAL 3' (DNA G, LNA A)
    'G+A/CT': {'dH': 0.444, 'dS': 2.898},
    'G+C/CG': {'dH': -0.925, 'dS': -1.111},
    'G+G/CC': {'dH': -0.943, 'dS': -0.933},
    'G+T/CA': {'dH': -0.635, 'dS': -0.342},

    # McTigue Notation: 5' TAL 3' (DNA T, LNA A)
    'T+A/AT': {'dH': 1.591, 'dS': 5.281},
    'T+C/AG': {'dH': 0.609, 'dS': 3.169},
    'T+G/AC': {'dH': 2.165, 'dS': 7.163},
    'T+T/AA': {'dH': 2.326, 'dS': 8.051},

    # --- 5' LNA to 3' DNA (Wing to Gap) ---
    # McTigue Notation: 5' ALA 3' (LNA A, DNA A)
    '+AA/TT': {'dH': 0.707, 'dS': 2.477},
    '+AC/TG': {'dH': 1.131, 'dS': 4.064},
    '+AG/TC': {'dH': 0.264, 'dS': 2.613},
    '+AT/TA': {'dH': 2.282, 'dS': 7.457},

    # McTigue Notation: 5' CLA 3' (LNA C, DNA A)
    '+CA/GT': {'dH': 1.049, 'dS': 4.320},
    '+CC/GG': {'dH': 2.096, 'dS': 7.996},
    '+CG/GC': {'dH': 0.785, 'dS': 3.709},
    '+CT/GA': {'dH': 0.708, 'dS': 4.175},

    # McTigue Notation: 5' GLA 3' (LNA G, DNA A)
    '+GA/CT': {'dH': 3.162, 'dS': 10.544},
    '+GC/CG': {'dH': -0.360, 'dS': -0.251},
    '+GG/CC': {'dH': -2.844, 'dS': -6.680},
    '+GT/CA': {'dH': -0.212, 'dS': 0.073},

    # McTigue Notation: 5' TLA 3' (LNA T, DNA A)
    '+TA/AT': {'dH': -0.046, 'dS': 1.562},
    '+TC/AG': {'dH': 1.893, 'dS': 6.685},
    '+TG/AC': {'dH': -1.540, 'dS': -3.044},
    '+TT/AA': {'dH': 1.528, 'dS': 5.298},
}