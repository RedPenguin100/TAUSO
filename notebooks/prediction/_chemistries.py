"""Data-validated ASO chemistries for the design scripts.

Each preset is a gapmer that is well represented in the training data, so the model interpolates
rather than extrapolates. Arbitrary backbones (partial PS) and off-grid wing/gap geometries are
out of distribution -- the model does not rank them reliably -- and are intentionally not offered.

Each entry: (aso_size, chemical_pattern, modification, output_subdir). Every preset uses a full
phosphorothioate backbone, and writes to its own output subfolder so chemistries stay separated.
"""

CHEMISTRIES = {
    # 5-10-5 2'-MOE 20-mer gapmer -- 94% of MOE gapmers in the data (the default "vanilla").
    "moe_5_10_5": (20, "MMMMMddddddddddMMMMM", "2'MOE/5-methylcytosines/deoxy", "2moe"),
    # 3-10-3 cEt 16-mer gapmer -- 98.5% of cEt gapmers in the data.
    "cet_3_10_3": (16, "CCCddddddddddCCC", "cEt/5-methylcytosines/deoxy", "cet"),
}


def build_config(default_config, chemistry):
    """Return (config, aso_size, output_subdir) for the named chemistry."""
    aso_size, chemical_pattern, modification, subdir = CHEMISTRIES[chemistry]
    config = default_config()
    config.standard_chemical_pattern = chemical_pattern
    config.standard_modification = modification
    config.standard_ps_pattern = "*" * (aso_size - 1)
    return config, aso_size, subdir
