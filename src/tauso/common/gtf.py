import logging

logger = logging.getLogger(__name__)

MITO_SEQIDS = ("chrM", "MT", "M")


def _first_attr(feature, keys):
    """Return the first present, non-empty attribute value from keys."""
    for key in keys:
        val = feature.attributes.get(key, [])
        if val:
            return val[0]
    return None


def _filter_genes(db, filter_mode, name_keys, biotype_keys, source_label):
    if filter_mode not in ("protein_coding", "non_mt"):
        raise ValueError(f"Unknown filter mode: {filter_mode}")

    allowed_names = set()
    logger.info(f"Querying {source_label} database for {filter_mode} genes...")

    for feature in db.features_of_type("gene"):
        gene_name = _first_attr(feature, name_keys)
        if not gene_name:
            continue

        is_mitochondrial = feature.seqid in MITO_SEQIDS
        if is_mitochondrial:
            continue

        if filter_mode == "protein_coding":
            biotype = _first_attr(feature, biotype_keys)
            if biotype != "protein_coding":
                continue

        allowed_names.add(gene_name)

    return allowed_names


def filter_gtf_genes(db, filter_mode):
    """Scan a GTF database and return a set of gene names matching the criteria."""
    return _filter_genes(
        db, filter_mode,
        name_keys=("gene_name", "name"),
        biotype_keys=("gene_type", "gene_biotype"),
        source_label="GTF",
    )


def filter_gff_genes(db, filter_mode):
    """Scan a GFF3 database and return a set of gene names matching the criteria."""
    return _filter_genes(
        db, filter_mode,
        name_keys=("gene_name", "Name", "gene"),
        biotype_keys=("gene_type", "biotype", "gene_biotype"),
        source_label="GFF3",
    )