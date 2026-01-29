def filter_gtf_genes(db, filter_mode):
    """
    Scans the GTF database and returns a set of gene names matching the criteria.
    """
    allowed_names = set()

    print(f"Querying GTF database for {filter_mode} genes...")

    # Iterate over all features of type 'gene'
    # Note: adjust 'gene' to 'transcript' if your GTF is transcript-level only,
    # but usually 'gene' is the standard feature type for this.
    for feature in db.features_of_type('gene'):

        # 1. Get Attributes safely (handles Gencode vs Ensembl naming)
        # 'gene_name' is standard, fallback to 'name'
        gene_name = feature.attributes.get('gene_name', feature.attributes.get('name', [None]))[0]

        if not gene_name:
            continue

        # 'gene_type' (Gencode) or 'gene_biotype' (Ensembl)
        biotype = feature.attributes.get('gene_type', feature.attributes.get('gene_biotype', [None]))[0]

        # 2. Check Mitochondrial Status
        # Standardize chrom names: 'chrM', 'MT', 'M'
        is_mitochondrial = feature.seqid in ('chrM', 'MT', 'M')

        # 3. Apply Filter Logic
        if filter_mode == 'protein_coding':
            # Strict: Must be Protein Coding AND Nuclear (usually implicit in "protein coding" filters for tools like this)
            if biotype == 'protein_coding' and not is_mitochondrial:
                allowed_names.add(gene_name)

        elif filter_mode == 'non_mt':
            # Loose: Everything except Mitochondria
            if not is_mitochondrial:
                allowed_names.add(gene_name)

        else:
            raise ValueError(f"Unknown filter mode: {filter_mode}")

    return allowed_names
