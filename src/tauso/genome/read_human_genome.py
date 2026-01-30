import bisect

from Bio.Seq import Seq

from ..data.data import load_db, load_genome
from .LocusInfo import LocusInfo
from ..timer import Timer


def cond_print(text, verbose=False):
    if verbose:
        print(text)


def get_locus_to_data_dict(include_introns=True, gene_subset=None, genome='GRCh38'):
    with Timer() as t:
        db = load_db(genome)
    print("Elapsed DB: ", t.elapsed_time)
    fasta_dict = load_genome(genome)
    print("Elapsed Fasta: ", t.elapsed_time)

    print("Length: ", len(fasta_dict))

    locus_to_data = dict()
    locus_to_strand = dict()

    basic_features = ['exon', 'gene', 'stop_codon', 'UTR']

    feature_types = list(basic_features)

    if include_introns:
        feature_types.append('intron')


    # --- OPTIMIZATION START ---
    iterator = []

    # If we have specific genes, don't scan the whole genome.
    if gene_subset:
        target_names = set(gene_subset)

        # 2. Iterate ONLY genes (Fast, ~60k items vs millions)
        # Note: If you have Gene IDs (ENSG...), you could do db[id] which is O(1).
        # Since you use names (KLKB1), we scan the genes.
        for gene in db.features_of_type('gene'):
            name_list = gene.attributes.get('gene_name', [])
            if name_list and name_list[0] in target_names:
                # Add the gene itself
                iterator.append(gene)

                # Add its children (exons, introns, etc.)
                # db.children() uses the DB index = INSTANT
                children = db.children(gene, featuretype=feature_types, order_by='start')
                iterator.extend(list(children))
    else:
        # Full genome scan
        iterator = db.features_of_type(feature_types, order_by='start')
    # --- OPTIMIZATION END ---

    for feature in iterator:
        chrom = feature.seqid
        if 'chrM' == chrom:
            continue
        locus_tags = feature.attributes['gene_name']
        if len(locus_tags) != 1:
            raise ValueError(f"Multiple loci: {locus_tags}")
        locus_tag = locus_tags[0]

        if gene_subset is not None:
            if locus_tag not in gene_subset:
                continue

        if locus_tag not in locus_to_data:
            locus_info = LocusInfo()
            locus_to_data[locus_tag] = locus_info
        else:
            locus_info = locus_to_data[locus_tag]

        if feature.featuretype == 'exon':
            exon = feature
            # seq = fasta_dict[chrom].seq[exon.start - 1: exon.end]
            seq = Seq(str(fasta_dict[chrom][exon.start - 1: exon.end]))
            if exon.strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()

            bisect.insort(locus_info.exons, (exon.start - 1, seq))
            bisect.insort(locus_info.exon_indices, (exon.start - 1, exon.end))
            locus_to_strand[locus_tag] = exon.strand

        elif feature.featuretype == 'intron' and include_introns:
            intron = feature
            # seq = fasta_dict[chrom].seq[intron.start - 1: intron.end]
            seq = Seq(str(fasta_dict[chrom][intron.start - 1: intron.end]))

            if intron.strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()

            bisect.insort(locus_info.introns, (intron.start - 1, seq))
            bisect.insort(locus_info.intron_indices, (intron.start - 1, intron.end))
            locus_to_strand[locus_tag] = intron.strand

        elif feature.featuretype == 'gene':
            gene = feature
            # seq = fasta_dict[chrom].seq[gene.start - 1: gene.end]
            seq = Seq(str(fasta_dict[chrom][gene.start - 1: gene.end]))

            if gene.strand == '-':
                seq = seq.reverse_complement()
            seq = str(seq).upper()

            locus_info.strand = gene.strand
            locus_info.gene_start = gene.start - 1
            locus_info.gene_end = gene.end
            locus_info.full_mrna = seq
            locus_to_strand[locus_tag] = gene.strand


        elif 'UTR' in feature.featuretype:
            utr = feature
            bisect.insort(locus_info.utr_indices, (utr.start - 1, utr.end))
        elif feature.featuretype == 'stop_codon':
            locus_info.stop_codons.append((feature.start, feature.end))
        else:
            print("Feature type: ", feature.featuretype)


    for locus_tag in locus_to_data:
        locus_info = locus_to_data[locus_tag]
        if locus_to_strand[locus_tag] == '-':
            locus_info.exons.reverse()
            if include_introns:
                locus_info.introns.reverse()
        locus_info.exons = [element for _, element in locus_info.exons]

        if include_introns:
            locus_info.introns = [element for _, element in locus_info.introns]

    return locus_to_data
