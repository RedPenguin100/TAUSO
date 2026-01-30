import requests

class TranscriptStructure:
    def __init__(self, transcript_id, full_mrna, exon_indices, intron_indices, utr_indices, cds_start):
        self.transcript_id = transcript_id
        self.full_mrna = full_mrna
        self.exon_indices = exon_indices
        self.intron_indices = intron_indices
        self.utr_indices = utr_indices
        self.cds_start = cds_start

    def __repr__(self):
        return (f"TranscriptStructure({self.transcript_id}, "
                f"{len(self.full_mrna)} nt, {len(self.exon_indices)} exons)")

# ------------------------------- #
# Main function to extract data
# ------------------------------- #
def get_transcript_structure(transcript_id: str):
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}

    # 1️⃣ Fetch transcript metadata (includes exons & UTRs)
    url = f"{server}/lookup/id/{transcript_id}?expand=1"
    r = requests.get(url, headers=headers)
    if not r.ok:
        raise RuntimeError(f"Failed to fetch transcript data: {r.text}")
    data = r.json()

    # 2️⃣ Fetch full pre-mRNA sequence
    seq_url = f"{server}/sequence/id/{transcript_id}?type=genomic"
    seq_r = requests.get(seq_url, headers=headers)
    if not seq_r.ok:
        raise RuntimeError(f"Failed to fetch sequence: {seq_r.text}")
    full_mrna = seq_r.json()["seq"]

    # 3️⃣ Extract exon coordinates relative to transcript
    exon_indices = []
    exons = sorted(data["Exon"], key=lambda e: e["start"])
    for e in exons:
        exon_indices.append((e["start"], e["end"]))

    # 4️⃣ Derive intron coordinates
    intron_indices = []
    for (prev_end, next_start) in zip([e["end"] for e in exons[:-1]],
                                      [e["start"] for e in exons[1:]]):
        intron_indices.append((prev_end + 1, next_start - 1))

    # 5️⃣ Extract CDS and UTR info
    cds_start = data.get("Translation", {}).get("start", None)
    utr_indices = []
    if "UTR" in data:
        utr_indices = [(u["start"], u["end"]) for u in data["UTR"]]

    # Create structure object
    return TranscriptStructure(
        transcript_id=transcript_id,
        full_mrna=full_mrna,
        exon_indices=exon_indices,
        intron_indices=intron_indices,
        utr_indices=utr_indices,
        cds_start=cds_start
    )


# ------------------------------- #
# Example for IDO1
# ------------------------------- #
if __name__ == "__main__":
    transcript_id = "ENST00000518237"
    IDO1_structure_data_object = get_transcript_structure(transcript_id)

    print("✅ Successfully retrieved transcript structure!")
    print(IDO1_structure_data_object)

    # You can now store it like:
    genes_u = ["IDO1"]
    gene_to_data = {"IDO1": IDO1_structure_data_object}


