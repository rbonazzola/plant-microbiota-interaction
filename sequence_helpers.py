from Bio import SeqIO

def get_contig_range(input_fasta, contig_name, start=None, end=None):
    """
    Retrieves a specific contig (and optionally a range within it) as a string.

    Parameters:
    - input_fasta: str, path to the input FASTA file.
    - contig_name: str, the name of the contig to query.
    - start: int or None, start position (1-based, inclusive). Default is None for full contig.
    - end: int or None, end position (1-based, inclusive). Default is None for full contig.

    Returns:
    - str: The extracted sequence as a string.
    - None: If the contig is not found.
    """
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id == contig_name:
            if start is not None and end is not None:
                return str(record.seq[start - 1:end])
            return str(record.seq)
    print(f"Contig '{contig_name}' not found in {input_fasta}")
    return None


# Example usage
# input_file = "data/genomes/reference_2/full_sequence.fasta"
# contig = "PRKV01000004.1"
# start_position = 394240 - 200
# end_position = 394240 + 200

# get_contig_range(input_file, contig, start_position, end_position)