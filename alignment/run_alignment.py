from pathlib import Path

from fasta_utils import read_fasta, print_alignment
from swa import smith_waterman

print('Importing BRCA1 sequences...')
data_fp = Path(__file__).parent.parent/"data"

sequences = read_fasta(data_fp / "BRCA1.fasta")

print('Loading two genes...')
X = sequences['>XM_058182896.1 PREDICTED: Ahaetulla prasina BRCA1 DNA repair associated (BRCA1), transcript variant X1, mRNA']
Y = sequences['>XM_026071272.1 PREDICTED: Apteryx rowi BRCA1, DNA repair associated (BRCA1), transcript variant X1, mRNA']

print('Calculating local similarity...')
result = smith_waterman(X, Y)

print(f"Similarity: {result.similarity:.2f}")
print(print_alignment(result.seq_x_aligned, result.seq_y_aligned))