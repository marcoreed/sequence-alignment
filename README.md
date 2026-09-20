# Medical Sequence Alignment 

A project exploring sequence alignment algorithms and their applications spanning classic DNA sequence alignment to hospital event sequence alignment. 

## Algorithms

[Smith Waterman Algorithm](docs/smith_waterman.md)

## Installation

```bash
git clone https://github.com/marcoreed/sequence-alignment
cd sequence-alignment
python -m venv swa_env
source swa_env/bin/activate
pip install -r requirements.txt
```


## Usage

Run pairwise alignment between two sequences from the FASTA file provided:

```bash
python run_alignment.py
```

Edit the sequence IDs at the top of `run_alignment.py` to select which genes to align.

## References

Smith TF, Waterman MS. Identification of common molecular subsequences. J Mol Biol. 1981 Mar 25;147(1):195-7. doi: 10.1016/0022-2836(81)90087-5. PMID: 7265238.


