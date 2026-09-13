from pathlib import Path

data_dir = Path(__file__).parent.parent/"Data"

def read_fasta(input_fp):
    record = {}
    with open(input_fp) as file:
        for line in file:
            if line.startswith('>'):
                header = line.strip()
                record[header] = ""
            else:
                record[header] += line.strip()

    return record

if __name__ == "__main__":
    test = read_fasta(data_dir/"BRCA1.fasta")
    print(next(iter(test.values())))

