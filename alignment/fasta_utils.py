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

def print_alignment(seq_a: str, seq_b: str, width: int = 60) -> str:
    middle = ''
    for a, b in zip(seq_a, seq_b):
        if a == b:
            middle += '|'
        elif a == '-' or b == '-':
            middle += ' '
        else:
            middle += '.'

    lines = []

    for i in range (0, len(seq_a), width):
        lines.append(seq_a[i:i+width])
        lines.append(middle[i:i+width])
        lines.append(seq_b[i:i+width])
        lines.append('')

    return '\n'.join(lines)

if __name__ == "__main__":
    test = read_fasta(data_dir/"BRCA1.fasta")
    print(next(iter(test.items())))

    alignment_test = print_alignment('AAA-AGA','AAATATA')
    print(alignment_test)

