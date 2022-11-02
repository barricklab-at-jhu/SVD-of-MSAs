import sys
import Bio
from Bio import SeqIO

file = sys.argv[1]

fasta_sequences = SeqIO.parse(open(file),'fasta')

names = []
sequences = []

for fasta in fasta_sequences:
    name, seq = fasta.id, str(fasta.seq)
    names.append(name)
    seq = seq.replace('X', '-')
    sequences.append(seq)

seqs_gap_removed = [seq.replace('-', '').replace('.', '') for seq in sequences]

with open(f'{file[:-4]}_gapsRemoved.txt', 'w') as f:
    for name, seq in zip(names, seqs_gap_removed):
        f.write(f'>{name}\n{seq}\n')
