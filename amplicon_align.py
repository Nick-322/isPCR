#!/usr/bin/env python3
from magnumopus import ispcr, needleman_wunsch
import argparse


parser = argparse.ArgumentParser(
    description="Perform in-silico PCR on two assemblies and align the amplicons",
    usage="synopsis: amplicon_align.py -1 ASSEMBLY1 -2 ASSEMBLY2 -p PRIMERS -m MAX_AMPLICON_SIZE --match MATCH --mismatch MISMATCH --gap GAP"
)
parser.add_argument('-1', '--ASSEMBLY1', required=True, help="Path to the first assembly file")
parser.add_argument('-2', '--ASSEMBLY2', required=True, help="Path to the second assembly file")
parser.add_argument('-p', '--PRIMERS', required=True, help="Path to the primer file")
parser.add_argument('-m', '--MAX_AMPLICON_SIZE', required=True, type=int,
                    help="Maximum amplicon size for isPCR")
parser.add_argument('--match', required=True, type=int, help="Match score to use in alignment")
parser.add_argument('--mismatch', required=True, type=int, help="Mismatch penalty to use in alignment")
parser.add_argument('--gap', required=True, type=int, help="Gap penalty to use in alignment")

args = parser.parse_args()

assembly1 = args.ASSEMBLY1
assembly2 = args.ASSEMBLY2
primers = args.PRIMERS
max_amplicon_size = args.MAX_AMPLICON_SIZE
match = args.match
mismatch = args.mismatch
gap = args.gap


amplicon1 = ispcr(primers, assembly1, max_amplicon_size)
lines = amplicon1.split('\n')
sequences = [line for line in lines if not line.startswith('>')]
sequence_1 = '\n'.join(sequences)
sequence_1 = sequence_1.replace('\n', '')


amplicon2 = ispcr(primers, assembly2, max_amplicon_size)
lines2 = amplicon2.split('\n')
sequences2 = [line for line in lines2 if not line.startswith('>')]
sequence_2 = '\n'.join(sequences2)
sequence_2 = sequence_2.replace('\n', '')

complement_table = {"A":"T", "T":"A", "G":"C", "C":"G"}
reverse_seq1 = sequence_1[::-1]
reverse_complement1 = "".join(complement_table[char] for char in reverse_seq1)

reverse_seq2 = sequence_2[::-1]
reverse_complement2 = "".join(complement_table[char] for char in reverse_seq2)

list1 = []
list2 = []
list1.append(sequence_1)
list1.append(reverse_complement1)
list1.append(sequence_1)

list2.append(sequence_2)
list2.append(sequence_2)
list2.append(reverse_complement2)

best_alignment = None
max_score = 0

#Check alignment and its score on all possible sequence combinations
#And output the best alignment and its score
for item1, item2 in zip(list1, list2):
    aln, score = needleman_wunsch(item1, item2, match, mismatch, gap)
    if score > max_score:
        max_score = score
        best_alignment = aln
print("\n".join(best_alignment) + "\n")
print(max_score)




