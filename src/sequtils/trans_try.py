# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


import regex as re
from Bio import SeqIO

# seqs = []
# records = SeqIO.parse('mtb_genome.fasta', 'fasta')
# for record in records:
#     seqs.append(record.seq)
seqs = ["ATGTGTGACCCGTGTGACGACATGCTGATGGACCATGATGCTATAAGGCCATGAC"]

starts = ['ATG', 'TGT']
stops = ['GAC']

# stop_pattern = ""
# for i in range(len(stops)):
#     if i != 0:
#         stop_pattern += f"|(?={stops[i]})"
#     else:
#         stop_pattern += f"(?={stops[i]})"
# start_pattern = ""
#
# for i in range(len(starts)):
#     if i != 0:
#         start_pattern += f"|({starts[i]})"
#     else:
#         start_pattern += f"({starts[i]})"

def __get_pattern(codons, codon_type=None):
    """ codon_type is either 'start' or 'stop'. 'codons' must refer to a list of start or stop codons."""
    pattern = ""
    for i in range(len(codons)):
        if i != 0:
            if codon_type == "stop":
                pattern += f"|(?={codons[i]})"
            elif codon_type == "start":
                pattern += f"|({codons[i]})"
        else:
            if codon_type == "stop":
                pattern += f"(?={codons[i]})"
            elif codon_type == "start":
                pattern += f"({codons[i]})"
    return pattern

start_pattern = __get_pattern(starts, codon_type="start")
stop_pattern = __get_pattern(stops, codon_type="stop")

for seq in seqs:

    # aa = re.finditer('%s(...)+?((?=GAC)|(?=CCG))' % start, "%s"%seq, overlapped=True)
    aa = re.finditer('(%s)(...)+?(%s)' % (start_pattern, stop_pattern), "%s"%seq, overlapped=True)

    for a in aa:
        print(a.group(), a.start(), a.end())

