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


import os
import sys
from Bio import SeqIO


class Decoy(object):
    def __init__(self, db, db_type):
        self.df = db
        self.seqs = self.__get_seqs()
        self.reversed = []
        self.type = db_type
        self.__create_dir()
        self.path = sys.path[0]

    def __create_dir(self):
        if not os.path.exists(f"{self.type}/Percolator"):
            cmd_dir = f'mkdir {self.type}/Percolator'
            os.system(cmd_dir)

    def __get_seqs(self):
        seqs = []
        records = SeqIO.parse(self.df, 'fasta')
        for record in records:
            seqs.append(record.seq)
        return seqs

    def reverse_sequences(self):
        """ Reverses the amino acid sequence of a protein, except for the aa in the c-terminal. """
        for seq in self.seqs:
            # cterminus = seq[len(seq)-1]
            to_reverse = seq[:-1]
            reversed = to_reverse[::-1]
            # reversed += cterminus
            self.reversed.append(reversed)
        return self

    def add_contaminants(self):
        seqs = []
        records = SeqIO.parse(f'{self.path}/seqlib/contaminants.txt', 'fasta')
        for record in records:
            seqs.append(f">{record.id}\n{record.seq}\n")
        return seqs

    def to_fasta(self):
        out = []
        for i in range(len(self.reversed)):
            string = f">{self.type}_decoy_Protein_{i+1}\n{self.reversed[i]}\n"
            out.append(string)
        seqs = self.add_contaminants()
        for seq in seqs:
            to_add = seq.replace("B", "")
            to_add = to_add.replace("X", "")
            to_add = to_add.replace("Z", "")
            out.append(to_add)
        with open(f"{self.type}/Percolator/{self.type}_decoy.fasta", 'w') as fa:
            fa.writelines(out)
        return self
        
#genome_decoy = Decoy(db="genome_database.fasta", db_type="Genome")
#genome_decoy.reverse_sequences().to_fasta()
# transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
# transcriptome_decoy.reverse_sequences().to_fasta()
