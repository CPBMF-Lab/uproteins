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

import pandas as pd
from .orflib import ORF, ORFCollection


class RefSeqGFF(object):
    def __init__(self, gff):
        self.gff = gff

    def __read_gff(self):
        df = pd.read_csv(self.gff, sep='\t', header=1)

        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

        df = df[df["source"].isin(["RefSeq", "ena"])]
        starts = df["start"].tolist()
        ends = df["end"].tolist()
        strands = df["strand"].tolist()
        attrs = df["attributes"].tolist()
        orf_list = []
        for i in range(len(attrs)):
            name = attrs[i].split(";")[1][12:]
            start = starts[i]
            end = ends[i]
            strand = strands[i]

            orf = ORF(name=name, start=start, end=end, strand=strand)
            orf_list.append(orf)
        return ORFCollection().add_orfs(orf_list)

    def get_dict(self):
        orfs = self.__read_gff()
        orf_dict = {}
        for orf in orfs:
            orf_dict[orf.name] = orf
        return orf_dict


class StringTieGFF(object):
    def __init__(self, gff):
        self.gff = gff
        self._fix_header()

    def _fix_header(self):
        with open(self.gff, 'r') as problemo, open('scapegoat.txt', 'w') as fixed_gff:
            fixed = []
            lines = problemo.readlines()
            for line in lines:
                if '# StringTie version' in line:
                    line = 'seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n'
                fixed.append(line)
            fixed_gff.writelines(fixed)
        os.system(f'mv {self.gff} gff_unfixed.gff')
        os.system(f'mv scapegoat.txt {self.gff}')

    def __read_gff(self):
        df = pd.read_csv(self.gff, sep='\t', header=1)
        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        df = df[df["feature"] == "transcript"]
        starts = df["start"].tolist()
        ends = df["end"].tolist()
        strands = df["strand"].tolist()
        attrs = df["attributes"].tolist()

        orf_list = []

        for i in range(len(attrs)):
            name = attrs[i].split(";")[1].split(" ")[-1].replace("\"", "")
            if 'uproteins' in name:
                name = '.'.join(name.split(".")[:2])
            # if 'ref_gene_id' in attrs[i]:
            #     name = attrs[i].split(";")[1][21:-1]
            #
            #     if 'rna-' in attrs[i]:
            #         name = attrs[i].split(";")[1][20:-1]

            start = starts[i]
            end = ends[i]
            strand = strands[i]
            orf = ORF(name=name, start=start, end=end, strand=strand)
            orf_list.append(orf)
        return ORFCollection().add_orfs(orf_list)

    def get_dict(self):
        orfs = self.__read_gff()
        orf_dict = {}
        for orf in orfs:
            orf_dict[orf.name] = orf
        return orf_dict


class GenomeCoordinatesRNA(object):
    """ For Transcriptome """
    def __init__(self, psm_table, string_dict, ref_dict=None):
        """ orf_dict is returned by get_dict() method from StringTieGFF class."""
        self.psm = pd.read_csv(psm_table, sep='\t')
        self.psm = self.psm[self.psm["proteinIds"].str.contains('ORF')]
        self.psm = self.psm[self.psm["proteinIds"].str.contains('|', regex=False) == False]
        self.ids = self.psm["proteinIds"].tolist()
        self.stringTieDict = string_dict
        self.refSeqDict = ref_dict
        self.coordinates = []
        # print(self.stringTieDict)

    def get_proteins(self):
        coordinates = []
        for i in self.ids:
            proteins = i.split(",")
            coord_set = ""
            for protein in proteins:
                # only for proteins in the transcriptome gff file
                if 'ORF' in protein:
                    if not ('gene' in protein or 'rna' in protein):
                        # pos1 = protein.find("_")
                        pos1 = protein.find("_") + 1
                        # pos2 = protein.find(".")
                        name = '.'.join(protein[pos1:].split(".")[:2])
                        # pos2 = protein.find("_", pos1)
                        # name = protein[pos1:pos2]
                    else:
                        # pos1 = protein.find("_gene-") + 6
                        pos1 = protein.find("_") + 1
                        # pos2 = protein.rfind("_", pos1, -1)
                        pos2 = protein.find("_", pos1)
                        # pos3 = protein.rfind("_", pos1, pos2)
                        # name = protein[pos1:pos3]
                        name = protein[pos1:pos2]
                    if name in self.stringTieDict:
                        start = self.stringTieDict[name].start
                        end = self.stringTieDict[name].end
                        orf_pos = protein.rfind("_", 0, -1) + 1
                        orf_coords = protein[orf_pos:-1].split("-")
                        orf_start = orf_coords[0]
                        orf_end = orf_coords [1]

                        genome_start = int(start) + int(orf_start)
                        genome_end = int(start) + int(orf_end)

                        if genome_start > genome_end:
                            start = genome_end
                            end = genome_start
                        else:
                            start = genome_start
                            end = genome_end
                        coords = f'{start}-{end}'
                        # print(coords)

                        if len(coord_set) > 0:
                            coord_set += f',{coords}'
                        else:
                            coord_set += f'{coords}'
                    else:
                        if len(coord_set) > 0:
                            coord_set += ',not found'
                        else:
                            coord_set += 'not found'
                # only for already annotated proteins not found in StringTie GFF file
                else:
                    if protein in self.refSeqDict:
                        orf = self.refSeqDict[protein]
                        start = orf.start - 1
                        end = orf.end
                        coords = f'{start}-{end}'
                        # print(orf, start, end, coords)
                        if len(coord_set) > 0:
                            coord_set += f',{coords}'
                        else:
                            coord_set += f'{coords}'
                    else:
                        if len(coord_set) > 0:
                            coord_set += ',not found'
                        else:
                            coord_set += 'not found'
            # print(coord_set)
            self.coordinates.append(coord_set)

        print(self.coordinates)
        return self

    def save_table(self, output):
        df = self.psm
        df.insert(5, "Genome Coordinates", self.coordinates)
        df.to_csv(f'{output}.txt', sep="\t", index=False)
        return self


def findnth(string, substring, n):
    parts = string.split(substring, n + 1)
    if len(parts) <= n + 1:
        return -1
    return len(string) - len(parts[-1]) - len(substring)


class GenomeCoordinates(object):
    """ For genome ORFs """
    def __init__(self, psm_table, ref_dict=None):
        """ ref_dict is returned by get_dict() method from RefSeqGFF class. """
        self.psm = pd.read_csv(psm_table, sep='\t')
        self.psm = self.psm[self.psm["proteinIds"].str.contains('|', regex=False) == False]
        self.ids = self.psm["proteinIds"].tolist()
        self.refSeqDict = ref_dict
        self.coordinates = []

    def get_coords(self):
        for i in range(len(self.ids)):
            proteins = self.ids[i].split(",")
            coord_set = ""
            for protein in proteins:
                if 'gORF' in protein:
                    # pos1 = protein.find("_") + 1
                    # print(protein)
                    pos1 = findnth(protein, '_', 2) + 1
                    pos2 = protein.rfind("_")
                    coords = protein[pos1:pos2].split("-")
                    start_c = coords[0]
                    end_c = coords[1]
                    if int(start_c) > int(end_c):
                        start = end_c
                        end = start_c
                    else:
                        start = start_c
                        end = end_c
                    coords = f'{start}-{end}'
                    if len(coord_set) > 0:
                        coord_set += f',{coords}'
                    else:
                        coord_set += coords
                else:
                    if protein in self.refSeqDict:
                        orf = self.refSeqDict[protein]
                        start = orf.start
                        end = orf.end
                        coords = f'{start}-{end}'
                        if len(coord_set) > 0:
                            coord_set += f',{coords}'
                        else:
                            coord_set += coords
                    else:
                        if len(coord_set) > 0:
                            coord_set += ',not found'
                        else:
                            coord_set += 'not found'
            self.coordinates.append(coord_set)
        return self

    def save_table(self, output):
        df = self.psm
        df.insert(5, "Genome Coordinates", self.coordinates)
        df.to_csv(f'{output}.txt', sep="\t", index=False)
        return self
