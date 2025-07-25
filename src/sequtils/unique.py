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


import pandas as pd

from .postsearch import Peptide


class PercolatorUTP(object):
    def __init__(self, coord_df=None, pep=0.01, qvalue=0.05):
        """ 'coord_df' must be a data frame created by get_coords() method from either GenomeCoordinates or
        GenomeCoordinatesRNA classes. """
        self.df = pd.read_csv(coord_df, sep="\t")
        self.df = self.df[self.df["Genome Coordinates"].str.contains('not found') == False]
        self.df = self.df[self.df["proteinIds"].str.contains('ORF')]
        self.df = self.df[self.df["Genome Coordinates"] != "Genome Coordinates"]
        # self.df = self.df[self.df["q-value"]]
        # self.df = self.df[""]
        self.df = self.df[(self.df["q-value"] < qvalue) & (self.df["posterior_error_prob"] < pep)]
        # self.df = self.df[self.df["posterior_error_prob"] <= 0.01]
        # self.df = self.df[self.df["proteinIds"].str.contains('lcl|') == False]
        self.df = self.df[self.df["proteinIds"].str.contains("|", regex=False) == False]
        self.df = self.df[self.df["proteinIds"].str.contains('Decoy', regex=False) == False]
        self.coordinates = self.df["Genome Coordinates"].tolist()
        self.peptides = self.df["peptide"].tolist()
        self.unique = []

    def get_utps(self):
        for pepseq, coords in zip(self.peptides, self.coordinates):
            coordinates = coords.split(",")
            peptide = Peptide(pepseq)
            if len(coordinates) < 10:
                for i in coordinates:
                    start = i.split("-")[0]
                    end = i.split("-")[1]
                    peptide.add_spec(start, end)
                    peptide.check_loci()
                # self.unique.append(peptide.unique)
                self.unique.append(True)

            else:
                self.unique.append(False)
        return self

    def save_table(self, output='unique_peptides', keep='all'):
        self.df.insert(6, "Unique Peptide", self.unique)
        if keep == 'unique':
            self.df = self.df[self.df["Unique Peptide"] == True]
        self.df.to_csv(f'{output}.txt', sep="\t", index=False)
        return self




