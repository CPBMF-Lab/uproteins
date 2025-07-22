# Copyright © 2025 Eduardo Vieira de Souza
# Copyright © 2025 Adriana Canedo
# Copyright © 2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
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

from argparse import Namespace
from importlib import resources as rsrc
import shutil

import pytest

from src.postprocess.specfilt import PostPercolator
from src import cli
from tests import resources

class TestPostPercolator:
    def setup(self):
        args = Namespace(
            mode='postms',
            genome=rsrc.files(resources).joinpath("genome.fasta"),
            proteome=rsrc.files(resources).joinpath("proteome.fasta"),
            gtf=rsrc.files(resources).joinpath("mtb.gtf"),
            mzml=rsrc.files(resources).joinpath("mzml"),
            rrna=rsrc.files(resources).joinpath("rrna.fna"),
            transcriptome=True,
        )

        self.gPostPerc = PostPercolator(args, 'Genome', 'genome')
        self.tPostPerc = PostPercolator(args, 'Transcriptome', 'transcriptome')

    def test_convert_output(self, class_dir, monkeypatch):
        monkeypatch.chdir(class_dir)
        gene_pout = rsrc.files(resources).joinpath(
            'results/'
            'Genome/'
            'Percolator/'
            'genome_results_psm.txt'
        )
        transc_pout = rsrc.files(resources).joinpath(
            'results/'
            'Transcriptome/'
            'Percolator/'
            'transcriptome_results_psm.txt'
        )
        shutil.copytree(gene_pout, class_dir)  # pyright: ignore
        shutil.copytree(transc_pout, class_dir)  # pyright: ignore
        
