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


import pathlib
from importlib import resources as rsrc

import pandas as pd
import pytest

from tests import resources
from src import uproteins


pytestmark = pytest.mark.ms


def test_ms_mode(database_folder):
    mzml = rsrc.files(resources).joinpath("mzml")

    args = [
        'ms',
        '--outdir', str(database_folder),
        '--Mass_spec', str(mzml),
        '--inst', '2',
        '--t', '800ppm',
        '--Transcriptome', 'YES'
    ]
    uproteins(args)

    genome_pin_path: pathlib.Path = (
        database_folder
        / 'Genome'
        / 'Percolator'
        / 'Genome_pin.txt'
    )
    transcriptome_pin_path: pathlib.Path = (
        database_folder
        / 'Transcriptome'
        / 'Percolator'
        / 'Transcriptome_pin.txt'
    )

    assert genome_pin_path.is_file()
    assert transcriptome_pin_path.is_file()

    with rsrc.path(resources, 'ms_results') as results:
        ok_genome_pin_path = (
            results
            / 'Genome'
            / 'Percolator'
            / 'Genome_pin.txt'
        )
        ok_transcriptome_pin_path = (
            results
            / 'Transcriptome'
            / 'Percolator'
            / 'Transcriptome_pin.txt'
        )

        genome_pin = pd.read_csv(genome_pin_path, comment='#', sep='\t')
        ok_genome_pin = pd.read_csv(ok_genome_pin_path, comment='#', sep='\t')

        transcriptome_pin = pd.read_csv(
            transcriptome_pin_path,
            comment='#',
            sep='\t'
        )
        ok_transcriptome_pin = pd.read_csv(
            ok_transcriptome_pin_path,
            comment='#',
            sep='\t'
        )

        pd.testing.assert_frame_equal(
            genome_pin,
            ok_genome_pin,
            check_like=True
        )
        pd.testing.assert_frame_equal(
            transcriptome_pin,
            ok_transcriptome_pin,
            check_like=True
        )
