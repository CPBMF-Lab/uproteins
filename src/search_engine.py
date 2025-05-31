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


import argparse
import subprocess
import sys
import shutil
import pathlib as p
import typing as t

import pandas as pd


class SearchEngine(t.Protocol):
    def __init__(self, args: argparse.Namespace):
        """Initialize the search engine with configurations.

        Arguments
        ---------
        folder : str
            Either the Genome or Transcriptome folder.
        args : Namespace
            An :obj:`argparse.Namespace` object used to get the SearchEngine's
            run configs from the cli arguments.
        """
        ...

    def run(
        self,
        folder: str,
        mzml_path: p.Path,
        database: p.Path
    ) -> subprocess.CompletedProcess[str]:
        """Execute the search and return results.

        Arguments
        ---------
        mzml_folder : Path
            The path to a folder containing the mzML files to execute the
            search.
        database : Path
            A path to the database fasta file to search against.

        Returns
        -------
        CompletedProcess
            :obj:`CompletedProcess` instance returned by the subprocess call.
        """
        ...

    def save_to_pin(
        self,
        folder: str,
        database: p.Path,
        decoy: p.Path
    ) -> None:
        """Convert search results to pin format.

        Arguments
        ---------
        result_files : list[Path]
            A list of search result files that should be converted to a single
            percolator pin file.
        database : Path
            The database file used for the search.
        decoy : Path
            The decoy database file used for the search.
        """
        ...


class CometMS:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.params = self._get_params()

    def run(
        self,
        folder: str,
        mzml_path: p.Path,
        database: p.Path
    ) -> subprocess.CompletedProcess[str]:
        self.params['database_name'] = str(database.absolute())
        comet_path = (
            f'{sys.path[0]}/dependencies/'
            f'CometMS/{folder}/comet.linux.exe'
        )
        self._build_params_file(folder)
        input_files = f'{mzml_path}/*.mzML'
        cmd = [comet_path, input_files]
        result = subprocess.run(cmd, text=True, capture_output=True)
        for file in mzml_path.iterdir():
            if file.suffix == '.pin':
                shutil.copy(file, folder)
        return result

    def save_to_pin(
        self,
        folder: str,
        database: p.Path,
        decoy: p.Path
    ) -> None:
        # We need database and decoy for compliance with the SearchEngine
        # Protocol
        pin_list = [
            pd.read_csv(file, sep='\t', header=1)
            for file in p.Path(f'{folder}/Percolator').iterdir()
            if file.name.endswith(f'_{folder}.pin')
        ]
        pd.concat(pin_list).to_csv(
            f'{folder}/Percolator/{folder}_pin.txt',
            sep='\t'
        )

    def _get_params(self) -> dict[str, str]:
        pipeline_args = (
            'outdir', 'mode', 'mass_spec', 'processes', 'transcriptome', 'xmx',
            'threads'
        )
        param_dict = {
            key: value
            for key, value
            in vars(self.args).items()
            if key not in pipeline_args
        }
        param_dict['output_pepxmlfile'] = '0'
        param_dict['output_percolatorfile'] = '1'
        if self.args.threads is not None:
            param_dict['num_threads'] = str(self.args.threads)
        return param_dict

    def _build_params_file(self, folder):
        with open(
            f'{sys.path[0]}/dependencies/'
            f'CometMS/{folder}/comet.params',
            'w'
        ) as params:
            for key, value in self.params.items():
                params.write(f'{key} = {value}')


class MSGFPlus:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.base_command = self._build_command()

    def run(
        self,
        folder: str,
        mzml_path: p.Path,
        database: p.Path
    ) -> subprocess.CompletedProcess[str]:
        output = p.Path(
            f'{folder}/{database.with_suffix(".mzid").name}'
        )
        cmd = [
            *self.base_command,
            '-d', str(database),
            '-s', str(mzml_path),
            '-o', str(output)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def save_to_pin(
        self,
        folder: str,
        database: p.Path,
        decoy: p.Path
    ) -> None:
        target, decoy = self._generate_metafiles(folder)
        msgf2pin_path = p.Path('msgf2pin')
        cmd_pin = [
            str(msgf2pin_path),
            str(target), str(decoy),
            '-o', f'{folder}/Percolator/{folder}_pin.txt',
            '-F', str(database) + ',' + str(decoy),
            '-c', '2'
        ]
        subprocess.run(cmd_pin, text=True, check=True)

    def _generate_metafiles(self, folder: str) -> tuple[p.Path, p.Path]:
        folder_path = p.Path(folder)

        decoy_meta = folder_path / f'{folder}_decoy_metafile.txt'
        target_meta = folder_path / f'{folder}_target_metafile.txt'

        decoys = [
            str(file.absolute())
            for file in folder_path.iterdir()
            if file.suffix == '.mzid'
            if 'decoy' in file.name
        ]
        targets = [
            str(file.absolute())
            for file in folder_path.iterdir()
            if file.suffix == '.mzid'
            if 'decoy' not in file.name
        ]

        with decoy_meta.open('w') as decoy:
            decoy.writelines('\n'.join(decoys))
        with target_meta.open('w') as target:
            target.writelines('\n'.join(targets))

        return target_meta, decoy_meta

    def _build_command(self) -> list[str]:
        cmd = [
            'java',
            f'-Xmx{self.args.xmx}',
            '-jar', f'{sys.path[0]}/dependencies/MSGF/MSGFPlus.jar',
            '-tda', '0',
            '-addFeatures', '1'
        ]
        pipeline_args = (
            'outdir', 'mode', 'mass_spec', 'processes', 'transcriptome', 'xmx',
            'threads'
        )
        for key, value in vars(self.args).items():
            if key not in pipeline_args and value is not None:
                cmd.extend([f'-{key}', value])
        if self.args.threads is not None:
            cmd.extend(['-thread', str(self.args.threads)])
        return cmd
