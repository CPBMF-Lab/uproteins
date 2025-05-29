# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
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


import os
import sys
import pathlib
import subprocess
import multiprocessing

from src import cli


class PeptideSearch(object):
    def __init__(
        self,
        database_type,
        ms_files_folder,
        database_file,
        decoy_file, args
    ):
        self.database_type = database_type
        self.ms_folder = pathlib.Path(ms_files_folder)
        self.database_file = pathlib.Path(database_file)
        self.decoy_file = pathlib.Path(decoy_file)
        self.args = args
        self.msgf_command = self._build_msgf_command()
        self._check_folder()
        self.path = sys.path[0]

    def peptide_identification(self):
        print(f"\nPerforming peptide search using {self.database_type} database\n")
        cmd_dir_ms = 'mkdir %s' % self.database_type
        os.system(cmd_dir_ms)
        cmd_copy_orf = 'cp %s %s/' % (self.database_file, self.database_type)
        os.system(cmd_copy_orf)
        cmd_copy_db = f'cp {self.database_file} {os.path.abspath(self.ms_folder)}/.'
        os.system(cmd_copy_db)

        self.parallel_search()
        cmd_move = 'mv %s/*.mzid %s/' % (self.ms_folder, self.database_type)
        os.system(cmd_move)

    def parallel_search(self):
        # A list of tasks to be run, where the first element represents whether
        # to use a decoy or not, and the second element is a strPath to a mzml
        # file.
        # This list has all possible combinations of decoy/not_decoy and each
        # mzml file. This allows both real and decoy runs to run in parallel
        # with the same pool.
        tasks = [
            (d, f) for d in (True, False)
            for f in self.ms_folder.iterdir() if f.suffix == '.mzML'
        ]
        max_workers = min(self.args.processes, len(tasks))
        with multiprocessing.Pool(max_workers) as pool:
            for i, stdout in enumerate(
                pool.imap_unordered(self._run_msgf, tasks)
            ):
                print(f"\n[MS-GF+ {i}]")
                print(stdout)
        return self

    def _run_msgf(self, data: tuple[bool, pathlib.Path]):
        is_decoy = data[0]
        ms_file = str(data[1].absolute())
        cmd = [
            *self.msgf_command,
            '-d', str(self.decoy_file if is_decoy else self.database_file),
            '-s', f'{self.args.mass_spec}/{ms_file}'
        ]

        if is_decoy:
            cmd.append('-o')
            cmd.append(f'{self.args.mass_spec}/{ms_file}_decoy.mzid')

        result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        if (code := result.returncode) != 0:
            err = f'MS-GF+ finished with non-zero return value: {code}'
            cli.exit(3, err)

        return result.stdout

    def _check_folder(self):
        mzid_folder = pathlib.Path(f'mzid/{self.database_type}')
        mzid_folder.mkdir(exist_ok=True, parents=True)

    def _build_msgf_command(self) -> list[str]:
        """Helper function that builds the reusable part of MSGF+ command
        based on the pipeline args.

        Returns
        -------
            list[str]
                A list of MSGF+ flags extracted from the pipeline args.
        """
        cmd: list[str] = [
            'java',
            f'-Xmx{self.args.xmx}',
            '-jar', f'{self.path}/dependencies/MSGF/MSGFPlus.jar',
            '-tda', '0',
            '-addFeatures', '1'
        ]
        pipeline_args = (
            'outdir', 'mode', 'mass_spec', 'processes', 'transcriptome'
        )
        for key, value in vars(self.args).items():
            if key not in pipeline_args and value is not None:
                cmd.append(f'-{key}')
                cmd.append(value)
        return cmd

    def peptide_filtering(self):
        """ DEPRECATED FOR NOW """
        print("Filtering peptide identification")

        # list of arguments passed by argparser
        arg_list = []
        for item in vars(self.args).items():
            item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode"]
            if item[0] not in item_list:
                arg_list.append("--%s %s" % (item[0], item[1]))
        arg_string = ""
        for i in range(len(arg_list)):
            arg_string += "%s " % arg_list[i]

        cmd_ms_filter = 'Rscript %s/workflow_filter.R %s--folder %s/' % (sys.path[0], arg_string, self.database_type)
        os.system(cmd_ms_filter)
        cmd_cat = 'cat %s/*pepseq.txt > %s/clustered_peptides.txt' % (self.database_type, self.database_type)
        os.system(cmd_cat)
