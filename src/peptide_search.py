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

from src.search_engine import SearchEngine
from src import utils


class PeptideSearch:
    def __init__(
        self,
        database_type,
        ms_files_folder,
        database_file,
        decoy_file,
        search_engine: SearchEngine,
        args
    ):
        self.database_type = database_type
        self.ms_folder = pathlib.Path(ms_files_folder)
        self.database_file = pathlib.Path(database_file)
        self.decoy_file = pathlib.Path(decoy_file)
        self.args = args
        self.search_engine = search_engine
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

        self.loop_search()

    def loop_search(self):
        if self.args.search_engine == 'CometMS':
            comet_db = pathlib.Path(
                f'{self.database_type}_comet_database.fasta'
            )
            utils.concat_fastas(
                comet_db,
                self.database_file,
                self.decoy_file
            )
            databases = [comet_db]
        else:
            databases = [self.database_file, self.decoy_file]

        for database in databases:
            self.search_engine.run(
                self.database_type,
                self.ms_folder,
                database
            )

        self.search_engine.save_to_pin(
            self.database_type,
            self.database_file,
            self.decoy_file
        )

        return self

    def _check_folder(self):
        mzid_folder = pathlib.Path(f'mzid/{self.database_type}')
        mzid_folder.mkdir(exist_ok=True, parents=True)

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
