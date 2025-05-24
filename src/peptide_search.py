import os
import sys
import pathlib
import subprocess
import multiprocessing

from src import cli


class PeptideSearch(object):
    def __init__(self, database_type, ms_files_folder, orf_file, args, decoy=False):
        self.database_type = database_type
        self.ms_files_folder = ms_files_folder
        self.orf_file = orf_file
        self.args = args
        self.path = sys.path[0]
        self.decoy = decoy

    def peptide_identification(self):
        print("\nPerforming peptide search using %s database\n" % self.database_type)
        cmd_dir_ms = 'mkdir %s' % self.database_type
        os.system(cmd_dir_ms)
        cmd_copy_orf = 'cp %s %s/' % (self.orf_file, self.database_type)
        os.system(cmd_copy_orf)
        cmd_copy_db = f'cp {self.orf_file} {os.path.abspath(self.ms_files_folder)}/.'
        os.system(cmd_copy_db)

        # list of arguments passed by argparser
        # arg_list = []
        # for item in vars(self.args).items():
        #     item_list = [None, "Mass_spec", "outdir", "Transcriptome", "mode"]
        #     if item[0] not in item_list:
        #         arg_list.append("--%s %s" % (item[0], item[1]))
        # arg_string = ""
        # for i in range(len(arg_list)):
        #     arg_string += "%s " % arg_list[i]

        # cmd_pep_search = 'Rscript %s/mzid_workflow.R %s--database %s --folder %s'\
        #                  % (sys.path[0], arg_string, os.path.abspath(self.orf_file), self.ms_files_folder)
        # os.system(cmd_pep_search)
        self.parallel_search()
        cmd_move = 'mv %s/*.mzid %s/' % (self.ms_files_folder, self.database_type)
        os.system(cmd_move)

    def parallel_search(self):
        self.msgf_command = self._build_msgf_command()
        files = [i for i in os.listdir(os.path.abspath(self.ms_files_folder)) if i.endswith('mzML')]
        max_workers = min(self.args.processes, len(files))
        pool = multiprocessing.Pool(max_workers)
        with multiprocessing.Pool(max_workers) as pool:
            for file, stdout in pool.imap_unordered(self._run_msgf, files):
                print(f"[MSGF+ {file}]")
                print(stdout)
        return self

    def _run_msgf(self, file):
        self._check_folder()
        database = pathlib.Path(self.orf_file).absolute()
        cmd = [
            *self.msgf_command,
            '-d', str(database),
            '-s', f'{self.args.mass_spec}/{file}'
        ]

        if self.decoy:
            cmd.append('-o')
            cmd.append(f'{self.args.mass_spec}/{file}_decoy.mzid')

        result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        if (code := result.returncode) != 0:
            err = f'MSGF+ finished with non-zero return value: {code}'
            cli.exit(3, err)

        return file, result.stdout

    def _check_folder(self):
        if not os.path.exists('mzid'):
            os.system('mkdir mzid')
        if self.args.transcriptome and not os.path.exists('mzid/Transcriptome'):
            os.system('mkdir mzid/Transcriptome')
        if not os.path.exists('mzid/Genome'):
            os.system('mkdir mzid/Genome')

    def _build_msgf_command(self) -> list[str]:
        """Helper function that builds the reusable part of MSGF+ command
        based on the pipeline args.

        Returns
        -------
            list[str]
                A list of MSGF+ flags extracted from the pipeline args.
        """
        database = pathlib.Path(self.orf_file).absolute()
        cmd: list[str] = [
            'java',
            '-Xmx48G',
            '-jar', f'{self.path}/dependencies/MSGF/MSGFPlus.jar',
            '-d', str(database),
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
