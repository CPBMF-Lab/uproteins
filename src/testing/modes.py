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

import pandas as pd

from src.database import database_generator as dg
from src import peptide_search as ps
from src import results_new_approach as pms
from .testargs import TestArgs
from ..assembly import TranscriptAssembly, CompareTranscripts, ReadMapper
from ..master import Archives
from ..database import Database
from ..percolator import Decoy
from ..postprocess import PostPercolator, ExtendedInformation, PercolatorProcessing, AllSub, TSVConverter
from ..sequtils.orflib import AltCodons
from ..upstream import SDInspection
from ..pipelines import PostMSPipeline, ValidatePipeline


class TestingOutput(object):
    def __init__(self, args):

        self.args = args
        self.boolDict = {"TRUE": "skipped", "FALSE": "not tested"}
        self.assembly = self.boolDict[self.args.skip_assembly]
        self.db = self.boolDict[self.args.skip_db]
        self.ms = self.boolDict[self.args.skip_ms]
        self.postms = self.boolDict[self.args.skip_postms]


    def db_checkpoint(self):
        """ Check whether the databases are empty or not. """
        out = os.path.abspath(self.args.outdir)
        content = None
        if os.path.exists(f'{out}/transcriptome_database.fasta'):
            if os.stat(f"{out}/transcriptome_database.fasta").st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content

    def ms_checkpoint(self):
        """ Checks whether the clustered_peptides file is at its proper place and is not empty. This is the final result
         of the MS step. """
        out = os.path.abspath(self.args.outdir)
        content = None
        file = f'{out}/Transcriptome/clustered_peptides.txt'
        if os.path.exists(file):
            if os.stat(file).st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content

    def postms_checkpoint(self):
        out = os.path.abspath(self.args.outdir)
        content = None
        file = f'{out}/Results/both_ORFs.fasta'
        if os.path.exists(file):
            if os.stat(file).st_size != 0:
                content = "Filled"
            else:
                content = "Empty"
        return content


    def testing_stdout(self, mode):

        stdout = f"======================== Testing uProteInS {mode}=======================\n" \
                                "|\n" \
                                f"|- Assembly: {self.assembly}\n" \
                                "|\n" \
                                f"|- Database: {self.db}\n" \
                                "|\n" \
                                f"|- MS: {self.ms}\n" \
                                "|\n" \
                                f"|- PostMS: {self.postms}\n" \
                                "================================================================================\n"
        return stdout


class PipelineTesting(object):
    def __init__(self, args):
        self.args = args
        self.path = sys.path
        self.testFolder = f"{self.path[0]}/testkit"
        self.assemblyState = 'Not tested'
        self.databaseState = 'Not tested'
        self.MSState = 'Not tested'
        self.postMSState = 'Not tested'
        self.validateState = 'Not tested'

    def print_status(self):
        print(f'\n--------------\nTesting uProteInS different modes\nAssembly: {self.assemblyState}\nDatabase:'
              f' {self.databaseState}\nPeptide Search (MS mode): {self.MSState}\nPostMS: {self.postMSState}'
              f'\nValidate mode: {self.validateState}\n--------------\n')

    def run(self):
        assembly = AssemblyTesting(self.args)
        self.assemblyState = assembly.test()
        if self.assemblyState != 'SKIPPED':
            self.print_status()
        database = DatabaseTesting(self.args)
        self.databaseState = database.test()
        if self.databaseState != 'SKIPPED':
            self.print_status()
        search = MSTesting(self.args)
        self.MSState = search.test()
        if self.MSState != 'SKIPPED':
            self.print_status()
        if self.postMSState != 'SKIPPED':
            postms = PostMSTesting(self.args)
            self.postMSState = postms.test()
            self.print_status()
        validate = ValidateTesting(self.args)
        self.validateState = validate.test()
        self.print_status()


class AssemblyTesting(PipelineTesting):
    def __init__(self, args):
        super().__init__(args)
        self._fix_args()

    def _fix_args(self):
        self.args.gtf = f'{self.testFolder}/test_gtf.gtf'
        self.args.single = f'{self.testFolder}/testing_reads.fastq'
        self.args.genome = f'{self.testFolder}/testing_genome.fasta'
        self.args.strandness = 'F'
        self.args.reads1 = None
        self.args.reads2 = None
        self.args.gff_compare_path = None
        self.args.gffread_path = None

    def test(self):
        if self.args.skip_assembly == 'FALSE':
            self._map_reads()
            self._assemble_transcriptome()
            self.check_assembly()
        else:
            self.assemblyState = 'SKIPPED'
        return self.assemblyState

    def _map_reads(self):
        mapping = ReadMapper(self.args)
        mapping.create_indexes()
        mapping.align_reads()

    def _assemble_transcriptome(self):
        arxivs = Archives()
        assembly = TranscriptAssembly(self.args)
        assembly.assemble()
        assembly.create_gtf_list()
        gtf = assembly.merge_transcripts()
        arxivs.assembledGTF = gtf
        comparisons = CompareTranscripts(self.args, arxivs.assembledGTF)
        compare_dir = comparisons.run_gffcompare()
        arxivs.comparisonsDirectory = compare_dir
        sequences = comparisons.extract_sequences()
        arxivs.RNASequences = sequences

    def check_assembly(self):
        self.assemblyState = 'FAILED'
        with open('assembled.gtf', 'r') as ass:
            lines = ass.readlines()
            for line in lines:
                if 'gene-Rv0001' in line:
                    self.assemblyState = 'OK'


class DatabaseTesting(PipelineTesting):
    def __init__(self, args):
        super().__init__(args)
        self._fix_args()
        self.transcriptome = f'{self.testFolder}/transcripts.fasta'
        self._get_transcriptome()

    def _get_transcriptome(self):
        if not os.path.exists('HISAT'):
            cmd_dir = 'mkdir HISAT'
            os.system(cmd_dir)
        cmd_mv = f'cp {self.transcriptome} HISAT/.'
        os.system(cmd_mv)

    def _fix_args(self):
        self.args.genome = f'{self.testFolder}/genome_for_database.fasta'
        self.args.proteome = f'{self.testFolder}/proteome_for_database.fasta'
        self.args.transcriptome = True
        self.args.minsize = 30
        self.args.maxsize = 300
        self.args.starts = ['ATG', 'ATT', 'TTG', 'GTG']
        self.args.stops = ['TAA', 'TAG', 'TGA']

    def test(self):
        if self.args.skip_db != 'TRUE':
            print("Generating the genome database.")
            db = Database(self.args)
            db.translate()
            # print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
            genome_db = dg.Database("genome_ORFs.fasta", self.args.proteome, "genome")
            # genome_db.mark_annotated()
            genome_db.unify()
            print("Genome database generated.")
            if self.args.transcriptome:
                print("Generating the transcriptome database.")
                transcriptome_db = dg.Database("transcriptome_ORFs.fasta", self.args.proteome, "transcriptome")
                # transcriptome_db.mark_annotated()
                transcriptome_db.unify()
                print("Transcriptome database generated.")
            self._check()
        else:
            self.databaseState = 'SKIPPED'
        return self.databaseState

    def _check(self):
        self.databaseState = 'FAILED'
        rna = 'FAILED'
        dna = 'FAILED'
        with open('transcriptome_database.fasta', 'r') as db:
            lines = db.readlines()
            for line in lines:
                if 'tORF_gene-Rv0001' in line:
                    rna = 'OK'
                    break
        with open('genome_database.fasta', 'r') as db:
            lines = db.readlines()
            for line in lines:
                if 'gORF__6_398-493' in line:
                    dna = 'OK'
                    break
        if rna == 'OK' and dna == 'OK':
            self.databaseState = 'OK'


class MSTesting(PipelineTesting):
    def __init__(self, args):
        super().__init__(args)
        self._fix_args()
        self.__move_test_files()

    def _fix_args(self):
        self.args.Mass_spec = f'{self.testFolder}/mzml/'
        self.args.transcriptome = True

    def __move_test_files(self):
        files = ['genome_database.fasta', 'transcriptome_database.fasta']
        for file in files:
            os.system(f'cp {self.testFolder}/{file} .')

    def test(self):
        if self.args.skip_ms == "FALSE":
            genome = ps.PeptideSearch("Genome", self.args.Mass_spec, "genome_database.fasta", self.args)
            genome.peptide_identification()
            genome_decoy = Decoy(db="genome_database.fasta", db_type="Genome")
            genome_decoy.reverse_sequences().to_fasta()
            genome_decoy_search = ps.PeptideSearch("Genome", self.args.Mass_spec, "Genome/Percolator/Genome_decoy.fasta",
                                                   self.args, decoy=True)
            genome_decoy_search.peptide_identification()
            if self.args.transcriptome:
                transcriptome = ps.PeptideSearch("Transcriptome", self.args.Mass_spec, "transcriptome_database.fasta", self.args)
                transcriptome.peptide_identification()
                transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
                transcriptome_decoy.reverse_sequences().to_fasta()
                transcriptome_decoy_search = ps.PeptideSearch("Transcriptome", self.args.Mass_spec,
                                                              "Transcriptome/Percolator/Transcriptome_decoy.fasta", self.args,
                                                              decoy=True)
                transcriptome_decoy_search.peptide_identification()
            self._check()

        else:
            self.MSState = 'SKIPPED'
        return self.MSState

    def _check(self):
        self.MSState = 'FAILED'
        target = '20140719_H37Rv_20140718_5ug_120min_top8_1.mzid'
        decoy = '20140719_H37Rv_20140718_5ug_120min_top8_1.mzML_decoy.mzid'

        def check_size(folder):
            # if os.path.getsize(f'{folder}/{target}') > 20000 and os.path.getsize(f'{folder}/{decoy}') > 20000:
            if os.path.exists(f'{folder}/{target}') and os.path.exists(f'{folder}/{decoy}'):
                state = 'OK'
            else:
                state = 'FAILED'
            return state

        genome_state = check_size('Genome')
        transcriptome_state = check_size('Transcriptome')
        if genome_state == 'OK' and transcriptome_state == 'OK':
            self.MSState = 'OK'


class PostMSTesting(PipelineTesting):
    def __init__(self, args):
        super().__init__(args)
        self.postMSKit = f'{self.testFolder}/postms'
        self.__check_folders()
        self._check_transcriptome_folder()
        self._fix_args()
        self._add_test_kit()

    def _fix_args(self):
        # self.args.genome = f'{self.testFolder}/genome_for_database.fasta'
        # self.args.genome = f'{self.testFolder}/testing_genome.fasta'
        self.args.genome = f'{self.postMSKit}/genome.fasta'
        self.args.proteome = f'{self.testFolder}/proteome_for_database.fasta'
        self.args.gff = f'{self.testFolder}/test_gff.gff'
        self.args.pep = 0.5
        self.args.qvalue = 0.5
        self.args.starts = ['ATG', 'TTG', 'ATT', 'GTG']
        self.args.stops = ['TAA', 'TAG', 'TGA']
        self.args.maxsize = 300
        self.args.rrna = f'{self.testFolder}/rrna.fna'
        self.args.transcriptome = True
        self.newArgs = self.args

    def __check_folders(self):
        folders = ['Genome', 'Transcriptome', 'HISAT', 'Genome/Percolator', 'Transcriptome/Percolator']
        for folder in folders:
            if not os.path.exists(folder):
                os.system(f'mkdir {folder}')

    def _add_test_kit(self):

        directions = {'genome_database.fasta': '.', 'Genome_decoy.fasta': 'Genome/Percolator/.',
                      'transcriptome_database.fasta': '.', 'Transcriptome_decoy.fasta': 'Transcriptome/Percolator/.',
                      'Transcriptome/20140719_H37Rv_20140718_5ug_120min_top8_1.mzid': 'Transcriptome/.',
                      'Transcriptome/20140719_H37Rv_20140718_5ug_120min_top8_1.mzML_decoy.mzid': 'Transcriptome/.',
                      'Genome/20140719_H37Rv_20140718_5ug_120min_top8_1.mzid': 'Genome/.',
                      'Genome/20140719_H37Rv_20140718_5ug_120min_top8_1.mzML_decoy.mzid': 'Genome/.',
                      'assembled.gtf': '.', 'transcripts.fasta': 'HISAT/.', }

        def move_cmd(item, destination):
            cmd = f'cp {self.postMSKit}/{item} {destination}'
            os.system(cmd)

        for file in directions:
            move_cmd(file, directions[file])

    @staticmethod
    def _check_transcriptome_folder():
        if not os.path.exists("HISAT"):
            os.system('mkdir HISAT')

    def test(self):
        if self.args.skip_postms == 'FALSE':
            genome = PostMSPipeline(args=self.args, filetype='genome', folder='Genome', qvalue=0.27, testing=True)
            genome.run()
            if self.args.transcriptome:
                transcriptome = PostMSPipeline(args=self.newArgs, filetype='transcriptome', folder='Transcriptome',
                                               qvalue=0.27, testing=True)
                transcriptome.run()
            # genome_perc = PercolatorProcessing("Genome", filetype="genome")
            # genome_perc.create_metafiles().convert_to_pin()
            # all_sub = AllSub(filetype='genome', folder="Genome")
            # all_sub.modify_decoy()
            # all_sub.remove_annotated()
            # genome_perc.percolate()
            # if self.args.transcriptome:
            #     rna_perc = PercolatorProcessing("Transcriptome", filetype="transcriptome")
            #     rna_perc.create_metafiles().convert_to_pin()
            #     rna_all_sub = AllSub(filetype='transcriptome', folder='Transcriptome')
            #     rna_all_sub.modify_decoy()
            #     rna_all_sub.remove_annotated()
            #     rna_perc.percolate()
            # genome_tsv = TSVConverter("Genome")
            # genome_tsv.convert_files()
            # if self.args.transcriptome:
            #     transcriptome_tsv = TSVConverter("Transcriptome")
            #     transcriptome_tsv.convert_files()
            # genome_filter = PostPercolator(self.args, "Genome", filetype='genome')
            # genome_filter.convert_output()
            # genome_filter.get_coordinates_genome()
            # genome_filter.filter_novel()
            # genome_filter.unique_peptides()
            # genome_filter.msgf_info()
            # genome_filter.protein_seqs()
            # genome_filter.add_coordinates()
            # genome_alts_pre_rf = AltCodons(file='Genome/post_perc/genome_results_02.txt', genome=self.args.genome,
            #                                maxsize=self.args.maxsize)
            # genome_alts_pre_rf.extend_orfs(args=self.args)
            # if self.args.rrna is not None:
            #     genome_rbs = SDInspection(self.args, filetype="genome", folder="Genome",
            #                               alternatives=genome_alts_pre_rf.alternatives)
            #     alts = genome_rbs.get_free_energy()
            #     genome_alts_pre_rf.alternatives = alts
            # genome_alts_pre_rf.sort_by_coordinates()
            # genome_alts_pre_rf.sort_by_atg()
            # genome_alts_pre_rf.sort_by_shine()
            # genome_alts_pre_rf.sort_by_peptides()
            # priorities = genome_alts_pre_rf.get_priorities()
            # ext = ExtendedInformation(folder='Genome', filetype='genome', alternatives=genome_alts_pre_rf.alternatives)
            # ext.filter_alternatives(priorities)
            # ext.extract_spectra()
            # if self.args.transcriptome:
            #     transcriptome_filter = PostPercolator(self.args, 'Transcriptome', filetype='transcriptome')
            #     transcriptome_filter.convert_output()
            #     transcriptome_filter.get_coordinates_rna()
            #     transcriptome_filter.filter_novel()
            #     transcriptome_filter.unique_peptides()
            #     transcriptome_filter.msgf_info()
            #     transcriptome_filter.protein_seqs()
            #     transcriptome_filter.add_coordinates()
            #     transcriptome_alts_pre_rf = AltCodons(file='Transcriptome/post_perc/transcriptome_results_02.txt',
            #                                           genome=self.args.genome, maxsize=self.args.maxsize,
            #                                           transcriptome_gff='assembled.gtf', assembly="HISAT/transcripts.fasta")
            #     transcriptome_alts_pre_rf.extend_orfs(args=self.args)
            #     if self.args.rrna is not None:
            #         transcriptome_rbs = SDInspection(self.args, filetype='transcriptome', folder='Transcriptome',
            #                                          alternatives=transcriptome_alts_pre_rf.alternatives,
            #                                          transcriptome_gff='assembled.gtf',
            #                                          transcripts="HISAT/transcripts.fasta")
            #         rna_alts = transcriptome_rbs.get_free_energy()
            #         transcriptome_alts_pre_rf.alternatives = rna_alts
            #     transcriptome_alts_pre_rf.sort_by_coordinates()
            #     transcriptome_alts_pre_rf.sort_by_atg()
            #     transcriptome_alts_pre_rf.sort_by_shine()
            #     transcriptome_alts_pre_rf.sort_by_peptides()
            #     rna_priorities = transcriptome_alts_pre_rf.get_priorities()
            #     rna_ext = ExtendedInformation(folder='Transcriptome', filetype='transcriptome',
            #                                   alternatives=transcriptome_alts_pre_rf.alternatives)
            #     rna_ext.filter_alternatives(rna_priorities)
            #     rna_ext.extract_spectra()
            self._check()
        else:
            self.postMSState = 'SKIPPED'
        return self.postMSState

    def _check(self):
        self.postMSState = 'FAILED'

        def inspect_result(pattern, folder, filetype):
            state = 'FAILED'
            results = f'{folder}/post_perc/{filetype}_results_04.txt'
            if os.path.exists(results):
                with open(results, 'r') as data:
                    lines = data.readlines()
                    for line in lines:
                        if pattern in line:
                            state = 'OK'
                            break
            return state

        transcriptome_state = inspect_result(pattern='controllerType=0', folder='Transcriptome',
                                             filetype='transcriptome')
        genome_state = inspect_result(pattern='controllerType=0', folder='Genome', filetype='genome')
        if transcriptome_state == 'OK' and genome_state == 'OK':
            self.postMSState = 'OK'


class ValidateTesting(PipelineTesting):
    def __init__(self, args):
        super().__init__(args)
        self.validationKit = f'{self.testFolder}/validate'
        self.__check_folders()
        self.__add_testing_files()

    @staticmethod
    def __check_folders():
        folders = ['Genome', 'Genome/Percolator', 'Genome/Percolator/for_predicting',
                   'Transcriptome', 'Transcriptome/Percolator', 'Transcriptome/Percolator/for_predicting',
                   "Genome/post_perc", "Transcriptome/post_perc"]
        for folder in folders:
            print(folder)
            # if not os.path.exists(folder):
            os.system(f'mkdir {folder}')

    def __add_testing_files(self):
        pin_files = os.listdir(f'{self.validationKit}/pin_files')
        files = {"genome_results_04.txt": "Genome/post_perc", "transcriptome_results_04.txt": "Transcriptome/post_perc"}
        for pin in pin_files:
            os.system(f'cp {self.validationKit}/pin_files/{pin} Transcriptome/Percolator/.')
            os.system(f'cp {self.validationKit}/pin_files/{pin} Genome/Percolator/.')
        for file in files:
            os.system(f'cp {self.validationKit}/{file} {files[file]}/.')

    def test(self):
        if self.args.skip_validation == 'FALSE':
            validation = ValidatePipeline(args=self.args, testing=True)
            validation.validate_genome()
            validation.validate_transcriptome()
            self._check()
        else:
            self.validateState = 'SKIPPED'
        return self.validateState

    def _check(self):
        if os.path.exists(f'Transcriptome/post_perc/transcriptome_results_05.txt'):

        # df = pd.read_csv(f'Transcriptome/post_perc/transcriptome_results_05.txt', sep='\t')
        # files = df["RenamedFiles"].tolist()
        # if '20150904_cdd_KO_CITOPL_50ug_MP_9' in files:
            self.validateState = "OK"
        else:
            self.validateState = "FAILED"

class PipelineTestingOld(object):
    def __init__(self, outdir, skip_assembly, skip_db, skip_ms, skip_postms):
        self.outdir = outdir
        self.path = sys.path
        self.testFolder = f"{self.path[0][:-3]}testing_kit"
        self.args = TestArgs(self.testFolder, self.outdir, skip_assembly, skip_db, skip_ms, skip_postms)

    def test_assembly(self):
        args = self.args
        rna_seq = dg.Assembly(args)
        rna_seq.denovo_assembly()

    def summarize_assembly(self):
        cmd_rename = f'mv {self.args.outdir}/trinity/Trinity.fasta {self.args.outdir}/trinity/Trinity_cut.fasta'
        os.system(cmd_rename)
        cmd_new_rna = f'cp {self.testFolder}/test_transcriptome.fasta {self.args.outdir}/trinity/Trinity.fasta'
        os.system(cmd_new_rna)

    def test_database(self):
        args = self.args
        genome = dg.OrfPrediction(args)
        genome.identify_orfs()
        print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
        genome_db = dg.Database("genome_ORFs.fasta", args.proteome, "genome")
        genome_db.blast_to_Proteome()
        genome_db.extract_entries()
        cmd_cat_entries = 'cat genome_entries.txt genome_entries_short.txt > genome_entries_cat.txt'
        os.system(cmd_cat_entries)
        genome_db.remove_entries()
        genome_db.create_custom()
        if args.transcriptome:
            print("\nGenome database generated. Now performing steps to generate the TRANSCRIPTOME database.\n")
            transcriptome_db = dg.Database("transcriptome_ORFs.fasta", args.proteome, "transcriptome")
            transcriptome_db.blast_to_Proteome()
            transcriptome_db.extract_entries()
            cmd_cat_entries_rna = 'cat transcriptome_entries.txt transcriptome_entries_short.txt > ' \
                                  'transcriptome_entries_cat.txt'
            os.system(cmd_cat_entries_rna)
            transcriptome_db.remove_entries()
            transcriptome_db.create_custom()
        print("Database generation step complete. Look for databases in %s" % args.outdir)

    def test_ms(self):
        """ Tests the peptide search mode of uProteInS. """
        args = self.args.msArgs
        ms_folder = args.get_ms_folder(self.args.path)
        genome = ps.PeptideSearch("Genome", ms_folder, "genome_database.fasta", args)
        # genome.peptide_identification()
        genome.peptide_filtering()
        if args.transcriptome:
            transcriptome = ps.PeptideSearch("Transcriptome", ms_folder, "transcriptome_database.fasta", args)
            # transcriptome.peptide_identification()
            transcriptome.peptide_filtering()

    def test_postms(self):
        args = self.args
        pms.merge_databases()
        genome_results = pms.Results("Genome", "genome", "genome_database.fasta")
        genome_results.rename_files()
        genome_results.write_results()
        genome_results.merge()
        genome_results.add_orf()
        genome_results.add_orf_sequence()
        genome_results.total_spec_count()
        genome_results.genome_coordinates()
        genome_results.check_unique(args)
        genome_results.coverage()
        genome_results.summarize_results()
        genome_results.summarize_dna_rna()
        # genome_results.create_fast("Genome")
        # genome_results.unique_orfs()

        # genome_results.minimal_results()
        # genome_results.minimum_runs_compact("15", "genome")
        if args.genbank is not None:
            genome_ncrnas = pms.GenomicContext(args, "genome", "Genome")
            genome_ncrnas.match_features()
            genome_ncrnas.add_to_df()
        if args.transcriptome:
            transcriptome_results = pms.Results("Transcriptome", "transcriptome",
                                                "transcriptome_database.fasta")
            transcriptome_results.rename_files()
            transcriptome_results.write_results()
            transcriptome_results.merge()
            transcriptome_results.add_orf()
            transcriptome_results.add_orf_sequence()
            transcriptome_results.total_spec_count()
            transcriptome_results.genome_coordinates()
            transcripts = pms.TranscriptLocalizer(args)
            transcripts.find_transcript()
            transcripts.fix_nucmer()
            transcripts.add_coords()
            transcriptome_results.check_unique(args)
            transcriptome_results.coverage()
            transcriptome_results.summarize_results()
            transcriptome_results.summarize_dna_rna()
            # pms.create_venn()
            # transcriptome_results.add_names()
            # genome_results.create_fast("Both")
            # transcriptome_results.create_fast("Rna")
            # transcriptome_results.merge_both_results()
            # transcriptome_results.unique_orfs()

            both_results = pms.Results("Genome", "Both", "_")
            """
            both_results.minimal_results()
            both_results.minimum_runs_compact("15", "both")
            transcriptome_results.minimal_results()
            transcriptome_results.minimum_runs_compact("15", "transcriptome")
            """
            # if args.genbank is not None:
            #     rna_features = pms.GenomicContext(args, "transcriptome", "Transcriptome")
            #     rna_features.match_features()
            #     rna_features.add_to_df()
            #     both_features = pms.GenomicContext(args, "both", "Genome")
            #     both_features.match_features()
            #     both_features.add_to_df()
