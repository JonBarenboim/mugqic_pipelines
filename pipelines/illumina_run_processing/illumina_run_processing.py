#!/usr/bin/env python

# Python Standard Modules
import logging
import math
import os
import re
import sys
import itertools
import xml.etree.ElementTree as XML
# TODO import future

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *

from bfx import bvatools
from bfx import bwa
from bfx import gatk
from bfx import igvtools
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import tools
from pipelines import common

log = logging.getLogger(__name__)

class RunInfoRead(object):
    """ Model of a read from the Illumina sequencer.
        Those attributes can be found in the RunInfo.xml file.
    """

    def __init__(self, number, nb_cycles, is_index):
        self._number = number
        self._nb_cycles = nb_cycles
        self._is_index = is_index

    @property
    def number(self):
        return self._number

    @property
    def nb_cycles(self):
        return self._nb_cycles

    @property
    def is_index(self):
        return self._is_index


class IlluminaRunProcessing(common.Illumina):
    """Main object of the Illumina Run Processing Pipeline."""

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_illumina_raw_readset_files(
                                                               self.output_dir,
                                                               "PAIRED_END" if self.is_paired_end else "SINGLE_END",
                                                               self.nanuq_readset_file,
                                                               self.casava_sheet_file,
                                                               self.args.lane_number,
                                                               config.param('global', 'default_species_genome'),
                                                               config.param('global', 'genomes_home')
                                                               )
        return self._readsets

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            self._is_paired_end = len([read_info for read_info in self.read_infos if (not read_info.is_index)]) > 1
        return self._is_paired_end

    @property
    def copy_inputs(self):
        if not hasattr(self, "_copy_inputs"):
            self._copy_inputs = []
        return self._copy_inputs

    @property
    def run_id(self):
        """ The run id from the run folder.
            Supports both default folder name configuration and GQ's globaly unique name convention.
        """
        if not hasattr(self, "_run_id"):
            if (re.search(".*_\d+HS\d\d[AB]", self.run_directory)):
                m = re.search(".*\/(\d+_[^_]+_\d+_[^_]+_(\d+)HS.+)", self.run_directory)
                self._run_id = m.group(2)
            elif (re.search(".*\d+_[^_]+_\d+_.+", self.run_directory)):
                m = re.search(".*\/(\d+_([^_]+_\d+)_.*)", self.run_directory)
                self._run_id = m.group(2)
            else:
                log.warn("Unsupported folder name: " + self.run_directory)

        return self._run_id

    @property
    def run_directory(self):
        return self.args.run_directory

    @property
    def lane_number(self):
        return self.args.lane_number

    @property
    def casava_sheet_file(self):
        return self.args.casava_sheet_file if (self.args.casava_sheet_file) else self.run_directory + os.sep + "SampleSheet.nanuq.csv"

    @property
    def nanuq_readset_file(self):
        return self.args.nanuq_readset_file if (self.args.nanuq_readset_file) else self.run_directory + os.sep + "run.nanuq.csv"

    @property
    def number_of_mismatches(self):
        return self.args.number_of_mismatches if (self.args.number_of_mismatches) else 1

    @property
    def first_index(self):
        return self.args.first_index if (self.args.first_index) else 1

    @property
    def last_index(self):
        return self.args.last_index if (self.args.last_index) else 999

    @property
    def steps(self):
        return [
            self.index_count,
            self.fastq,
            self.md5,
            self.qc_graphs,
            self.blast,
            self.bwa_mem_picard_sort_sam,
            self.picard_mark_duplicates,
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.bam_md5,
            self.start_copy_notification,
            self.copy,
            self.end_copy_notification
        ]

    @property
    def read_infos(self):
        if not hasattr(self, "_read_infos"):
            self._read_infos = self.parse_run_info_file()
        return self._read_infos

    @property
    def argparser(self):
        if not hasattr(self, "_argparser"):
            self._argparser = super(common.Illumina, self).argparser

            # overiding r parameter from readset filename to run dir
            self._argparser.add_argument("-r", "--run", help="run directory", required=True, dest="run_directory")
            self._argparser.add_argument("-p", "--lane", help="lane number", type=int, required=True, dest="lane_number")
            self._argparser.add_argument("-n", "--readsets", help="nanuq run sheet", type=file, required=False, dest="nanuq_readset_file")
            self._argparser.add_argument("-i", help="illumina casava sheet", type=file, required=False, dest="casava_sheet_file")
            self._argparser.add_argument("-x", help="first index base to use for demultiplexing", type=int, required=False, dest="first_index")
            self._argparser.add_argument("-y", help="last index base to use for demultiplexing", type=int, required=False, dest="last_index")
            self._argparser.add_argument("-m", help="number of index mistmaches allowed for demultiplexing", type=int, required=False, dest="number_of_mismatches")

        return self._argparser


    def index_count(self):
        """ Generate a file with all the indexes found in the index-reads of the run.
            The file nammed "RUNFOLDER_LANENUMBER.metrics" will be in saved in the output directory.
        """
        jobs = []

        mask = ""
        index_length = self.getSequencerIndexLength()

        for read in self.read_infos:
            if (read.is_index):
                mask += str(index_length) + 'B'
                break
            else:
                mask += str(read.nb_cycles) + 'T'

        if (index_length == 0):
            log.info("No Indexes, *NOT* Generating index counts")
        else:
            input = self.run_directory + os.sep + "RunInfo.xml"
            output = self.run_directory + os.sep + os.path.basename(self.run_directory) + "_" + str(self.lane_number) + '.metrics'
            job = Job([input], [output], [["index_count", "module_java"]], name="index_" + self.run_id + str(self.lane_number))
            job.command = """\
java -Djava.io.tmpdir={tmp_dir}\\
 {java_other_options}\\
 -Xmx{ram}\\
 -jar {jar}\\
 MAX_MISMATCHES={mistmaches}\\
 NUM_PROCESSORS={threads}\\
 BARCODE_FILE={barcode_file}\\
 BASECALLS_DIR={basecalls_dir}\\
 LANE={lane_number}\\
 READ_STRUCTURE={read_structure}\\
 METRICS_FILE={output}\\
 TMP_DIR={tmp_dir}""".format(
                tmp_dir = config.param('index_count', 'tmp_dir'),
                java_other_options = config.param('index_count', 'java_other_options'),
                ram = config.param('index_count', 'ram'),
                jar = config.param('index_count', 'jar'),
                mistmaches = self.number_of_mismatches,
                threads = config.param('index_count', 'threads'),
                barcode_file = config.param('index_count', 'barcode_file'),
                basecalls_dir = self.run_directory + os.sep + config.param('index_count', 'basecalls_dir'),
                lane_number = self.lane_number,
                read_structure = mask,
                output = output
            )
            jobs.append(job)

        return jobs

    def fastq(self):
        """ Launch fastq generation from Illumina raw data using BCL2FASTQ conversion software.
            Demultiplexing, if needed, takes place in this step.
        """
        jobs = []

        input = self.casava_sheet_file

        outputs = [readset.fastq1 for readset in self.readsets]
        if (self.is_paired_end) :
            outputs += [readset.fastq2 for readset in self.readsets]

        output_dir = self.output_dir + os.sep + "Unaligned." + str(self.lane_number)

        job = Job([input], outputs, name="fastq_" + self.run_id + "_" + str(self.lane_number))
        job.command = "cd {unaligned_folder} && make -j {threads}".format(
            threads = config.param('index_count', 'threads'),
            unaligned_folder = output_dir
        )
        jobs.append(job)

        # TODO
        #if (not job.is_up2date()):
        if (True):
            self.validateBarcodes()
            mask = self.get_mask()
            casava_sheet_prefix = config.param('fastq', 'casava_sample_sheet_prefix')
            demultiplexing = False
            command = """\
module load {module} &&
configureBclToFastq.pl
 --input-dir {input_dir}\\
 --output-dir {output_dir}\\
 --tiles {tiles}\\
 --sample-sheet {sample_sheet}\\
 --fastq-cluster-count 0\\""".format(
                module = config.param('fastq', 'module_bcl_to_fastq'),
                input_dir = self.run_directory + os.sep + config.param('fastq', 'basecalls_dir'),
                output_dir = output_dir,
                tiles = "s_" + str(self.lane_number),
                sample_sheet = self.output_dir + os.sep + casava_sheet_prefix + str(self.lane_number) + ".csv"
            )
            if (re.search("I", mask)):
                demultiplexing = True
                command += " --mismatches {number_of_mismatches} --use-bases-mask {mask}".format(
                    number_of_mismatches = self.number_of_mismatches,
                    mask = mask
                )
            # TODO output command directly in bash
            log.info(command)

            notification_command = config.param('fastq', 'notification_command', required=False)
            if (notification_command):
                notification_command =notification_command.format(
                        output_dir = self.output_dir,
                        number_of_mismatches = self.number_of_mismatches if (demultiplexing) else "-",
                        lane_number = self.lane_number,
                        mask = mask if (demultiplexing) else "-",
                        technology = config.param('fastq', 'technology'),
                        run_id = self.run_id
                )
                # TODO output notification_command directly in bash
                log.info(notification_command)
        return jobs

    def md5(self):
        """ Create md5 checksum files for the fastq using the system 'md5sum' util.
            One checksum file will be created for each fastq.
        """
        jobs = []
        for readset in self.readsets:
            inputs = [readset.fastq1]
            outputs = [readset.fastq1 + ".md5"]
            command = "md5sum -b " + readset.fastq1 + " > " + readset.fastq1 + ".md5"

            # Second read in paired-end run
            if (readset.fastq2):
                inputs.append(readset.fastq2)
                outputs.append(readset.fastq2 + ".md5")
                command += " && md5sum -b " + readset.fastq2 + " > " + readset.fastq2 + ".md5"
            job = Job(inputs, outputs, name="md5_" + readset.name, command=command)
            jobs.append(job)
        return jobs

    def qc_graphs(self):
        """ Generate some QC Graphics and a summary XML file using BVATools.
            Files are created in a 'qc' subfolder of the fastq directory.
            Examples of output graphic:
                - Per cycle qualities, sequence content and sequence length;
                - Known sequences (adaptors);
                - Abundant Duplicates;
        """
        jobs = []

        for readset in self.readsets:
            output_dir = os.path.dirname(readset.fastq1) + os.sep + "qc"
            region_name = readset.name + "_" + readset.index + "_L00" + readset.lane

            job = concat_jobs([
                Job(command="mkdir -p " + output_dir),
                bvatools.readsqc(
                    readset.fastq1, 
                    readset.fastq2,
                    "FASTQ",
                    region_name, 
                    output_dir
                )]
            )

            job.name = "qc." + readset.name
            jobs.append(job)

        return jobs

    def blast(self):
        """ Run blast on a subsample of the reads to find the 20 most frequent hits.
            The "runBlat.sh" util from MUGQIC Tools is used.
            The number of reads to subsample can be configured by sample or for the whole lane.
            The output will be in the "Blast_sample" folder, under the Unaligned folder.
        """
        jobs = []

        nb_blast_to_do = config.param('blast', 'nb_blast_to_do', type="posint")
        is_nb_blast_per_lane = config.param('blast', 'is_nb_for_whole_lane', type="boolean")

        if (is_nb_blast_per_lane):
            nb_blast_to_do = int(nb_blast_to_do) / len(self.readsets)


        for readset in self.readsets:
            output_prefix = os.path.join(self.output_dir,
                                         "Unaligned." + readset.lane,
                                         "Blast_sample",
                                         readset.name + "_" + readset.index + "_L00" + readset.lane)
            inputs = [readset.fastq1, readset.fastq2]
            output = output_prefix + '.R1.RDP.blastHit_20MF_species.txt'
            command = "runBlast.sh " + str(nb_blast_to_do) + " " + output_prefix + " " + readset.fastq1 + " "
            if (readset.fastq2):
                command += readset.fastq2
            job = concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(output)),
                Job(inputs, [output], [["blast", "module_mugqic_tools"]], command=command)
            ], name= "blast_" + readset.name)

            jobs.append(job)

        return jobs

    def bwa_mem_picard_sort_sam(self):
        """
        """
        jobs = []

        for readset in [readset for readset in self.readsets if (readset.bam)]:
            output = readset.bam + ".sorted.bam"

            job = concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(output)),
                pipe_jobs([
                    bwa.mem(
                        readset.fastq1,
                        readset.fastq2,
                        read_group="'@RG" + \
                            "\tID:" + readset.name + \
                            "\tSM:" + readset.sample.name + \
                            ("\tLB:" + readset.library if readset.library else "") + \
                            ("\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                            ("\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=readset.aligner_reference_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        output,
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam." + readset.name)

            jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        """ Runs Picard mark duplicates util on the previously generated bam file."""
        jobs = []
        for readset in [readset for readset in self.readsets if (readset.bam)]:
            input_file_prefix = readset.bam + '.sorted.'
            input =  input_file_prefix + "bam"
            output = input_file_prefix + "dup.bam"
            metrics_file = readset.bam + "dup.metrics"

            job = picard.mark_duplicates([input], output, metrics_file)
            job.name = "picard_mark_duplicates." + readset.name
            jobs.append(job)
        return jobs

    def metrics(self):
        """"""
        jobs = []
        for readset in [readset for readset in self.readsets if (readset.bam)]:
            input_file_prefix = readset.bam + '.sorted.'
            input = input_file_prefix + "bam"

            job = picard.collect_multiple_metrics(input, input_file_prefix + "all.metrics", reference_sequence=readset.reference_file)
            job.name = "picard_collect_multiple_metrics." + readset.name
            jobs.append(job)

            if readset.beds:
                job = bvatools.depth_of_coverage(
                    input, 
                    input_file_prefix + "coverage.tsv", 
                    bvatools.resolve_readset_coverage_bed(readset), 
                    config.param('bvatools_depth_of_coverage', 'other_options', required=False),
                # reference_genome=readset.reference_file
                )

                job.name = "bvatools_depth_of_coverage." + readset.name
                jobs.append(job)

        return jobs

    def picard_calculate_hs_metrics(self):
        """"""
        jobs = []

        created_interval_lists = []

        for readset in [readset for readset in self.readsets if (readset.bam)]:
            coverage_bed = bvatools.resolve_readset_coverage_bed(readset)
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

                input_file_prefix = readset.bam + ".sorted.dup."
                job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "onTarget.tsv", interval_list)
                if not interval_list in created_interval_lists:
                    job = concat_jobs([tools.bed2interval_list(None, coverage_bed, interval_list), job])
                    created_interval_lists.append(interval_list)

                job.name = "picard_calculate_hs_metrics." + readset.name
                jobs.append(job)
        return jobs

    def bam_md5(self):
        """ Create md5 checksum files for the previously generated BAM (and BAI) files using the system 'md5sum' util."""
        jobs = []
        for readset in [readset for readset in self.readsets if (readset.bam)]:
            input_bai = readset.bam + ".sorted.bai"
            input_bam = readset.bam + ".sorted.bam"
            output_bai = input_bai + ".md5"
            output_bam = input_bam + ".md5"
            command = "md5sum -b " + input_bai + " > " + output_bai + " && md5sum -b " + input_bam + " > " + output_bam

            job = Job([input_bam], [output_bai, output_bam], name="bmd5_" + readset.name, command=command)
            jobs.append(job)
        return jobs

    def start_copy_notification(self):
        """ Send an optional notification for the processing completion. """
        jobs = []

        input = self.copy_inputs

        output1 = "notificationProcessingComplete." + str(self.lane_number) + ".out"
        output2 = "notificationCopyStart." + str(self.lane_number) + ".out"

        notification_command = config.param('start_copy_notification', 'notification_command', required=False)
        if (notification_command):
            job = Job([input], [output1, output2], name="start_copy_" + self.run_id + "_" + str(self.lane_number))
            job.command = notification_command.format(
                technology = config.param('start_copy_notification', 'technology'),
                output_dir = self.output_dir,
                run_id = self.run_id,
                output1 = output1,
                output2 = output2,
                lane_number = self.lane_number
            )
            jobs.append(job)        

        return jobs

    def copy(self):
        """Copy processed files to another place where they can be served or loaded into a LIMS."""
        jobs = []

        input = self.copy_inputs
        output = self.output_dir + "copyCompleted." + str(self.lane_number) + ".out"

        # TODO: the actual copy
        job = Job([input], [output], name="copy_" + self.run_id + "_" + str(self.lane_number), command="touch " + output)
        jobs.append(job)
        return jobs

    def end_copy_notification(self):
        """ Send an optional notification to notify that the copy is finished. """
        jobs = []

        input = self.output_dir + "copyCompleted." + str(self.lane_number) + ".out"
        output = self.output_dir + "notificationAssociation." + str(self.lane_number) + ".out"

        notification_command = config.param('end_copy_notification', 'notification_command', required=False)
        if (notification_command):
            job = Job([input], [output], name="end_copy_" + self.run_id + "_" + str(self.lane_number))
            job.command = notification_command.format(
                technology = config.param('end_copy_notification', 'technology'),
                output_dir = self.output_dir,
                run_name = os.path.basename(self.run_directory),
                output = output,
                lane_number = self.lane_number
            )
            jobs.append(job)

        return jobs

    def getSequencerIndexLength(self):
        """ Returns the total number of index cycles of the run. """
        return sum(index_read.nb_cycles for index_read in [read for read in self.read_infos if (read.is_index)])

    def validateBarcodes(self):
        """ Validate all index sequences against each other to ensure they aren't in collision according to the chosen number of mismatches parameter."""
        min_allowed_distance = (2 * self.number_of_mismatches) + 1

        validated_indexes = []
        collisions = []

        for readset in self.readsets:
            current_index = readset.index.replace('-', '')

            for candidate_index in validated_indexes:
                if (distance(current_index, candidate_index) < min_allowed_distance):
                    collisions.append("'" + current_index + "' and '" + candidate_index + "'")
            validated_indexes.append(current_index)

        if (len(collisions) > 0):
            raise Exception("Barcode collisions: " + ";".join(collisions));

    def get_mask(self):
        """ Returns a BCL2FASTQ friendly mask of the reads cycles.
            The mask is calculated using:
                - first base and last base of index;
                - the index length in the sample sheet;
                - the number of index cycles on the sequencer;
        """
        mask = ""
        index_lengths = self.get_smallest_index_length()
        index_read_count = 0
        nb_total_index_base_used = 0

        for read_info in self.read_infos:
            if (len(mask) > 0):
                mask += ','
            if (read_info.is_index):
                if (read_info.nb_cycles >= index_lengths[index_read_count]):
                    if (index_lengths[index_read_count] == 0 or self.last_index <= nb_total_index_base_used):
                        # Don't use any index bases for this read
                        mask += 'n' + str(read_info.nb_cycles)
                    else:
                        nb_n_printed = 0

                        # Ns in the beginning of the index read
                        if (self.first_index > (nb_total_index_base_used + 1)):
                            nb_n_printed = min(read_info.nb_cycles, self.first_index - nb_total_index_base_used - 1)
                            if (nb_n_printed >= index_lengths[index_read_count]):
                                nb_n_printed = read_info.nb_cycles
                            mask += 'n' + str(nb_n_printed)

                        # Calculate the number of index bases
                        nb_index_bases_used = max(index_lengths[index_read_count] - nb_n_printed, 0)
                        nb_index_bases_used = min(self.last_index - nb_total_index_base_used - nb_n_printed, nb_index_bases_used)
                        nb_total_index_base_used += nb_index_bases_used + min(nb_n_printed, index_lengths[index_read_count])
                        if (nb_index_bases_used > 0):
                            mask += 'I' + str(nb_index_bases_used)

                        # Ns at the end of the index read
                        remaining_base_count = read_info.nb_cycles - nb_index_bases_used - nb_n_printed
                        if (remaining_base_count > 0):
                            mask += 'n' + str(remaining_base_count)
                index_read_count += 1
            else:
                # Normal read
                mask += 'Y' + str(read_info.nb_cycles)
        return mask

    def generate_illumina_lane_sample_sheet(self):
        """ Create a sample sheet to use with the BCL2FASTQ software.
            Only the samples of the chosen lane will be in the file.
            The sample indexes are trimmed according to the mask used.
        """
        read_masks = self.get_mask().split(",")

        csv_headers = ["FCID", "Lane", "SampleID" , "SampleRef", "Index", "Description",
                       "Control", "Recipe", "Operator", "SampleProject"]
        has_single_index = self.has_single_index();

        for readset in self.readsets:
            index_to_use = ""

            if len(readset.index) > 0 and len(self.readsets) > 1:
                indexes = readset.index.split("-")
                nb_index = len(indexes)

                if has_single_index:
                    # we have a mixed of index in the sample, there are samples with 1 or 2 index,
                    # ignore the second index in the samplesheet
                    nb_index = 1

                for i in range(0, nb_index):
                    nb_ignored_leading_bases = 0
                    nb_of_index_bases = 0

                    m = re.match("(n\d+)?(I\d+)(n\d+)?", read_masks[i+1])
                    if m :
                        if m.group(1):
                            nb_ignored_leading_bases = m.group(1)[1:]
                        if m.group(2):
                            nb_of_index_bases = m.group(2)[1:]

                    # remove ignored leading bases and trim index to smallest lane index
                    index = indexes[i][nb_ignored_leading_bases:nb_ignored_leading_bases + nb_of_index_bases]

                    if i > 0 and len(index) > 0:
                        index_to_use += "-"
                    index_to_use += index

            csv_dict = {
                "FCID" : readset.flow_cell,
                "Lane" : self.lane_number,
                "SampleID" : readset.name,
                "SampleRef" : readset.species,
                "Index" : index_to_use,
                "Description" : readset.description,
                "Control" : readset.control,
                "Recipe" : readset.recipe,
                "Operator" : readset.operator,
                "SampleProject" : readset.project
            }

        csv_file = self.output_dir + os.sep + config.param('global', 'casava_sample_sheet_prefix') + str(self.lane_number) + ".csv"
        writer = csv.writer(open(csv_file, 'wb'),
                            delimiter=',',
                            quotechar="",
                            fieldnames=csv_headers)
        writer.writeheader()
        for key, value in csv_dict.items():
            writer.writrerow([key, value])

    def has_single_index(self):
        """ Returns True when there is at least one sample on the lane that doesn't use double-indexing. """
        return len([readset for readset in self.readset if ("-" not in readset.index)]) > 0

    def get_smallest_index_length(self):
        """ Returns a list (for each index read of the run) of the minimum between the number of index cycle on the sequencer and all the index lengths."""
        run_index_lengths = [r.nb_cycles for r in self.read_infos if r.is_index]

        if (len(run_index_lengths) == 0 and len(self.readsets) > 1):
            raise Exception("Multiple samples on a lane, but no indexes were read from the sequencer.")

        for i in range(0, len(run_index_lengths)):
            min_sample_index_length = min(len(readset.index.split("-")[i])
                                            for readset in
                                            self.readsets
                                            if (len(readset.index.split("-")) > i and len(readset.index.split("-")[i]) > 0)
                                          )
            run_index_lengths[i] = 0 if (min_sample_index_length is None) else min(min_sample_index_length, run_index_lengths[i])

        return run_index_lengths

    def parse_run_info_file(self):
        """ Parse the RunInfo.xml file of the run and load the reads data in the memory """
        reads = XML.parse(self.run_directory + os.sep + "RunInfo.xml").getroot().find('Run').find('Reads')
        return [RunInfoRead(int(r.get("Number")), int(r.get("NumCycles")), r.get("IsIndexedRead") == "Y") for r in reads.iter('Read')]

def distance(str1, str2):
    """ Returns the hamming distance. http://code.activestate.com/recipes/499304-hamming-distance/#c2 """
    return sum(itertools.imap(str.__ne__, str1, str2))


if __name__ == '__main__': 
    pipeline = IlluminaRunProcessing();
    pipeline.submit_jobs();