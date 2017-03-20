#!/usr/bin/env python
"""
Epigenetics pipeline for RRBS/WGBS data
=======================================

This simple pipeline uses Bismark, Trim Galore!, R, and Picard to process data. See the README file for more information
about the steps in the pipeline. See the accompanying documentation regarding tutorials and running the pipeline.

Many objects are defined by the superclass, note that the __init__ function calls on the init function of the parent.
Refer to the MUGQIC modules that are imported below. Looking at their source code will also help you understand some
of the underlying infrastructure that keeps this pipeline working.

Background
----------z
Epi-Seq is a differential analysis pipeline for BS-seq sequencing. Currently only RRBS and WGBS datasets are
tested to work with this pipeline. Similar to the other MUGQIC pipeline series, EpiSeq uses two metadata files to
set up the pipeline. The design file is used to group samples into case vs control. The readsets files maps input files
to names and, if applicable, which pairs of input files correspond to paired reads.
A readset is considered to be one instance/lane/segment of sequencing data. This is often used when multiple libraries
are used for a given sample, or if multiplexing was done. These techniques tend to generate multiple sets of data
for a given sample. The readset file allows users to specify these relationships.

The main outputs are the results of the differential analysis. In this pipeline, 2 different analysis are done. One is
for marking regions, with the other marks positions. These two files are in the form of a `.csv` file, which could be
parsed to produce helpful visualizations. However, the two files are unlikely to be sufficient in an experiment. Indeed,
some form of quality metrics needs to be collected throughout the pipeline to ensure the the final results are accurate.
This can be found in the reports at various steps and are listed [below](#output-reports).

Input
-----
- `FASTQ` or `BAM` files containing methylation sequencing data
- Reference genome in `FASTA` format
- MUGQIC formatted `.design` file
- MUGQIC formatted `.readset` file
- Epi-seq pipeline's `.ini` file

Output Data
-----------
- Bam file:
    - dedup/{sample.name}/{sample.name}.merged.deduplicated.bam
- Methylation Graph Plot:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated.bedGraph.gz
- Methylation Coverage:
    - methyl_calls{sample.name}/{sample.name}.merged.deduplicated.bismark.cov.gz
- CpG Methylation List:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated.CpG_report.txt.gz
- Differentially Methylated Positions:
    - differential_methylated_positions/{contrast.name}_RRBS_differential_methylated_pos.csv
- Differentially Methylated Regions:
    - differential_methylated_regions/{contrast.name}_RRBS_differential_methylated_regions.csv

Output Reports
--------------
- Pre-QC Metrics Report:
    - pre_qc_check/{sample.name}/{readset.name}/{readset.name}_fastqc.html (multiple readsets)
    - pre_qc_check/{sample.name}/{readset.name}_fastqc.html (1 readset in sample)
- Trimming Report:
    - trimmed/{sample.name}/{readset.name}/{readset.name}/{readset.bam}_trimming_report.txt (multiple readsets)
    - trimmed/{sample.name}/{readset.name}/{readset.bam}_trimming_report.txt (1 readset)
- Post-Trim QC Report:
    - trimmed/{sample.name}/{readset.name}/{readset.name}/{readset.bam}_fastqc.html (multiple readsets)
    - trimmed/{sample.name}/{readset.name}/{readset.bam}_fastqc.html (1 readset)
- Alignment Report:
    - merged/{sample.name}/{sample.name}.merged_aligned_PE_report.txt
- Deduplication Report:
    - dedup/{sample.name}/{sample.name}.merged.deduplication_report.txt
- Nucleotide Coverage:
    - dedup/{sample.name}/{sample.name}.merged.deduplicated.nucleotide_stats.txt
    - merged/{sample.name}/{sample.name}.merged.nucleotide_stats.txt
- Methylation Extraction Report:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated_splitting_report.txt
- M-bias Plot:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated.M-bias.txt
- Bismark HTML Sample Summary:
    - bismark_summary_report/{sample.name}_final_bismark_report.html

Workflow
--------
Below is the steps defined in this pipeline. The alphabetical designation shows which steps run concurrently in the
pipeline. It shows the order of the pipeline and the dependencies for each step.

 1. (a) Genome Methylation Conversion (bismark_prepare_genome)
 2. (a) Pre-Trim Quality Check (pre_qc_check)
 3. (a) Read and Adapter Trimming (trim_galore)
 4. (b) Read Alignment (bismark_align)
 5. (c) Merge Alignment Statistics (merge_bismark_alignment_report)
 6. (c) Merge Aligned BAMs (picard_merge_sam_files)
 7. (d) Recalculate nucleotide coverage (merged_nuc_stats)
 8. (d) Remove Duplicate Reads (bismark_deduplication)
 9. (e) Recalculate nucleotide coverage 2 (calc_dedup_nucleotide_coverage)
10. (e) Methylation Calling and Analysis (bismark_methylation_caller)
11. (f) Bismark Sample-Level Report Generator (bismark_html_report_generator)
12. (g) Position Specific Differential Analysis (differential_methylated_pos)
13. (g) Regional Differential Analysis (differential_methylated_regions)

Creating new Steps and Jobs
---------------------------------------
The pipeline is arranged to gradually decrease the number of input files when it make sense to do so. The trimming
and QC step would be too slow if we merge readsets together, so we decided to trim each file individually. Then,
the alignment phase allows us to convert out fastQ files to BAM files. This allows us to merge paired datasets into
one BAM file per readset without explicitly doing so. After this, readsets are combined to run the analysis steps.

Every step in the pipeline is simply a method within the Episeq class. These methods all have the following
signature: (None) -> List(Job). In general, every step should first look at the following.
This structure/framework is something that should be constant between each job declaration of a step.
1. What readsets/samples are involved
2. Additional user configuration from .ini file
3. Required environment variables
4. Required resources

Next, jobs may require some additional information. These points need to be re-evaluated between samples/jobs:
1. Determine which data is eligible to run
2. Determine input and output file paths.
3. Determine what formatting should be done prior to running (creating directories, etc.)

The above highlights the many factors when creating a new step in the pipeline. There is a fair amount of
preparation required to obtain every bit of information. In addition, these information may require different
sources. For example, readset and sample information comes automatically from the superclasses (which is from
the readset file as `self.readsets`), but the path information needs to be referenced from other jobs within the
pipeline. Input files would need to be obtained from the previous step, which will require recreating this
information. Keep in mind that any fixed strings may require refactoring in other related steps.

With the above information generated and gathered, creating Job objects should be straightforward. Filling all params
in a Job object will help prompt you for the correct information. Please look at the examples below to help identify
tricks to format a bash script or to join separate jobs into one job object.

Notes about the final report
----------------------------
Developer's note on generating report stubs. The final report document requires a small bit of data from each job that
ran. Unfortunately, not all jobs complete and jobs don't necessarily run in the same order, so idempotency is a
valuable property to have in the steps below. You will find that the following 3 bits are important:
1. Identify what output files are needed to display the final results. (Result representation)
2. Identify how to organize all jobs onto a report page. (Template creation)
3. Identify what operations need to be done to visualize the results. (Input formatting)

Code Organization
-----------------
A simpler system is required to avoid repetitive code that will easily cause errors when refactoring. Otherwise, the
code will be large and difficult to maintain. As an effort to maintain readability in the code, here are some metrics
produce by the python package radon as of Dec. 8/16:

Maintainability Index

episeq.py - A (29.92)  - including doc strings

episeq.py - B (16.41)  - without doc strings

Cyclomatic Complexity score

episeq.py

    M 316:4 Episeq.trim_galore - D (26)
    M 531:4 Episeq.bismark_align - D (23)
    M 751:4 Episeq.picard_merge_sam_files - C (15)
    M 210:4 Episeq.pre_qc_check - C (14)
    M 1090:4 Episeq.differential_methylated_pos - C (12)
    M 991:4 Episeq.bismark_html_report_generator - B (8)
    C 154:0 Episeq - B (7)
    M 710:4 Episeq.merge_bismark_alignment_report - B (7)
    M 1198:4 Episeq.differential_methylated_regions - B (7)
    M 848:4 Episeq.bismark_deduplicate - A (5)
    M 933:4 Episeq.bismark_methylation_caller - A (4)
    M 828:4 Episeq.merged_nuc_stats - A (3)
    M 913:4 Episeq.calc_dedup_nucleotide_coverage - A (3)
    M 160:4 Episeq.__init__ - A (1)
    M 164:4 Episeq.merge_py - A (1)
    M 172:4 Episeq.bismark_prepare_genome - A (1)
    M 1296:4 Episeq.bam2nuc_job - A (1)
    M 1323:4 Episeq.steps - A (1)

"""
# Python Standard Modules
import os

# MUGQIC Modules
from pipelines.common import Illumina, Job, concat_jobs, config, logging
from bfx import metrics

# Use this logger to print warning messages to the debug log. (Global, imported from common.py)
log = logging.getLogger(__name__)


# Locally defined. May be moved to "bfx" package.
def bam2nuc_job(output_dir, sample_name, suffix, in_bam):
    """
    Generates jobs for Bismark's bam2nuc script. Used to calculate nucleotide coverage statistics.

    :param output_dir: A specified output directory for the report file
    :type output_dir: str
    :param sample_name: The sample of the sample to run
    :type sample_name: str
    :param suffix: A suffix to add to the filename before '.nucleotide_stats.txt'
    :type suffix: str
    :param in_bam: The bam file to analyse.
    :type in_bam: str
    :return: A Job object to run bam2nuc
    :rtype: Job
    """
    output_file = os.path.join(output_dir, sample_name + suffix + '.nucleotide_stats.txt')
    coverage_calc = Job(
        input_files=[in_bam],
        output_files=[output_file],
        module_entries=[['bismark_deduplicate', 'module_samtools'],
                        ['bismark_deduplicate', 'module_perl'],
                        ['bismark_deduplicate', 'module_bismark']],
        command='bam2nuc --dir ' + output_dir + ' --genome_folder bismark_prepare_genome ' + in_bam,
        removable_files=[output_file],
        name='bam2nuc.' + sample_name)
    return coverage_calc


class EpiSeq(Illumina):
    """
    The Episeq pipeline takes FASTQ or BAM files (unsorted) as input as well as two metadata files and a configuration
    file. Refer to the user guide for more information on running the pipeline.
    """

    def __init__(self):
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(EpiSeq, self).__init__()

    @property
    def merge_py(self):
        """
        :return: A full path to the bismark report merge script, for easy reference
        :rtype: str
        """
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bismark_merge_reports.py')

    @property
    def steps(self):
        """
        These are the steps to the pipeline and are listed in the order that should be processed. Don't forget to
        change this as needed!

        :return: list
        """
        return [
            self.bismark_prepare_genome,
            self.pre_qc_check,
            self.trim_galore,
            self.bismark_align,
            self.merge_bismark_alignment_report,
            self.picard_merge_sam_files,
            self.merged_nuc_stats,
            self.bismark_deduplicate,
            self.calc_dedup_nucleotide_coverage,
            self.bismark_methylation_caller,
            self.bismark_html_report_generator,
            self.methylation_values,
            self.differential_methylated_pos,
            self.differential_methylated_regions,
            self.prepare_annotations,
            self.annotate_regions,
            self.annotate_positions
        ]

    # Pipeline steps start here
    def bismark_prepare_genome(self):  # Step 1
        """
        Bismark requires a processed reference genome to compare with the epigenome. This step can take several hours,
        depending on the size of the genome. It creates an index of bisulfite conversions and takes make space than
        the genome itself. This module will only create the output in the same directory as the genome file, so
        a symlink is needed to create a "copy" in the desired output directory. Only runs once and generates one job.

        Note: The genome should be in fasta format.

        Note: main_job can be modified to run bam2nuc, but there's a bug in Bismark that breaks it, fix when possible.

        Input: A reference sequence file as specified by user. Configuration is set in the episeq.ini file.
        Output: A directory called Bisulfite_Genome and a dinucleotide composition report.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Variable preparations
        ref_seq = config.param('bismark_prepare_genome', 'genome_file', type='filepath')

        local_ref_seq = os.path.join('bismark_prepare_genome', os.path.basename(ref_seq))
        output_idx = "bismark_prepare_genome/Bisulfite_Genome"
        modules = [['bismark_prepare_genome', 'module_bowtie2'],
                   ['bismark_prepare_genome', 'module_samtools'],
                   ['bismark_prepare_genome', 'module_perl'],
                   ['bismark_prepare_genome', 'module_bismark']]

        # Job creation - If ref_seq appears to be processed already, just make a symlink
        if os.path.isdir(os.path.join(os.path.dirname(ref_seq), 'Bisulfite_Genome')):
            log.info("Found converted genome from given genome file. Creating symlink instead.")
            job = Job(input_files=[ref_seq],
                      output_files=[output_idx],
                      removable_files=[output_idx],
                      command='cp -sL ' + os.path.dirname(ref_seq) + ' ' + 'bismark_prepare_genome')
        else:  # Process genome file
            log.info("Unable to find any signs of preprocessing. Will prepare genome for pipeline.")
            mkdir_job = Job(command='mkdir -p bismark_prepare_genome')
            link_job = Job(input_files=[ref_seq],
                           command='cp -sfu ' + os.path.abspath(ref_seq) + ' ' +
                                   os.path.abspath(os.path.join(self.output_dir, local_ref_seq)),
                           removable_files=[local_ref_seq])
            main_job = Job(output_files=[output_idx], module_entries=modules,
                           command="bismark_genome_preparation --verbose bismark_prepare_genome/",
                           removable_files=[output_idx])
            nuc_count = Job(output_files=[], module_entries=modules,
                            command="bam2nuc --genomic_composition_only --genome_folder bismark_prepare_genome/ "
                                    "--dir bismark_prepare_genome/",
                            removable_files=[])
            job = concat_jobs([mkdir_job, link_job, main_job, nuc_count],
                              name='bismark_prepare_genome.' + os.path.basename(ref_seq))
        return [job]

    def pre_qc_check(self):  # Step 2
        """
        Runs FastQC on the unprocessed fastq files to generate a baseline report. Helpful when comparing to post-trim
        metrics.

        Note: While FastQC does use perl, do not load any perl modules until the sha-bang statement uses /usr/bin/env.
        Failing to do so will result in incorrect libraries being used as /usr/bin/perl is not likely the same version
        as what you load.

        A single report file is generated with the help of a intermediate file
        that lets tasks update the table format as they complete. The drawback
        is that pandoc will have to run for each task, but is to allow some
        jobs to fail and still have a working report.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []  # List of jobs to run
        requires = [['pre_qc_check', 'module_java'],
                    ['pre_qc_check', 'module_fastqc']]
        tmpdir = config.param('pre_qc_check', 'temp_dir', required=False) or config.param('DEFAULT', 'tmp_dir')
        # Must set file name in the bash script and not in python, else the job is not up to date just because of day change
        #template_var_hold = 'pre_qc_check/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        report_file = 'report/EpiSeq.pre_qc_check.md'
        report_data = 'data/pre_qc_check'
        fill_in_templ = "| {sample} | {readset} | [Report 1]({fq1_report})<br>[Download 1]({fq1_download}) | " \
                        "{fq2_report}<br>{fq2_download} | {start} | {completion} |"
        # Generate jobs
        for sample in self.samples:
            for readset in sample.readsets:
                # Check inputs for readset
                raw_fq = filter(None, [readset.fastq1, readset.fastq2])  # size: [0,2]
                if not raw_fq:
                    log.info('No fastq files found for readset for ' + readset.name + '. Skipping...')
                    if readset.bam:  # Check if we have a bam file
                        continue
                    else:
                        raise ValueError('No sequencing data found for ' + readset.name)

                # Determine output directory format
                # My method of collapsing unneeded directories. If we only have 1 readset, don't make nested directories
                if len(sample.readsets) == 1:
                    out_dir = os.path.join('pre_qc_check', sample.name)
                else:
                    out_dir = os.path.join('pre_qc_check', sample.name, readset.name)

                # Generate output path
                # Obtain the output of the filename. ex. 'SRS23542_1' + '...'
                id_name = [os.path.basename(nom).split('.gz')[0].split('.fastq')[0] for nom in raw_fq]  # Basename
                file_names = [nom + '_fastqc' + ext for nom in id_name for ext in ['.html', '.zip']]  # Output Name
                output = [os.path.join(out_dir, name) for name in file_names]  # Output file path

                # Job creation
                mkdir_job = Job(command='mkdir -p ' + out_dir)
                job = Job(input_files=raw_fq,
                          output_files=output,
                          module_entries=requires,
                          command="fastqc -o {out_dir} {others} -d {tmpdir} {inputs}".format(
                              inputs=' '.join(raw_fq),
                              out_dir=out_dir,
                              others=config.param('pre_qc_check', 'other_options', required=False),
                              tmpdir=tmpdir),
                          removable_files=[output])

                # Fill in the blanks to generate a row in the report table. Will list the available jobs and results
                report_body = fill_in_templ.format(
                    sample=sample.name,
                    readset=readset.name,
                    fq1_report=os.path.join(report_data, id_name[0]) + '_fastqc.html',
                    fq1_download=os.path.join(report_data, id_name[0] + '.zip'),
                    fq2_report='[Report 2](' + os.path.join(report_data, id_name[1]) + '_fastqc.html' + ')'
                               if len(id_name) == 2 else '',
                    fq2_download='[Download 2](' + os.path.join(report_data, id_name[1] + '.zip') + ')'
                                 if len(id_name) == 2 else '',
                    start="$START",
                    completion="$(date \"+%Y-%m-%d %H:%M:%S\")"
                )
                # Code to fill in template.
                command = """\
TEMPLATE_STR_FILE=pre_qc_check/$(date +%F)_template_var_strings.txt && \\
flock -x "${{TEMPLATE_STR_FILE}}.lock" -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
mkdir -p {data_loc} && \\
cp -f {output_file} {data_loc} && \\
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable data_table="$table" \\
  --to markdown \\
  > {report_file}""".format(
                    entry=report_body,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    output_file=" ".join(output),
                    data_loc=os.path.join('report', report_data),
                    report_file=report_file)
                update_template = Job(
                    output_files=[os.path.join('report', report_data, os.path.basename(out_log))
                                  for out_log in output] + [report_file],
                    module_entries=[['pre_qc_check', 'module_pandoc']],
                    command=command,
                    report_files=[report_file])

                # Add to list of jobs
                jobs.append(concat_jobs([Job(command="START=$(date '+%Y-%m-%d %H:%M:%S')"),  # Get time when running
                                         mkdir_job, job, update_template], name='pre_qc_check.' + readset.name))
        return jobs

    def trim_galore(self):  # Step 3
        """
        This step trims raw FASTQ files for quality control using Trim Galore!
        This is a pre-processing step to ensure quality control.

        To run Trim Galore in paired mode, two fastq files must be specified in the readset file with the following
        pairewise naming convention: file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz

        Note: Do not load perl modules because Trim Galore and FastQC only uses /usr/bin/perl. It will not use the
        correct perl binary and thus the libraries are set up for the wrong version.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Prep
        report_file = 'report/EpiSeq.trim_galore.md'
        report_data = 'data/trim_galore'
        #template_string_file = 'trimmed/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {readset} | {trim1_view} {trim1_download}<br>{trim2_view} {trim2_download} | ' \
                        '{qc1_view} {qc1_download}<br>{qc2_view} {qc2_download} | {start} | {completion} |'
        # Awk script to parse trim reports view easy viewing
        # It will read the start of each line to determine if it is something that we want to keep, or discard.
        # Headers have to be reformatted to a lower level, making the document make sense within the template.
        # Additional code is added to properly add table syntax when extra blank lines are found. This should create
        # *.md files that we can link to from the pipeline report.
        output_parser = """\
BEGIN {table=0;}
{
    if ($0 ~ /^SUMMARISING/) {
        print "##### " $0;
        table=0;
    } else if ($0 ~ /^RUN STATISTICS/) {
        print "###### " $0;
        table=0;
    } else if ($0 ~ /^=+$/) {
        table=0;
        next;
    } else if (match($0, /^=== ([[:alnum:][:blank:]]+) ===$/, a)) {
        print "###### " a[1];
        table=0;
    } else if ($0 ~ /^\|Input filename.+/) {
        print ""; print "|Attribute|Value|";
        print "|----|----|"; print $0;
        table=1;
    } else if ($0 ~ /^\|Total reads.+/) {
        print "";
        print "|Metric|Value|";
        print "|----|----|"; print $0;
        table=1;
    } else if ($0 ~ /^length\|.+/) {
        print "";
        print $0;
        print "----:|----:|----:|----:|----:";
        table=1;
    } else if ($0 ~ /^Bases preceding.+/) {
        print $0;
        print "";
        print "|Base|Value|";
        print "|:----:|----:|";
        table=1;
    } else if ($0 ~ /^Running/ || $0 ~/^Output/ || $0 ~ /^This is cut.+/) {
        next;
    } else if (/^\|/) {
        print $0;
        table=1;
    } else if ($0 ~ /^$/ && table == 1) {
        next;
    } else if ($0 ~ /^$/ && table == 0) {
        print $0;
    } else {
        print $0;
    }
}"""
        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:  # Iterate through each readset in project
                if readset.bam:  # Input check!
                    log.info('User-supplied bam file found for ' + readset.name + '. Skipping...')
                    continue

                # Setup output directory structure
                trim_directory = os.path.join("trimmed", sample.name)
                if len(sample.readsets) > 1:
                    trim_directory = os.path.join(trim_directory, readset.name)

                # Some variables about the readset and current settings
                output_reports = list()
                key_report = list()
                run_type = readset.run_type
                protocol = readset.library
                input_files = filter(None, [readset.fastq1, readset.fastq2])
                file_basename = [os.path.join(trim_directory, os.path.basename(in_file).split('.')[0])
                                 for in_file in input_files]  # Might be a bit too agressive of a split.
                add_options = config.param('trim_galore', 'other_options').split()  # options that may affect output
                run_qc = '--fastqc_args' in add_options or '--fastqc' in add_options
                report_out = '--no_report_file' not in add_options

                # List all possible outputs
                if len(input_files) == 2:
                    input1_logs = [trim_directory + '/' + os.path.basename(readset.fastq1) + '_trimming_report.txt',
                                   file_basename[0] + "_val_1_fastqc.html",
                                   file_basename[0] + "_val_1_fastqc.zip",
                                   file_basename[0] + "_val_1.fq.gz"]
                    input2_logs = [trim_directory + '/' + os.path.basename(readset.fastq2) + '_trimming_report.txt',
                                   file_basename[1] + "_val_2_fastqc.html",
                                   file_basename[1] + "_val_2_fastqc.zip",
                                   file_basename[1] + "_val_2.fq.gz"]
                else:  # Different names, so have to have this block of code
                    input1_logs = [trim_directory + '/' + os.path.basename(readset.fastq1) + '_trimming_report.txt',
                                   file_basename[0] + "_trimmed_fastqc.html",
                                   file_basename[0] + "_trimmed_fastqc.zip",
                                   file_basename[0] + "_trimmed.fq.gz"]

                # Build list of input reports/data
                try:
                    if run_type == "PAIRED_END":
                        if report_out:
                            output_reports += [input1_logs[0]] + [input2_logs[0]]  # Trim report
                            key_report = [input1_logs[0]] + [input2_logs[0]]
                        if run_qc:
                            output_reports += [input1_logs[1]] + [input2_logs[1]]  # FastQC html
                            output_reports += [input1_logs[2]] + [input2_logs[2]]  # FastQC zip
                        output_files = [input1_logs[3]] + [input2_logs[3]]
                    else:
                        if report_out:
                            output_reports += [input1_logs[0]]
                            key_report = [input1_logs[0]]
                        if run_qc:
                            output_reports += [input1_logs[1]] + [input1_logs[2]]
                        output_files = [input1_logs[3]]
                except NameError:  # When fastq2 is None
                    log.error("Second fastq file is missing. Expecting files as follows file1_1.fq, file1_2.fq.")
                    log.error("Files in readset " + readset.name + " :" + " ".join(input_files))
                    raise

                # Define jobs
                mkdir_job = Job(command="mkdir -p " + trim_directory)
                trim_job = Job(input_files=input_files,
                               output_files=output_reports + output_files,
                               module_entries=[['trim_galore', 'module_fastqc'],
                                               ['trim_galore', 'module_java'],
                                               ['trim_galore', 'module_trim_galore'],
                                               ['trim_galore', 'module_cutadapt']],
                               command="trim_galore {protocol} {library_type} {other} {directory} {fastq}".format(
                                   library_type="--paired" if run_type == "PAIRED_END" else "",
                                   protocol='--rrbs' if protocol == 'RRBS' else '',
                                   other=config.param("trim_galore", "other_options"),
                                   directory='--output_dir ' + trim_directory,
                                   fastq=' '.join(input_files)),
                               removable_files=output_files)

                entry = fill_in_entry.format(
                    sample=sample.name,  # For report jobs.
                    readset=readset.name,
                    trim1_view='[View details 1](' +
                               os.path.join(report_data, os.path.basename(input1_logs[0])) + '.html)',
                    trim1_download='[Download details 1](' +
                                   os.path.join(report_data, os.path.basename(input1_logs[0])) + ')',
                    trim2_view='[View details 2](' +
                               os.path.join(report_data, os.path.basename(input2_logs[0])) + '.html)'
                    if input2_logs else '',
                    trim2_download='[Download details 2](' +
                                   os.path.join(report_data, os.path.basename(input2_logs[0])) + ')'
                    if input2_logs else '',
                    qc1_view='[View report 1](' +
                             os.path.join(report_data, os.path.basename(input1_logs[1])) + ')',
                    qc1_download='[Download report 1](' +
                                 os.path.join(report_data, os.path.basename(input1_logs[2])) + ')',
                    qc2_view='[View report 2](' +
                             os.path.join(report_data, os.path.basename(input2_logs[1])) + ')'
                    if input2_logs else '',
                    qc2_download='[Download report 2](' +
                                 os.path.join(report_data, os.path.basename(input2_logs[2])) + ')'
                    if input2_logs else '',
                    start="$START",
                    completion="$(date '+%Y-%m-%d %H:%M:%S')")
                command = """\
TEMPLATE_STR_FILE=trimmed/$(date +%F)_template_var_strings.txt && \\
flock -x "${{TEMPLATE_STR_FILE}}.lock" -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
mkdir -p {data_loc} && \\
cp -f {output_file} {data_loc}; \\
for i in {individual_page}; do
    sed -r 's%^([^0-9S][a-z ].+):(\s+.+)%|\\1|\\2|%g' {data_loc}/$i | \\
        sed 's/\t/|/g' | \\
        awk '{script}' | \\
        pandoc --output "{data_loc}/$i.html"; \\
done
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
{report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown \\
> {report_file}""".format(
                    entry=entry,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    output_file=" ".join(output_reports),
                    individual_page=" ".join([os.path.basename(loc) for loc in key_report]),
                    data_loc=os.path.join('report', report_data),
                    script=output_parser,
                    report_file=report_file)

                report_job = Job(
                    # removed key_report files from output_files as they never get written or used by other jobs
                    output_files=[report_file] +
                                 [os.path.join('report', report_data,
                                               os.path.basename(item)) for item in output_reports],
                    command=command,
                    module_entries=[['trim_galore', 'module_pandoc']],
                    report_files=[report_file])

                jobs.append(concat_jobs([Job(command="START=$(date '+%Y-%m-%d %H:%M:%S')"),
                                         mkdir_job, trim_job, report_job], name='trim_galore.' + readset.name))

        return jobs

    def bismark_align(self):  # Step 4
        """
        This step aligns trimmed reads to a bisulfite converted reference genome using Bismark. This create
        BAM files and will only be compressed if the input is also compressed (which usually will be the case).
        All readsets are aligned individually, but combines paired files into one output file. The output files
        are all placed in the same sample directory. No sub-directories are made.

        This step requires bismark_prepare_genome and the relevant trim_galore step.

        Note: Despite what the manual says, the source code shows that -un and --ambiguous produces fq files, not txt.

        Input: Trimmed version of input files as a fastq file. (trimmed/*)
        Output: A BAM/SAM file in aligned/<sample_name>/*

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Generation Setup
        report_file = 'report/EpiSeq.bismark_align.md'
        report_data = 'data/bismark_align'
        #template_string_file = 'aligned/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {readset} | {report_view} | {coverage_view} | {start} | {completion} |'
        awk_script = """\
BEGIN { table=0; };
{
    if ($0 ~ /^$/) {
        if (table==1) {
            next;}
        else {
            print $0;
        }
    }
    else if ($0 ~ /^[^|]+|[^|]+$/ && table==0) {
            print ""; print "|Item|Value|"; print "|-----|----:|"; print $0;
            table=1;
    }
    else if($0 ~ /^Final.+/) {
            table=0;
            print "";
            print "##### " $0;
    }
    else if($0 ~ /^Number of sequence.+:$/) {
        table=0''
        print "";
        print $0;
    }
    else {
        print $0;
    }
}"""
        jobs = []
        for sample in self.samples:  # Process samples
            for readset in sample.readsets:
                # Check if the readset has an additional level of directories
                trim_prefix = os.path.join("trimmed", sample.name)
                if len(sample.readsets) > 1:
                    trim_prefix = os.path.join(trim_prefix, readset.name)

                # Common Variables
                run_type = readset.run_type
                align_directory = os.path.join("aligned", sample.name)
                output_basename = os.path.join(align_directory, readset.name)
                user_options = config.param('bismark_align', 'other_options').split()
                input_basename = [os.path.join(trim_prefix, os.path.basename(in_file).split('.')[0])
                                  for in_file in filter(None, [readset.fastq1, readset.fastq2])]

                # Check what input files are found
                if not input_basename:
                    if readset.bam:
                        continue
                    else:
                        raise IOError("""{readset} has no input files!""".format(readset=readset.name))

                # Again, the suffix is hardcoded into the script. So we have to match it too. PE and SE have diff names
                if run_type == "PAIRED_END" and readset.fastq2:
                    input_files = [input_basename[0] + '_val_1.fq.gz', input_basename[1] + "_val_2.fq.gz"]
                    cmd_in = '-1 {fastq1} -2 {fastq2}'.format(fastq1=input_files[0], fastq2=input_files[1])
                    out_files = [output_basename + "_aligned_pe.bam"]
                    report_log = [output_basename + "_aligned_PE_report.txt"]
                    # Optional output files, depending on flags specified
                    if '--nucleotide_coverage' in user_options:
                        report_log += [output_basename + "_aligned_pe.nucleotide_stats.txt"]
                    if '-un' in user_options or '--unmapped' in user_options:
                        out_files += [output_basename + "_aligned_unmapped_reads_1.fq.gz",
                                      output_basename + "_aligned_unmapped_reads_2.fq.gz"]
                    if '--ambiguous' in user_options:
                        out_files += [output_basename + "_aligned_ambiguous_reads_1.fq.gz",
                                      output_basename + "_aligned_ambiguous_reads_2.fq.gz"]
                    if '--ambig_bam' in user_options:
                        out_files += [output_basename + '_aligned_pe.ambig.bam']
                elif run_type == "SINGLE_END":
                    input_files = [input_basename[0] + "_trimmed.fq.gz"]
                    cmd_in = '--single_end {fastq1}'.format(fastq1=input_files[0])
                    out_files = [output_basename + "_aligned.bam"]
                    report_log = [output_basename + "_aligned_SE_report.txt"]
                    # Optional output files - depends on flag specified.
                    if '--nucleotide_coverage' in user_options:
                        report_log += [output_basename + "_aligned.nucleotide_stats.txt"]
                    if '-un' in user_options or '--unmapped' in user_options:
                        out_files += [output_basename + "_unmapped_reads.fq.gz"]
                    if '--ambiguous' in user_options:
                        out_files += [output_basename + "_ambiguous_reads.fq.gz"]
                    if '--ambig_bam' in user_options:
                        out_files += [output_basename + '_aligned.ambig.bam']
                else:
                    raise AttributeError("Unknown run_type or unknown file output name for " + sample.name)

                # Job creation
                mkdir_job = Job(command="mkdir -p " + align_directory)
                job = Job(
                    input_files + ["bismark_prepare_genome/Bisulfite_Genome"],
                    out_files + report_log,
                    [["bismark_align", "module_bowtie2"],
                     ["bismark_align", "module_samtools"],
                     ['bismark_align', 'module_perl'],
                     ['bismark_align', 'module_bismark']],
                    command="""\
bismark -q {other} --temp_dir {tmpdir} --output_dir {directory} \
    --basename {basename} --genome_folder bismark_prepare_genome {input}""".format(
                        directory=align_directory,
                        other=config.param("bismark_align", "other_options"),
                        tmpdir=config.param('bismark_align', 'tmp_dir', required=False) or
                            config.param('DEFAULT', 'tmp_dir', required='True'),
                        input=cmd_in,
                        basename=readset.name + '_aligned'),
                    removable_files=out_files
                )

                # Generate report stub
                new_logs = [os.path.join('report', report_data, os.path.basename(txt)) for txt in report_log]
                report_entry = fill_in_entry.format(
                    sample=sample.name,
                    readset=readset.name,
                    report_view='[View HTML]({a})<br>[Download Raw]({b})'.format(
                        a=os.path.join(report_data, os.path.basename(report_log[0]) + '.html'),
                        b=os.path.join(report_data, os.path.basename(report_log[0]))),
                    coverage_view='[View HTML]({a})<br>[Download Text]({b})'.format(
                        a=os.path.join(report_data, os.path.basename(report_log[1]) + '.html'),
                        b="(" + os.path.join(report_data, os.path.basename(report_log[1])) + ")"
                          if '--nucleotide_coverage' in user_options else ""),
                    start="$START",
                    completion="$(date '+%Y-%m-%d %H:%M:%S')")
                command = """\
TEMPLATE_STR_FILE=aligned/$(date +%F)_template_var_strings.txt && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
mkdir -p {data_dir} && \\
cp -f {reports} {data_dir}; \\
for i in {logs}; do
    sed -r 's/:\s+/|/g' $i | egrep -v "^Option" | egrep -v "^=+$" | awk '{script}' | pandoc --output $i".html"; \\
done
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
{report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown \\
> {report_file}""".format(
                    entry=report_entry,
                    data_dir=os.path.join('report', report_data),
                    reports=' '.join(report_log),
                    logs=' '.join(new_logs),
                    script=awk_script,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                )
                report_job = Job(output_files=[report_file] + new_logs,
                                 report_files=[report_file],
                                 module_entries=[['bismark_align', 'module_pandoc']],
                                 command=command)

                jobs.append(concat_jobs([Job(command="START=$(date '+%Y-%m-%d %H:%M:%S')"),
                                         mkdir_job, job, report_job], name="bismark_align." + readset.name))
                # To next readset
        return jobs

    def merge_bismark_alignment_report(self):  # Step 5
        """
        This steps takes all of Bismark's alignment reports for a sample and merges them with a custom script.
        Some stats are recalculated to match the total read population and specific settings are lost due to the
        aggregation. Outputs the file to the merge directory, which will contain other merged results and outputs.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            log_reports = []
            align_directory = os.path.join("aligned", sample.name)  # Previous step's output dir

            # Get all log reports for this sample
            for readset in sample.readsets:
                if readset.bam or not (readset.fastq1 or readset.fastq2):
                    continue
                # Get correct name values
                log_basename = os.path.join(align_directory, readset.name)
                if readset.run_type == "PAIRED_END":
                    log_reports.append(log_basename + "_aligned_PE_report.txt")
                    output_report = os.path.join("merged", sample.name, sample.name + ".merged_aligned_PE_report.txt")
                else:
                    log_reports.append(log_basename + "_aligned_SE_report.txt")
                    output_report = os.path.join("merged", sample.name, sample.name + ".merged_aligned_SE_report.txt")

            # Job creation
            mkdir_job = Job(command="mkdir -p merged/" + sample.name)
            merge_job = Job(log_reports, [output_report],
                            [['merge_bismark_alignment_report', 'module_python']],
                            command="""python {script_loc} -o {output} -n {name} {logs}""".format(
                                script_loc=self.merge_py,
                                output=output_report,
                                name=sample.name,
                                logs=' '.join(log_reports)))

            job = concat_jobs([mkdir_job, merge_job], name="merge_align_reports." + sample.name)
            jobs.append(job)
        return jobs

    def picard_merge_sam_files(self):  # Step 6
        """
        This step merges all readsets of each sample into one handy bam file. Here, if a readset is defined by a
        bam file, it will finally be used to add with other readsets. Because merging multiple alignments together
        can dramatically change the coverage content, it is recalculated to reflect the combined reads.

        While already defined in the bfx module, I want to use Picard Tools v2.0.1. This (and newer versions) has
        a different syntax than before, requiring me to rewrite the job definition.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []
        for sample in self.samples:
            # Bam files from pipeline (FASTQ)
            processed_fastq_pe = [os.path.join('aligned', sample.name, readset.name + "_aligned_pe.bam") for
                                  readset in sample.readsets if readset.run_type == 'PAIRED_END' and not readset.bam]
            processed_fastq_se = [os.path.join('aligned', sample.name, readset.name + "_aligned.bam") for
                                  readset in sample.readsets if readset.run_type == 'SINGLE_END' and not readset.bam]
            processed_fastq = processed_fastq_pe + processed_fastq_se
            # Bam files from user, if specified by readset file. Not exclusive with having fastq
            listed_bam_files = [readset.bam for readset in sample.readsets if readset.bam != '']

            input_files = filter(None, processed_fastq + listed_bam_files)  # All bam files that belong to the sample
            merge_prefix = os.path.join('merged', sample.name)
            output_bam = os.path.join(merge_prefix, sample.name + '.merged.bam')
            mkdir_job = Job(command='mkdir -p ' + merge_prefix)

            # I want to use Picard Tools v2.0.1, which has a different syntax than v1.x
            if len(input_files) > 1:
                picard_v2 = Job(
                    input_files,
                    [output_bam],
                    [
                        ['picard_merge_sam_files', 'module_java'],
                        ['picard_merge_sam_files', 'module_picard']
                    ],
                    command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} \\
  -jar $PICARD_HOME/picard.jar MergeSamFiles \\
  VALIDATION_STRINGENCY=SILENT \\
  TMP_DIR={tmp_dir} \\
  {inputs} \\
  OUTPUT={output} \\
  USE_THREADING=true \\
  SORT_ORDER=queryname \\
  MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
                        tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
                        java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
                        ram=config.param('picard_merge_sam_files', 'ram'),
                        inputs=" \\\n  ".join(["INPUT=" + in_put for in_put in input_files]),
                        output=output_bam,
                        max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', type='int')))
                job = concat_jobs([mkdir_job, picard_v2], name="picard_merge_sam_files." + sample.name)

            elif len(input_files) == 1:  # Save time and resources by just copying the single data source
                input_nuc_stats = os.path.join('aligned', sample.name, sample.readsets[0].name)
                if processed_fastq_pe:
                    input_nuc_stats += '_aligned_PE_report.txt'
                else:
                    input_nuc_stats += '_aligned_SE_report.txt'
                target_readset_bam = input_files[0]
                job = concat_jobs([
                    mkdir_job,
                    Job([input_files[0]], [output_bam],
                        command="cp -s -L -f " + os.path.join(os.getcwd(), 'output', target_readset_bam) +
                                " " + os.path.join(os.getcwd(), 'output', output_bam),
                        removable_files=[output_bam])],
                    name="symlink_readset_sample_bam." + sample.name)
            else:
                raise ValueError('Sample ' + sample.name + ' has no readsets!')
            jobs.append(job)
        return jobs

    def merged_nuc_stats(self):  # Step 7
        """
        This step calculates the new nucleotide coverge for the merged bam file.
        This independent step is to reduce the likelihood of failing mid step.

        The output of this step is found in the relevant merge sub-directory.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            output_dir = os.path.join('merged', sample.name)
            if os.path.exists(os.path.join(output_dir, sample.name) + '.merged.nucleotide_stats.txt'):
                continue
            merged_bam = os.path.join(output_dir, sample.name) + '.merged.bam'
            job = bam2nuc_job(output_dir, sample.name, '.merged', merged_bam)
            jobs.append(job)
        return jobs

    def bismark_deduplicate(self):  # Step 8
        """
        Calls the de-duplication module from Bismark. The module operates on the output alignment files from
        Bismark, creating a new output file (SAM default, BAM possible) with the suffix *.deduplicated.bam.
        The output file appears in the same directory as the input file, but this method will move the output file
        to its own directory.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []
        for sample in self.samples:
            work_dir = os.path.join('dedup', sample.name)
            in_file = os.path.join('merged', sample.name, sample.name + '.merged.bam')
            out_file = os.path.join(work_dir, sample.name + '.merged.deduplicated.bam')
            report_file = os.path.join(work_dir, sample.name + '.merged.deduplication_report.txt')
            protocol = sample.readsets[0].library
            run_type = sample.readsets[0].run_type

            # Get I/O
            if run_type == 'PAIRED_END':
                in_report_file = os.path.join('merged', sample.name, sample.name + '.merged_aligned_PE_report.txt')
                copy_report = os.path.join(work_dir, sample.name + '.merged.deduplication_aligned_PE_report.txt')
            else:
                in_report_file = os.path.join('merged', sample.name, sample.name + '.merged_aligned_SE_report.txt')
                copy_report = os.path.join(work_dir, sample.name + '.merged.deduplication_aligned_SE_report.txt')

            mkdir_job = Job(command='mkdir -p ' + work_dir)
            copy_job = Job([in_report_file],
                           [copy_report],
                           command="cp -fu " + in_report_file + " " + copy_report)

            # A job name that is different from the heading will not use the params listed. (Uses default, instead)
            if protocol == 'RRBS':  # Deduplication is not recommended for RRBS datatypes. Keep what we have
                # You can only make a relative link in the current directory, so use absolute paths.
                abs_in_file = os.path.join(self.output_dir, in_file)
                abs_out_file = os.path.join(self.output_dir, out_file)
                job = concat_jobs([mkdir_job,
                                   Job([in_file], [out_file], command="cp -Ls -f " + abs_in_file + " " + abs_out_file),
                                   copy_job],
                                  name="skip_rrbs_deduplicate." + sample.name)
            else:  # WGBS
                merge_job = Job([in_file],
                                module_entries=[['bismark_deduplicate', 'module_samtools'],
                                                ['bismark_deduplicate', 'module_perl'],
                                                ['bismark_deduplicate', 'module_bismark']],
                                command="""deduplicate_bismark {type} --bam {other} {input}""".format(
                                    type='--paired' if run_type == 'PAIRED_END' else '--single',
                                    input=in_file,
                                    other=config.param('bismark_deduplicate', 'other_options', required=False)))
                move_bam = Job(output_files=[out_file],
                               command='mv -fu ' + os.path.join(os.path.dirname(in_file),
                                                                os.path.basename(out_file)) + ' ' + out_file,
                               removable_files=[out_file])
                move_log = Job(output_files=[report_file],
                               command='mv -fu ' + os.path.join(os.path.dirname(in_file),
                                                                os.path.basename(report_file)) + ' ' + report_file)
                job = concat_jobs([mkdir_job, merge_job, move_bam, move_log, copy_job],
                                  name='bismark_deduplicate.' + sample.name)
            jobs.append(job)
        return jobs

    def calc_dedup_nucleotide_coverage(self):  # Step 9
        """
        This step recalculates the nucleotide coverage values using Bismark's bam2nuc script. This recalculation is
        done after deduplication. The removal of some reads will likely alter this value, so this is helps update it.

        The output file is found in the related folder under dedup/.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            if sample.readsets[0].library == 'RRBS':
                continue  # No deduplication done for RRBS samples
            output_dir = os.path.join('dedup', sample.name)
            merged_bam = os.path.join(output_dir, sample.name) + '.merged.deduplicated.bam'
            job = bam2nuc_job(output_dir, sample.name, '.merged.deduplicated', merged_bam)
            jobs.append(job)
        return jobs

    def bismark_methylation_caller(self):  # Step 10
        """
        This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
        The following input files are accepted:
            1.	Bismark result files from previous alignment step
            2.	BAM files (unsorted) from readset file

        Input: Merged sample files (merged/)
        Output: Methylation calls in BedGraph format. (methyl_calls/)

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []

        for sample in self.samples:
            # Either select aligned sample from previous alignment step or aligned BAM/SAM files in readset file
            merged_sample = self.select_input_files([[readset.bam for readset in sample.readsets], [
                os.path.join("dedup", sample.name, sample.name + '.merged.deduplicated.bam')]])
            report_files = [
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bedGraph.gz"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.M-bias.txt"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated_splitting_report.txt"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.CpG_report.txt.gz")]
            other_files = [
                os.path.join("methyl_calls", sample.name, "CHG_OB_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHG_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHH_OB_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHH_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CpG_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CpG_OB_" + sample.name + ".merged.deduplicated.txt.gz")]
            run_type = sample.readsets[0].run_type

            job = Job(
                merged_sample + ['bismark_prepare_genome/Bisulfite_Genome'],
                report_files + other_files,
                [['bismark_methylation_caller', 'module_samtools'],
                 ['bismark_methylation_caller', 'module_perl'],
                 ['bismark_methylation_caller', 'module_bismark']],
                command="""\
mkdir -p {directory}
bismark_methylation_extractor {library_type} {other} --multicore {core} --output {directory} \
--bedGraph --cytosine_report --gzip --genome_folder {genome} {sample}""".format(
                    directory=os.path.join("methyl_calls", sample.name),
                    library_type="--paired-end" if run_type == "PAIRED_END" else "--single-end",
                    other=config.param("bismark_methylation_caller", "other_options"),
                    core=config.param('bismark_methylation_caller', 'cores'),
                    sample=" ".join(merged_sample),
                    genome=os.path.join(self.output_dir, 'bismark_prepare_genome')),
                removable_files=other_files,
                name="bismark_methylation_caller." + sample.name)
            jobs.append(job)
        return jobs

    def bismark_html_report_generator(self):  # Step 11
        """
        Generates the Bismark Report page by combining data from alignment, deduplication, methylation, and
        nucleotide coverage reports. The alignment report is a requirement while all others are optional.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Generation Setup
        report_file = 'report/EpiSeq.bismark_html_report_generator.md'
        report_data = 'data/bismark_html_report_generator'
        #template_string_file = 'bismark_summary_report/' + \
        #                       datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {report} | {methyl_calls} | {start} | {completion} |'

        jobs = []
        module_list = [['bismark_html_report_generator', 'module_samtools'],
                       ['bismark_html_report_generator', 'module_perl'],
                       ['bismark_html_report_generator', 'module_bismark']]
        for sample in self.samples:
            report_list = ['', '', '', '', '']
            # Required file
            if sample.readsets[0].run_type == 'PAIRED_END':
                report_list[0] = os.path.join("merged", sample.name,
                                              sample.name + ".merged_aligned_PE_report.txt")
            else:
                report_list[0] = os.path.join("merged", sample.name,
                                              sample.name + ".merged_aligned_SE_report.txt")
            # Optional output depending on run type
            if sample.readsets[0].library != 'RRBS':
                report_list[1] = os.path.join('dedup', sample.name, sample.name + '.merged.deduplication_report.txt')
                report_list[4] = os.path.join('dedup', sample.name, sample.name +
                                              '.merged.deduplicated.nucleotide_stats.txt')
            else:
                report_list[1] = None
                report_list[4] = os.path.join('merged', sample.name, sample.name +
                                              '.merged.nucleotide_stats.txt')
            # Based on our pipeline, these files are always generated, so make it manditory
            report_list[2] = os.path.join("methyl_calls", sample.name, sample.name +
                                          ".merged.deduplicated_splitting_report.txt")
            report_list[3] = os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.M-bias.txt")

            # Output file!
            html_report = os.path.join('bismark_summary_report', sample.name + '_final_bismark_report.html')

            # Job creation
            mkdir_job = Job(command='mkdir -p ' + os.path.dirname(html_report))
            job = Job(input_files=report_list, output_files=[html_report],
                      module_entries=module_list,
                      command="""bismark2report -o {out} --verbose --alignment_report {align} \
                      {dedup} {split} {mbias} {nt}""".format(
                          out=html_report,
                          align=report_list[0],
                          dedup=' --dedup_report ' + report_list[1] if report_list[1] else '',
                          split=' --splitting_report ' + report_list[2] if report_list[2] else '',
                          mbias=' --mbias_report ' + report_list[3] if report_list[3] else '',
                          nt=' --nucleotide_report ' + report_list[4] if report_list[4] else ''))

            # Report Creation
            zip_file = os.path.join(report_data, sample.name + '_methyl.zip')
            report_entry = fill_in_entry.format(
                sample=sample.name,
                report='[View HTML](' + os.path.join(report_data, os.path.basename(html_report)) + ')',
                methyl_calls='[Download Raw Data](' + zip_file + ')',
                start="$START",
                completion="$(date '+%Y-%m-%d %H:%M:%S')")
            command = """\
            TEMPLATE_STR_FILE=bismark_summary_report/$(date +%F)_template_var_strings.txt && \\
            flock -x "${{TEMPLATE_STR_FILE}}.lock" -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
            mkdir -p {data_dir} && \\
            cp -f {report} {data_dir} && \\
            zip {zip_file} {methyl_calls} && \\
            table=$(cat $TEMPLATE_STR_FILE) && \\
            pandoc \\
            {report_template_dir}/{basename_report_file} \\
                --template {report_template_dir}/{basename_report_file} \\
                --variable data_table="$table" \\
                --to markdown \\
            > {report_file}""".format(
                entry=report_entry,
                data_dir=os.path.join('report', report_data),
                methyl_calls='methyl_calls/' + sample.name + '/' + sample.name + '.merged.deduplicated.*',
                zip_file=os.path.join('report', zip_file),
                report=html_report,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            )
            # join zip_file path to report/
            report_job = Job(output_files=[report_file,
                                           os.path.join('report', report_data, os.path.basename(html_report)),
                                           os.path.join('report', report_data, os.path.basename(zip_file))],
                             report_files=[report_file],
                             module_entries=[['bismark_html_report_generator', 'module_pandoc']],
                             command=command)
            jobs.append(concat_jobs([Job(command="START=$(date '+%Y-%m-%d %H:%M:%S')"),
                                     mkdir_job, job, report_job], name='bismark_report.' + sample.name))
        return jobs

    def methylation_values(self):

        report_file = 'report/EpiSeq.methylation_values.md'
        beta_metrics_file = os.path.join("methylation_values", "methylation_metrics.tsv")
        beta_beanplot_file = os.path.join('report', 'beta_distribution.png')
        beta_file = os.path.join('methylation_values', 'methylation_values.csv')
        report_data = 'report/data/methylation_values'

        cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz")
                     for sample in self.samples]  # Input files
        sample_group = [sample.name for sample in self.samples]

        command="""\
mkdir -p {directory} && \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(BiSeq))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(minfi))
rrbs <- readBismark(c{samples}, colData=DataFrame(group=factor(c{group}), row.names=c{sample_names}))
rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any (x > {coverage})),]
beta <- methLevel(rawToRel(rrbs.filtered))

metrics <- data.frame(total.pos=apply(beta, 2, function(x) length(x[!is.na(x)]) ), row.names=colnames(beta))
metrics$num.na <- apply(beta, 2, function(x) length(x[is.na(x)]) )
cat(kable(metrics, row.names=TRUE), file="{beta_metrics_file}", sep="\n")

write.csv(beta, file="{beta_file}", quote=FALSE, row.names=FALSE)

png('{beta_beanplot_file}')
densityBeanPlot(as.matrix(beta))
dev.off()

EOF

mkdir -p {data_dir} && \\
cp -f {beta_file} {data_dir}; \\

table=$(cat {beta_metrics_file}) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable metrics_table="$table" \\
    --to markdown > {report_file}""".format(
            directory=os.path.dirname(beta_file),
            samples=tuple(cov_files),
            group=tuple(sample_group),
            sample_names=tuple([sample.name for sample in self.samples]),
            coverage=config.param("methylation_value_metrics", "read_coverage"),
            beta_file=beta_file,
            beta_metrics_file=beta_metrics_file,
            data_dir=report_data,
            report_template_dir=self.report_template_dir,
            basename_report_file=os.path.basename(report_file),
            beta_beanplot_file=beta_beanplot_file,
            report_file=report_file)

        job = Job(
            cov_files,
            [beta_file],
            [
                ["methylation_value_metrics", "module_R"],
                ["methylation_value_metrics", "module_pandoc"]
            ],
            command = command,
            report_files=[report_file],
            name="methylation_value_metrics")

        return [job]

    def differential_methylated_pos(self):
        """
    This step finds a list of differentially methylated CpG sites with respect to a categorical
    phenotype (controls vs. cases). The BedGraph files from the previous methylation calling step are first combined
    to a BSRaw object with the R package BiSeq. Then, the dmpFinder function from the R package minfi is used to
    compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
    for each site with the option of correcting them for multiple testing. Differential analysis is done for each
    contrast specified in the design file

    The values from the design files dictates how the samples are treated and compared.

    Input: Methylation data (methyl_calls/)
    Output: A CSV file in differential_methylated_positions/

    :return: A list of jobs that needs to be executed in this step.
    :rtype: list(Job)
    """

        # Report file variables
        report_file = 'report/EpiSeq.differential_methylated_pos.md'
        report_data = 'report/data/differential_methylated_pos/'
        beta_file = os.path.join('methylation_values', 'methylation_values.csv')
        fill_in_entry = '| {contrast_name} | {contrast_view} |'

        jobs = []
        for contrast in self.contrasts:
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]

            # Get file paths
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz")
                     for sample in contrast_samples]  # Input files
            dmps_file = os.path.join("differential_methylated_positions",
                contrast.name + "_RRBS_differential_methylated_pos.csv")  # Output file
            metrics_file = os.path.join("differential_methylated_positions",
                contrast.name + "_metrics.tsv")

            # Abort analysis if not enough samples (Will cause dmpFinder to throw an error)
            if len(contrast.controls) == 0 or contrast.treatments == 0 or len(contrast_samples) <= 2:  # No 1v1 or less
                log.warn("Insufficient sample size to compare case and control. Skipping contrast: " + contrast.name)
                continue

            report_entry = fill_in_entry.format(
                contrast_name = contrast.name,
                contrast_view = '[download csv]({link})'.format(
                    link = os.path.join(report_data, os.path.basename(dmps_file))
                )
            )

            command = """\
TEMPLATE_STR_FILE=differential_methylated_positions/$(date +%F)_template_var_strings.txt && \\
mkdir -p {directory} && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(BiSeq))
suppressPackageStartupMessages(library(minfi))

rrbs <- readBismark(c{samples}, colData=DataFrame(group=factor(c{group}), row.names=c{sample_names}))
rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
beta <- methLevel(rawToRel(rrbs.filtered))

#Use M values to do statistical tests because they are more reliable
#dmpFinder does not work with M values that are 0 or INF so the beta values must be shifted slightly
#Although there is no such thing as a beta value > 1, it will not matter in this step because only
#the average beta values are shown to the user
beta[beta == 0] = 0.000001
beta[beta == 1] = 0.999999
M <- log2(beta/(1-beta))

dmp <- dmpFinder(M, pheno=colData(rrbs.filtered)[,"group"], type="categorical")
dmp["pval.adjusted"] <- p.adjust(dmp[,"pval"], method = "{padjust_method}")
dmp <- dmp[dmp["pval.adjusted"] < {pvalue},][c("pval", "pval.adjusted")]

controls <- c({controls})
cases <- c({cases})
result = as.data.frame(rowRanges(rrbs.filtered))[1:4]
result["Avg Control Beta"] = rowMeans(beta[,controls])
result["Avg Case Beta"] = rowMeans(beta[,cases])
result["Avg Delta Beta"] = result[,"Avg Case Beta"] - result[,"Avg Control Beta"]
result <- cbind(result, beta)
result <- merge(result, dmp, by=0)

result <- result[abs(result["Avg Delta Beta"]) > {delta_beta_threshold},]

write.csv(result, file="{dmps_file}", quote=FALSE, row.names=FALSE)

EOF
mkdir -p {data_dir} && \\
cp -f {dmps_file} {data_dir}; \\
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown > {report_file}""".format(
                entry=report_entry,
                directory=os.path.dirname(dmps_file),
                samples=tuple(cov_files),
                group=tuple(sample_group),
                sample_names=tuple([sample.name for sample in contrast_samples]),
                coverage=config.param("differential_methylated_pos", "read_coverage"),
                beta_file=beta_file,
                controls=', '.join(["'" + sample.name + "'" for sample in contrast.controls]),
                cases=', '.join(["'" + sample.name + "'" for sample in contrast.treatments]),
                padjust_method=config.param("differential_methylated_pos", "padjust_method"),
                pvalue=config.param("differential_methylated_pos", "pvalue", type="float"),
                delta_beta_threshold=config.param("differential_methylated_pos", "delta_beta_threshold",
                    type="float"),
                dmps_file=dmps_file,
                data_dir=report_data,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file,
                contrast_name=contrast.name)

            job = Job(
                cov_files,
                [dmps_file],
                [
                    ["differential_methylated_pos", "module_R"],
                    ["differential_methylated_pos", "module_mugqic_R_packages"],
                    ["differential_methylated_pos", "module_pandoc"]
                ],
                command=command,
                report_files=[report_file],
                name="differential_methylated_pos." + contrast.name)

            jobs.append(job)

            # metrics and report
            cases = [sample.name for sample in contrast.treatments]
            controls = [sample.name for sample in contrast.controls]
            metrics_job = metrics.dmp_metrics(dmps_file, beta_file, cases, controls, 'report', contrast.name)
            metrics_report = 'report/EpiSeq.dmp_metrics.md'
            metrics_job.report_files = [metrics_report]

            metrics_table_file = os.path.join('report', 'dmp.metrics.table')
            graphs = metrics_job.output_files
            img_list = "\n\n".join(['![](' + os.path.basename(graph) + ')' for graph in graphs])

            pandoc_command = """
TEMPLATE_STR_FILE=differential_methylated_positions/$(date +%F)_metrics_template_strings.txt && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{img_list}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
metrics_table=$(cat {metrics_table_file}) && \\
img_list=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable metrics_table="$metrics_table" \\
    --variable img_list="$img_list" \\
    --variable pval_cutoff={pvalue} \\
    --variable delta_cutoff={delta_cutoff} \\
    --to markdown > {report_file}""".format(
            report_template_dir=self.report_template_dir,
            basename_report_file=os.path.basename(metrics_report),
            metrics_table_file=metrics_table_file,
            img_list=img_list if len(self.contrasts) == 1 else "####" + contrast.name + img_list + "\n",
            report_file=metrics_report,
            pvalue=config.param("differential_methylated_pos", "pvalue", type="float"),
            delta_cutoff=config.param("differential_methylated_pos", "delta_beta_threshold",
                    type="float")
        )
            pandoc_job = Job(command=pandoc_command, module_entries=[['differential_methylated_pos', 'module_pandoc']])
            jobs.append(concat_jobs([metrics_job, pandoc_job], name="dmp_metrics." + contrast.name))

        return jobs


    def differential_methylated_regions(self):
        """
    Similar to differential_methylated_positions, this step looks at methylation patterns on a larger, regional
    level. This step compares large-scale differences in methylation as opposed to comparing local methylation
    sites.

    Input: Methylation data (methyl_calls/)
    Output: A CSV file in differential_methylated_regions/

    :return: A list of jobs that needs to be executed in this step.
    :rtype: list(Job)
    """
        # Report file variables
        report_file = 'report/EpiSeq.differential_methylated_regions.md'
        report_data = 'report/data/differential_methylated_regions/'
        fill_in_entry = '| {contrast_name} | {contrast_view} |'
        beta_file = os.path.join("methylation_values", "methylation_values.csv")

        jobs = []
        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]

            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz")
                     for sample in contrast_samples]  # Input files
            dmrs_file = os.path.join("differential_methylated_regions",
                contrast.name + "_RRBS_differential_methylated_regions.csv")

            report_entry = fill_in_entry.format(
                contrast_name = contrast.name,
                contrast_view = '[download csv]({link})'.format(
                    link = os.path.join(report_data, os.path.basename(dmrs_file))
                )
            )

            command = """\
TEMPLATE_STR_FILE=differential_methylated_regions/$(date +%F)_template_var_strings.txt && \\
mkdir -p {directory} && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}"; \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(bumphunter))
suppressPackageStartupMessages(library(BiSeq))
library(doParallel)
registerDoParallel(cores={cores})

rrbs <- readBismark(c{samples}, colData=DataFrame(group=c{group}, row.names=c{sample_names}))
rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
beta <- methLevel(rawToRel(rrbs.filtered))
chr <- as.character(seqnames(rowRanges(rrbs.filtered)))
pos <- start(ranges(rowRanges(rrbs.filtered)))
pheno <- colData(rrbs.filtered)[,"group"]
designM <- model.matrix(~pheno)

dmrs <- bumphunterEngine(beta,
                         chr=chr,
                         pos=pos,
                         design=designM,
                         cutoff={delta_beta_threshold},
                         pickCutoffQ=0.99,
                         null_method=c("permutation","bootstrap"),
                         smooth=FALSE,
                         smoothFunction=locfitByCluster,
                         B={permutations},
                         verbose=TRUE,
                         maxGap=500)

dmrs <- na.omit(dmrs)

write.csv(dmrs$tab, "{dmrs_file}", quote=FALSE, row.names=FALSE)

EOF
mkdir -p {data_dir} && \\
cp -f {dmrs_file} {data_dir}; \\
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown > {report_file}""".format(
                entry=report_entry,
                directory=os.path.dirname(dmrs_file),
                samples=tuple(cov_files),
                group=tuple(sample_group),
                coverage=config.param("differential_methylated_regions", "read_coverage", type="int"),
                sample_names=tuple([sample.name for sample in contrast_samples]),
                cores=config.param('bismark_methylation_caller', 'cluster_cpu').split('=')[-1],
                delta_beta_threshold=config.param("differential_methylated_regions", "delta_beta_threshold",
                    type="float"),
                permutations=config.param("differential_methylated_regions", "permutations", type="int"),
                dmrs_file=dmrs_file,
                data_dir=report_data,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file)

            job = Job(
                cov_files,
                [dmrs_file],
                [
                    ["differential_methylated_regions", "module_R"],
                    ["differential_methylated_regions", "module_mugqic_R_packages"],
                    ["differential_methylated_regions", "module_pandoc"]
                ],
                command=command,
                report_files=[report_file],
                name="differential_methylated_regions." + contrast.name)


            jobs.append(job)

            metrics_job = metrics.dmr_metrics(dmrs_file, beta_file, 'report', contrast.name)
            metrics_report = 'report/EpiSeq.dmr_metrics.md'
            metrics_job.report_files = [metrics_report]

            graphs = metrics_job.output_files
            img_list = "\n\n".join(['![](' + os.path.basename(graph) + ')' for graph in graphs])

            pandoc_command = """\
TEMPLATE_STR_FILE=differential_methylated_regions/$(date +%F)_metrics_template_stings.txt && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{img_list}\\" >> ${{TEMPLATE_STR_FILE}}" && \\
img_list=$(cat $TEMPLATE_STR_FILE) &&
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable img_list="$img_list" \\
    --to markdown > {report_file}""".format(
                img_list=img_list if len(self.contrasts)==1 else "####" + contrast.name + img_list + "\n",
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(metrics_report),
                report_file=metrics_report
            )
            pandoc_job = Job(command=pandoc_command, module_entries=[['differential_methylated_pos', 'module_pandoc']])

            jobs.append(concat_jobs([metrics_job, pandoc_job], name="dmr_metrics." + contrast.name))


        return jobs

    def prepare_annotations(self):
        annotations_file = "prepare_annotations/annotations.Rdata"
        gtf_file = config.param('prepare_annotations', 'gtf')

        command = """\
mkdir -p {directory} && \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(bumphunter))

txdb <- makeTxDbFromGFF('{gtf_file}', format="gtf", organism='{organism}')

# If not present, add 'chr' prefix to sequence names in txdb to match dmr/dmp format
regexp <- '^([XYM]|MT|[1-9][0-9]*)$'
s <- seqlevels(txdb)
s <- unlist(lapply(s, function(x) if (grepl(regexp, x)) paste('chr', x, sep='') else x))
txdb <- renameSeqlevels(txdb, s)

annotations <- annotateTranscripts(txdb, by='tx', codingOnly=FALSE)

save(annotations, file='{annotations_file}')
EOF""".format(
            directory=os.path.dirname(annotations_file),
            gtf_file=gtf_file,
            organism=config.param('prepare_annotations', 'scientific_name').replace('_', ' '),
            annotations_file=annotations_file
        )

        return [Job(
            [gtf_file],
            [annotations_file],
            [['prepare_annotations', 'module_R']],
            command=command,
            name="perpare_annotations"
        )]

    def annotate_regions(self):
        annotations_file = "prepare_annotations/annotations.Rdata"
        report_file = 'report/EpiSeq.annotate_regions.md'
        report_data = 'report/data/annotate_regions'
        fill_in_entry = '| {contrast_name} | [download csv]({contrast_data}) |'

        jobs = []
        for contrast in self.contrasts:
            dmr_file = os.path.join("differential_methylated_regions",
                                    contrast.name + "_RRBS_differential_methylated_regions.csv")
            matched_file = os.path.join("annotate_regions", contrast.name + ".matched_genes.csv")

            report_entry = fill_in_entry.format(
                contrast_name = contrast.name,
                contrast_data = os.path.join(report_data, os.path.basename(matched_file))
            )

            command = """\
TEMPLATE_STR_FILE=annotate_regions/$(date +%F)_template_var_strings.txt && \\
mkdir -p {directory} && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}"; \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(bumphunter))
library(doParallel)
registerDoParallel(cores={cores})

dmrs <- read.csv('{dmr_file}')

load('{annotations_file}')

# remove any regions that are on a sequence not in the annotations
# required as matchGenes crashes if any regions fail to map to a gene
dmrs <- dmrs[dmrs$chr %in% seqlevelsInUse(annotations), ]

matched.genes <- matchGenes(dmrs, annotations,
    promoterDist={promoterDist}, type='{type}', skipExons={skipExons})

annotated <- cbind(dmrs, matched.genes)
write.csv(annotated, file='{matched_file}')
EOF
mkdir -p {data_dir} && \\
cp -f {matched_file} {data_dir} && \\
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown > {report_file}""".format(
                directory=os.path.dirname(matched_file),
                cores=config.param('annotate_regions', 'cluster_cpu').split('=')[-1],
                dmr_file=dmr_file,
                annotations_file=annotations_file,
                promoterDist=config.param('annotate_regions', 'promoter_distance'),
                type=config.param('annotate_regions', 'distance_type'),
                skipExons=str(config.param('annotate_regions', 'skip_exons', type='boolean')).upper(),
                matched_file=matched_file,
                entry=report_entry,
                data_dir=report_data,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file)

            job = Job(
                [dmr_file, annotations_file],
                [matched_file],
                [
                    ['annotate_regions', 'module_R'],
                    ['annotate_regions', 'module_pandoc']
                ],
                command=command,
                name="annotate_regions." + contrast.name
            )

            jobs.append(job)

        return jobs

    def annotate_positions(self):
        annotations_file = "prepare_annotations/annotations.Rdata"
        report_file = 'report/EpiSeq.annotate_positions.md'
        report_data = 'report/data/annotate_positions'
        fill_in_entry = '| {contrast_name} | [download csv]({contrast_data}) |'

        jobs = []
        for contrast in self.contrasts:
            dmp_file = os.path.join("differential_methylated_positions",
                                    contrast.name + "_RRBS_differential_methylated_pos.csv")
            matched_file = os.path.join("annotate_positions", contrast.name + ".matched_genes.csv")

            report_entry = fill_in_entry.format(
                contrast_name = contrast.name,
                contrast_data = os.path.join(report_data, os.path.basename(matched_file))
            )

            command = """\
TEMPLATE_STR_FILE=annotate_positions/$(date +%F)_template_var_strings.txt && \\
mkdir -p {directory} && \\
flock -x ${{TEMPLATE_STR_FILE}}.lock -c "echo \\"{entry}\\" >> ${{TEMPLATE_STR_FILE}}"; \\
mkdir -p {directory} && \\
R --no-save --no-restore <<-'EOF'
suppressPackageStartupMessages(library(bumphunter))
library(doParallel)
registerDoParallel(cores={cores})

dmps <- read.csv('{dmp_file}')

load('{annotations_file}')

# remove any positions that are on a sequence not in the annotations
# required as matchGenes crashes if any regions fail to map to a gene
dmps <- dmps[dmps$seqnames %in% seqlevelsInUse(annotations), ]

matched.genes <- matchGenes(dmps, annotations,
    promoterDist={promoterDist}, type='{type}', skipExons={skipExons})

annotated <- cbind(dmps, matched.genes)
write.csv(annotated, file='{matched_file}')
EOF
mkdir -p {data_dir} && \\
cp -f {matched_file} {data_dir} && \\
table=$(cat $TEMPLATE_STR_FILE) && \\
pandoc \\
    {report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown > {report_file}""".format(
                directory=os.path.dirname(matched_file),
                cores=config.param('annotate_positions', 'cluster_cpu').split('=')[-1],
                dmp_file=dmp_file,
                annotations_file=annotations_file,
                promoterDist=config.param('annotate_positions', 'promoter_distance'),
                type=config.param('annotate_positions', 'distance_type'),
                skipExons=str(config.param('annotate_positions', 'skip_exons', type='boolean')).upper(),
                matched_file=matched_file,
                entry=report_entry,
                data_dir=report_data,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file)

            job = Job(
                [dmp_file, annotations_file],
                [matched_file],
                [
                    ['annotate_positions', 'module_R'],
                    ['annotate_positions', 'module_pandoc']
                ],
                command=command,
                name="annotate_positions." + contrast.name
            )

            jobs.append(job)

        return jobs


if __name__ == '__main__':
    EpiSeq()
