#!/usr/bin/env python

from os.path import join

# MUGQIC Modules
from core.config import config
from core.job import Job


def merge_overlapping_reads(fastq1, fastq2, output_dir, output_prefix):
    return Job(
        input_files=[fastq1, fastq2],
        output_files=[join(output_dir, output_prefix + '.notCombined_1.fastq'),
                      join(output_dir, output_prefix + '.notCombined_2.fastq'),
                      join(output_dir, output_prefix + '.extendedFrags.fastq')],
        module_entries=['flash', 'module_flash'],
        command='flash -d {output_dir} -o {output_prefix} {options} ' \
                '{fastq1} {fastq2}'.format(output_dir=output_dir,
                                           output_prefix=output_prefix,
                                           options=config.param('flash', 'options'),
                                           fastq1=fastq1,
                                           fastq2=fastq2))
