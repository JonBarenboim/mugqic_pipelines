Episeq Pipeline
===============

The Episeq pipeline takes methylation sequencing data and conducts a differential methylation analysis among the given samples. The pipeline takes FASTQ or BAM files (unsorted) as input and produces an differential analysis in the methylome. Currently, only WGBS and RRSB data are supported.

You should have filled out the readset, design and configuration (`*.ini`) file. Please refer to the Episeq User Manual for more information and instructions.

[TOC]

Quick Start
-----------
1. Make a symbolic symlink to the episeq pipeline folder from your working directory. 
    
    `cp -sL <path-to-mugqic-pipelines>/pipelines/episeq .`
    
1. Either create a directory **or** symbolic symlink named `data` that contains your sample data to your working directory.

1. Copy sample configuration files from the episeq folder to your working directory

    `cp episeq/example/episeq.ini episeq/example/episeq.design episeq/example/episeq.readset episeq/example/generate-qsub.sh .`

1. Using your prefered text editor, modify `./episeq.readset` and `./episeq.design`  to add your samples and define your case/control groups, respectively.
    - If you have a manifest file with the following columns, you may instead use the `episeq/readset_gen.py` utility to generate personalised readset and design files.
        - SRA_Sample_s
        - Run_s
        - LibraryLayout_s
        - SRA_Study_s
        - Assay_Type_s
        - LibrarySelection_s
    
        `./episeq/readset_gen.py <manifest_file> <output_readset_file> <output_design_file> <root data directory>`
        
        **Note:** If using script, the data files should be organized as `<data_dir>/<sample_name>/<sample_name>*.{bam,fastq}`

1. Run `./generate-qsub.sh` and check the `./debug.log` and `./qsub.sh` files for any errors.

1. When everything is set, run `./qsub.sh` to start the pipeline.

1. Once the entire pipeline has completed successfully, add `--report` in `generate-qsub.sh` and repeat the previous two steps to generate a report for the entire pipeline. It will be made in the reports directory, under the output directory.

Usage
-----
```
usage: episeq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [-d DESIGN]
                 [-r READSETS] [-v]

Version: 2.2.0

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch}, --job-scheduler {pbs,batch}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
 1- bismark_prepare_genome,
 2- pre_qc_check,
 3- trim_galore,
 4- bismark_align,
 5- merge_bismark_alignment_report,
 6- picard_merge_sam_files,
 7- merged_nuc_stats,
 8- bismark_deduplicate,
 9- calc_dedup_nucleotide_coverage,
10- bismark_methylation_caller,
11- bismark_html_report_generator,
12- methylation_values
13- differential_methylated_pos,
14- differential_methylated_regions
15- prepare_annotations
16- annotate_positions
17- annotate_regions
18- position_enrichment_analysis
19- region_enrichment_analysis
```

## Reports
The pipeline generates a massive amount of data, many of which is not immediately visualized. Reports help provide a summary of results that is more readily understood. As such, it is one of the easier ways towards understanding the results. Thus, it is our responsibility to ensure that the data is organized in a way that others can look over the results. Essential information for each step can be found [here](#epi-seq-pipeline-steps).


### Processing and Quality Control
The [Bismark HTML report](#11-bismark_html_report_generator) provides an abundance of information the pre-processing steps that each sample went through. This visual report helps summarize the output for steps 4-10 of the pipeline. From step 11 alone, this report gives the user a good understanding on any possible problems with the input data.

## Epi-Seq Pipeline Steps
**NOTE:** For each step below, the logs containing stdout and stderr messages can be found at `job_output/<output_directory_name>`.
### 1. bismark_prepare_genome
This step takes in a reference genome fasta file and convert the sequence for methylation alignment and sequencing. This is a pre-processing step for [bismark_align](#4-bismark_align). The step will link the reference genome to the output directory (if needed), and create the methylome sequence in the directory called `Bisulfite_Genome`, which contains two subdirectories within it. Each subdirectory contains a converted sequence of either C->T and G->A conversions. Each of these converted sequences are used to align the bisulfite-treated reads. Since bisulfite sequencing causes unmethylated C to covert to U and later interpreted as T, this step allows alignment to be made without excessive mismatches due to the bisulfite treatment. This step only needs to be done once within a project's output folder.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `bismark_prepare_genome` |
| Job name prefix | `bismark_prepare_genome` |
| Requires | None
| Blocks | [bismark_align](#4-bismark_align) <br /> [bismark_methylation_caller](#10-bismark_methylation_caller) |


__Note:__ Depending on the size of the genome, this step can take several hours to complete.

### 2. pre_qc_check
This step gives the user information about the input data before any processing steps are done. This helps everyone determine if a problem comes from the pipeline or the sequencing experiment. This output is recorded in the final report.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `pre_qc_check` |
| Job name prefix | `pre_qc_check` |
| Requires | None |
| Blocks | None |

### 3. trim_galore
This step helps improve the alignment efficiency by performing quality trimming via the open source package Trim Galore!. Briefly, this package addresses many of the issues that occur when analyzing bisulfite sequencing libraries. The pipeline does trimming using only the default options in Trim Galore. Additional options such as stricter or more relaxed trimming can be entered through the `other_options` parameter in the configuration file. Output files are gzipped `.fq` files. This step can be ignored if the readset has a `.BAM` file. This step is recorded in the final pipeline report.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `trimmed` |
| Job name prefix | `trim_galore` |
| Requires | None |
| Blocks | [bismark_align](#4-bismark_align) |

### 4. bismark_align
This step aligns the trimmed reads to a reference genome from `bismark_prepare_genome`. The alignment is done by the open source package Bismark using the default stringency settings. By default, the settings can be somewhat strict, but is essential to avoid mismatches from sequencing error. Additional options can be entered through the "`other_options`" parameter in the configuration file. Output files are `.bam` files, but may be configured to output `cram` or `sam` files. This step can be ignored if the readset contains `.BAM` files. This step is recorded in the final pipeline report.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `aligned` |
| Job name prefix | `bismark_align` |
| Requires | [bismark_prepare_genome](#1-bismark_prepare_genome) <br /> [trim_galore](#3-trim_galore) |
| Blocks | [merge_bismark_alignment_report](#5-merge_bismark_alignment_report) <br /> [picard_merge_sam_files](#6-picard_merge_sam_files) |


### 5. merge_bismark_alignment_report
So far, the pipeline has been handling data at a readset level. However, our analysis requires information on a sample level basis. Thus, we begin to collate all of the information we have collected so for. This step focuses on merging the alignment report that yields information about the mapping quality of the reads.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `merge_bismark_alignment_report` |
| Requires | [bismark_align](#4-bismark_align) |
| Blocks | [bismark_html_report_generator](#11-bismark_html_report_generator) |

### 6. picard_merge_sam_files
This step combines all readsets together to form a single `.bam` file for each sample. This is often required when multiple libraries or multiplexing is done before sequencing. Merging at this step allows us to parallelize the pipeline as much as we can before aggregating. `bismark_align` will merge any paired reads into a single file, which takes care some of the work for us.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `picard_merge_sam_files`<br/> `symlink_readset_sample_bam` |
| Requires | [bismark_align](#4-bismark_align) |
| Blocks | [bismark_deduplicate](#8-bismark_deduplicate) <br/> [merged_nuc_stats](#7-merged_nuc_stats) |

### 7. merged_nuc_stats
There isn't a dedicated merge script for the nucleotide coverage information, so we have to rerun the analysis here to generate a new report file. This step is (again) for QC and troubleshooting purposes. Looking at the nucleotide coverage after all readsets are merged may reveal biases that are found in some readsets, but not others.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `bam2nuc` |
| Requires | [picard_merge_sam_files](#6-picard_merge_sam_files) |
| Blocks | None |

### 8. bismark_deduplicate
Merging readsets reads can cause all sorts of artifacts and errors. Depending on the experiment, you may have duplicate reads that are frequently found across multiple readsets. Thus, deduplication will help eliminate some background noise that may occur. This yields a processed bam file that is ready for calling and downstream analysis.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `dedup` |
| Job name prefix | `bismark_deduplicate`<br/> `skip_rrbs_deduplicate` |
| Requires | [picard_merge_sam_files](#6-picard_merge_sam_files) |
| Blocks | [calc_dedup_nucleotide_coverage](#9-calc_dedup_nucleotide_coverage) <br/> [bismark_methylation_caller](#10-bismark_methylation_caller) <br/> [bismark_html_report_generator](#11-bismark_html_report_generator) |

### 9. calc_dedup_nucleotide_coverage
One more calculation is needed because the bam file has been modified since the last time this calculation has been done. The reason to keep each report file is for diagnostic and troubleshooting purposes. These calculations yield a metric that can be compared as the data moves through the pipeline. For this pipeline, this is the final read coverage information.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `dedup` |
| Job name prefix | `bam2nuc` |
| Requires | [bismark_deduplicate](#8-bismark_deduplicate) |
| Blocks | [bismark_html_report_generator](#11-bismark_html_report_generator) |

### 10. bismark_methylation_caller
This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
The following input files are accepted:

1.	Bismark result files from previous alignment step
1.	`BAM` files (unsorted) from readset file

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `methyl_calls` |
| Job name prefix | `bismark_methylation_caller` |
| Requires | [bismark_prepare_genome](#1-bismark_prepare_genome) <br/> [bismark_deduplicate](#8-bismark_deduplicate) |
| Blocks | [bismark_html_report_generator](#11-bismark_html_report_generator) <br/> [differential_methylated_pos](#12-differential_methylated_pos) <br/>[differential_methylated_regions](#13-differential_methylated_regions) |

### 11. bismark_html_report_generator
This job summarizes all data from steps 5-10 into one HTML report file. It contains diagrams that summarizes the various report files Bismark creates in it's proccessing toolkit. This can serve as an excellent overview for the quality of the sample data and could make all other Bismark output reports redundant. This step is recorded in the final pipeline report.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `bismark_summary_report` |
| Job name prefix | `bismark_report` |
| Requires | [merge_bismark_alignment_report](#5-merge_bismark_alignment_report) <br/> [bismark_deduplicate](#8-bismark_deduplicate) <br/> [calc_dedup_nucleotide_coverage](#9-calc_dedup_nucleotide_coverage) <br /> [bismark_methylation_caller](#10-bismark_methylation_caller)|
| Blocks | None |

### 12. methylation_values
This step reads reads methylation values for each sample. The `BedGraph` files from the previous methylation calling step are combined to a `BSRaw` object with the R package `BiSeq`. The `BSRaw` object is then converted to a `BSRel` object and saved as an R data file.

| Job Attribute | Value |
|:--------------|:----- |
| Output directory name: | `methylation_values` |
| Job name prefix | `methylation_values` |
| Requires | `bismark_methylation_caller` |
| Blocks | [differential_methylated_pos](#13-differential_methylated_pos) <br/> [differential_methylated_regions](#14-differential_methylated_regions) <br/> [position_enrichment_analysis](#18-position_enrichment_analysis) <br/> [region_enrichment_analysis](#19-region_enrichment_analysis) |

### 13. differential_methylated_pos
This step finds a list of differentially methylated CpG sites with respect to a categorical
phenotype (controls vs. cases). The `BSRel` object from the previous `methylation_values` step is loaded. Then, the `dmpFinder` function from the R package `minfi` is used to compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned for each site with the option of correcting them for multiple testing. Differential analysis is done for each contrast specified in the design file. This step is recorded in the final pipeline report. 

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `differential_methylated_pos` |
| Job name prefix | `differential_methylated_pos` <br/> `dmp_metrics` |
| Requires | [methylation_values](#12-methylation_values) |
| Blocks | [annotate_positions](#16-annotate_positions) <br/> [position_enrichment_analysis](#17-position_enrichment_analysis) |

### 14. differential_methylated_regions
This step finds a list of differentially methylated regions with respect to a categorical phenotype (controls vs. cases). The `BSRel` object from the previous `methylation_values` step is loaded, and then the bumphunting algorithm is run to locate regions of differential methylation. 

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `differential_methylated_regions` |
| Job name prefix | `differential_methylated_regions` <br/> `dmr_metrics` |
| Requires | [methylation_values](#12-methylation_values) |
| Blocks | [annotate_regions](#17-annotate_regions) <br/> [region_enrichment_analysis](#17-region_enrichment_analysis) |

__WARNING:__ This step is slow and requires large amounts of memory!

### 15. prepare_annotations
This step prepares gene and transcript annotations. A `GTF` file is used to prepare the annotations, but the pipeline can be modified to use annotations data from the Bioconductor repository. The transcripts are saved to an R data file in the form of a `GRanges` object

| Job Attribute | Value |
|:--------------|:------|
| Output directory name: | `prepare_annotations` |
| Job name prefix | `prepare_annotations` |
| Requires | None |
| Blocks | [annotate_positions](#16-annotate_positions) <br/> [annotate_regions](#17-annotate_regions) |

### 16. annotate_positions
This step annotates the CpG sites found in the previous `differential_methylated_pos` step using the annotations prepared in the previous `prepare_annotations` step. The CpG sites are matched to the annotated regions using the `matchGenes` function of the R package `bumphunter`

| Job Attribute | Value |
|:------------- |:----- |
| Output directory name: | `annotate_positions` |
| Job name prexix | `annotate_positions` | 
| Requires | [differential_methylated_pos](#13-differential_methylated_pos) <br/> [prepare_annotations](#15-prepare_annotations) |
| Blocks | None |

### 17. annotations_regions
This step annotates the regions found in the previous `differential_methylated_regions` step using the same methods as the `annotate_positions` step

| Job Attribute | Value |
|:------------- |:----- |
| Output directory name: | `annotate_regions` |
| Job name prefix | `annotate_regions` |
| Requires | [differential_methylated_regions](#14-differential_methylated_regions) <br/> [prepare_annotations](#15-prepare_annotations) |
| Blocks | None |

### 18. position_enrichment_analysis
This step tests overlap of CpG sites identified in the `differential_methylated_pos` step against regions sets. The R ackage `LOLA` is used to test for enrichment and order the region sets by statistical significance. The region sets are taken from two places. First, `bed` files may be supplied by the user. Second, the R script `LOLAsearch.r` is used to select region sets from the curated LOLAcore database by searching for keywords in metadata. The set of all CpG sites with a methylation value calculated in the `methylation_values` step is used as a background universe for the enrichment analysis.

| Job Attribute | Value | 
|:------------- |:----- |
| Output directory name | `enrichment_analysis` |
| Job name prefix | `position_enrichment_analysis` |
| Requires | [methylation_values](#12-methylation_values) <br/> [differential_methylated_pos](#13-differential_methylated_pos) |
| Blocks | None |

### 19 region_enrichment_analysis
This step performs enrichment analysis of regions identified in the `differential_methylated_regions` step using the methods as the `position_enrichment_analysis` step

| Job Attribute | Value |
|:------------- |:----- |
| Output directory name | `enrichment_analysis` |
| Job name prefix | `region_enrichment_analysis` |
| Requires | [methylation_values](#12-methylation_values) <br/> [differential_methylated_regions](#14-differential_methylated_regions) |
| Blocks | None |





