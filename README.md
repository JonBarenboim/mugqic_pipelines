# MUGQIC PIPELINES
-------------

This repository holds perl libs, wrappers and scripts of several bioinformatics pipelines.

The main repo dir organization is:

mugqic_pipeline  -  root dir

root/lib       - dir containing all libraries of all pipelines. 

root/pipelines - dir containing the pipelines itself. An addtional dir should be created in this dir with the pipeline name \(root/pipeline/Pipeline_name\)

Perl documentation on *.pm and *.pl files should \(as much as possible\) be created using POD \(Pod::Usage\). 


# Documentation 
-------------

Visit our [wiki](https://biowiki.atlassian.net/wiki/display/PS/Pipeline+Space+Home) for an overview of the pipelines included in the mugqic_pipeline repository. To automatically generate documentation for perl wrappers and libraries, use [pod2html](http://perldoc.perl.org/Pod/Html.html). For instance, 

`pod2html --infile=src/mugqic_pipeline/pipelines/chipseq/chipSeq.pl --outfile=chipSeq.html`

will generate HTML documentation for the chipSeq.pl script (The chipseq pipeline wrapper)

MUGQIC pipelines are perl programs that write to the standard output a list of bash commands intended to be run in a cluster computing environment. Dependencies between steps are controlled from the program. When executed, commands will submit jobs to a batch server. The cluster submit commands and the parameters associated to each job can be modified through a configuration file.


## Download and setup

In order to make pipelines that work anywhere, regardless of user path settings, and to control for versions of third party software and data associated to the analysis, the pipelines depend on environment modules. The modules created by the bioinformatics team in Compute Canada&#39;s Guillimin and Mammouth clusters are located in the &#36;MUGCIC_INSTALL_HOME/modulefiles path and use the notation: mugqic/&#60software&#62/&#60version&#62 . Click [here](https://biowiki.atlassian.net/wiki/display/CS/Software+and+Data+Dependencies) to setup modules on Guillimin and Mammouth clusters.


To download the pipeline use git to obtain the most recent development version. Mugqic pipelines are hosted on github, and can be obtained via:

`git clone https://bitbucket.org/mugqic/mugqic_pipeline.git`


If you encounter issues with the development version (master), or simply wish to obtain the most recent stable revision (1.0) then use:

` git fetch && git checkout 1.0`

Whenever the latest snapshot from github is needed, use the command pull

` git pull origin`

### Setup

Add the mugqic directory library to your PERL5LIB path.

`export PERL5LIB=&#36;{PERL5LIB}:/home/user/src/mugqic_pipeline/lib/`

## Usage

In its general operation all the mugqic pipelines require two input files: a project&#39;s read set sheet and a configuration (ini) file. Fastq files must be properly setup using a specific naming and directory structure. Additionally, the RNAseq and CHIPseq pipelines require a design file and the path for the job output logs. The command line options of the main pipelines are described below. Summaries of usage are printed when a command is run with no arguments. 

### The project&#39;s read set sheet 

Is the project read set sheet, generated from [NANUQ](http://gqinnovationcenter.com/index.aspx). See [this](https://biowiki.atlassian.net/wiki/display/PS/Read+Set+Files+%28FastQ%29+Setup) page to learn how to properly setup your fastq files and your project read set sheet.


### The configuration (ini) file. 
Is the standard configuration file for the pipeline. It&#39;s a plain text file composed of sections, keys and values. Section name appears on a line in square brackets ([default]). Every key has a name and a value, separated by an equals sign (=), i.e. clusterSubmitCmd=msub. Semicolons (;) or number signs (#) at the beginning of the line indicate a comment. Comment lines are ignored. Generally sections are associated to specific steps of the pipeline. If a property is associated to a specific step, the program will search for keys in the respective section. If the key is not found, the values in the default section will be used. Templates for Compute Canada&#39;s Guillimin and Mammouth clusters may already be available in the pipeline repository.


### DNAseq pipeline

`perl dnaSeq.pl -c dnaSeq.abacus.ini -s 1 -e 21 -n project.nanuq.csv > toRun.sh`

will generate a bash script for steps 1 to 21. This script can then be executed:

`sh toRun.sh`

#### Options

    -c dna.abacus.ini the standard configuration file for the pipeline. Templates for some cluster systems like Abacus or Guillimin may already be available  [here](https://bitbucket.org/mugqic/mugqic_pipeline/src/31db440f92036e3bae7e7a9676fdb57c3a9f3112/pipelines/dnaseq/?at=master)
    -s The start step
    -e The end step
    -n project.nanuq.csv the read set sheet, prepared as described above.

### RNAseq pipeline

`perl rnaSeq.pl -c rnaSeq.abacus.ini -s 1 -e 14 -n project.nanuq.csv -d design.txt -w  `pwd` > toRun.sh`

will generate a bash script for steps 1 to 14. This script can then be executed:

`sh toRun.sh `

#### Options

    -c rnaSeq.abacus.ini the standard configuration file for the pipeline. Templates for some cluster systems like Abacus or Guillimin may already be available  [here](https://bitbucket.org/mugqic/mugqic_pipeline/src/31db440f92036e3bae7e7a9676fdb57c3a9f3112/pipelines/rnaseq/?at=master)
    -s The start step
    -e The end step
    -n project.nanuq.csv  the NANUQ project read set sheet, prepared as described above.
    -d design.txt the design file. A tab separated value file that specifies the experimental design information of the project. The first column lists the sample names, which should match elements the column Name in the read set sheet. Subsequent columns specify all the pairwise comparisons which should be undertaken: values should be either "2" (nominator), "1" (denominator) or "0" (exclude from comparison). See example here. 
    -w The project&#39;s working directory. All job outputs will be sent to this directory.

### CHIPseq pipeline

`perl bin/chipSeq.pl -c chipSeq.abacus.ini -n project.nanuq.csv -d design.csv -w  ``pwd`` -s 1 -e 11 > toRun.sh

will generate a bash script for steps 1 to 11. This script can then be executed:

`sh toRun.sh`

#### Options

    -c chipSeq.abacus.ini the standard configuration file for the pipeline. Templates for some cluster systems like Abacus or Guillimin may already be available  [here](https://bitbucket.org/mugqic/mugqic_pipeline/src/31db440f92036e3bae7e7a9676fdb57c3a9f3112/pipelines/chipseq/?at=master)
    -s The start step
    -e The end step
    -n project.nanuq.csv the NANUQ Project sample file, generated automatically from the Read Set Search tab.
    -d design.csv the design file. A tab separated value file that specifies the experimental design information of the project.
    -w The current working directory. All job outputs will be sent to this directory.
