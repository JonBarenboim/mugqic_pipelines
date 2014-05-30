#!/usr/bin/env perl

=head1 NAME

I<rnaSeqDeNovoAssembly>

=head1 SYNOPSIS

Usage: perl rnaSeqDeNovoAssembly.pl -h | -c CONFIG_FILE -r step_range [-s SAMPLE_FILE] [-d DESIGN_FILE] [-o OUTPUT_DIR]
  -h  help and usage
  -c  .ini config file
  -r  step range e.g. "1-5", "3,6,7", "2,4-8"
  -s  sample file
  -d  design file
  -o  output directory (default: current)

=head1 DESCRIPTION

B<rnaSeqDeNovoAssembly.pl> is the main RNA-Seq De Novo assembly pipeline.

=head1 AUTHORS

B<David Morais> - I<dmorais@cs.bris.ac.uk>

B<Mathieu Bourgey> - I<mbourgey@genomequebec.com>

B<Joel Fillon> - I<joel.fillon@mcgill.ca>

=cut

# Strict Pragmas
#---------------
use strict qw(vars subs);
use warnings;
#---------------------

# Add the mugqic_pipeline/lib/ path relative to this Perl script to @INC library search variable
use FindBin;
use lib "$FindBin::Bin/../../lib";

# Dependency modules
#-------------------
use Cwd 'abs_path';
use File::Basename;
use Getopt::Std;
use GqSeqUtils;
use LoadConfig;
use Metrics;
use Picard;
use Pipeline;
use Trimmomatic;
use Trinity;
use Version;

# Steps array: each step is run globally or per read set (bam/paired fastq/single fastq) and has a list of parent steps defining step job dependencies
#-------------
my @A_steps = (
  {
    'name'   => 'bamToFastq',
    'loop'   => 'readSet',
    'parentSteps' => []
  },
  {
    'name'   => 'trim',
    'loop'   => 'readSet',
    'parentSteps' => ['bamToFastq']
  },
  {
    'name'   => 'normalization',
    'loop'   => 'readSet',
    'parentSteps' => ['trim']
  },
  {
    'name'   => 'normalizationMergeResults',
    'loop'   => 'global',
    'parentSteps' => ['normalization']
  },
  {
    'name'   => 'trinity',
    'loop'   => 'global',
    'parentSteps' => ['normalizationMergeResults']
  },
  {
    'name'   => 'blastSplitQuery',
    'loop'   => 'global',
    'parentSteps' => ['trinity']
  },
  {
    'name'   => 'blast',
    'loop'   => 'global',
    'parentSteps' => ['blastSplitQuery']
  },
  {
    'name'   => 'blastMergeResults',
    'loop'   => 'global',
    'parentSteps' => ['blast']
  },
  {
    'name'   => 'rsemPrepareReference',
    'loop'   => 'global',
    'parentSteps' => ['trinity']
  },
  {
    'name'   => 'rsem',
    'loop'   => 'sample',
    'parentSteps' => ['rsemPrepareReference']
  },
  {
    'name'   => 'differentialGeneExpression',
    'loop'   => 'global',
    'parentSteps' => ['rsem', 'blastMergeResults']
  },
  {
    'name'   => 'metrics',
    'loop'   => 'global',
    'parentSteps' => ['trim']
  },
  {
    'name'   => 'deliverable',
    'loop'   => 'global',
    'parentSteps' => ['metrics', 'differentialGeneExpression']
  }
);

# Create step hash indexed by step name for easy retrieval
my %H_steps =  map {$_->{'name'} => $_} @A_steps;

# Global variables passed as script options
my $configFile;
my $sampleFile;
my $designFile;
my $pipeline;

# Main call
main();

# General pipeline functions
#---------------------------
sub debug {
  my $message = shift;

  my $debug = 1;    # Set to 1 to display debug messages, 0 to keep output silent
  
  $debug and print STDERR "[DEBUG] $message\n";
}

sub getUsage {
  my $usage = <<END;
MUGQIC Pipeline RNA-Seq De Novo Assembly Version: $Version::version

Usage: perl $0 -h | -c CONFIG_FILE -r STEP_RANGE -s SAMPLE_FILE [-d DESIGN_FILE] [-o OUTPUT_DIR]
  -h  help and usage
  -c  .ini config file
  -r  step range e.g. "1-5", "3,6,7", "2,4-8"
  -s  sample file
  -d  design file
  -o  output directory (default: current)

Steps:
END

  # List and number step names
  for (my $i = 1; $i <= @A_steps; $i++) {
    $usage .= $i . "- " . $A_steps[$i - 1]->{'name'} . "\n";
  }

  return $usage;
}

sub main {
  # Check options
  my %opts;
  getopts('hc:r:s:d:o:', \%opts);

  if (defined($opts{'h'}) ||
     !defined($opts{'c'}) ||
     !defined($opts{'r'}) ||
     !defined($opts{'s'})) {
    die (getUsage());
  }

  # Assign options
  my $stepRange = $opts{'r'};
  my $outputDirectory = $opts{'o'};
  $configFile = $opts{'c'};
  $sampleFile = $opts{'s'};
  $designFile = $opts{'d'};

  # Get config values
  unless (defined $configFile) {die "Error: configuration file is not defined! (use -c option)\n" . getUsage()};
  unless (-f $configFile) {die "Error: configuration file $configFile does not exist!\n" . getUsage()};
  my %cfg = LoadConfig->readConfigFile($configFile);

  $pipeline = Pipeline->new(\@A_steps, $sampleFile, $outputDirectory);

  # Go through steps and create global or read-set jobs accordingly
  foreach my $step ($pipeline->getStepsByRange($stepRange)) {
    my $stepName = $step->getName();
    debug "main: processing step $stepName";

    # ReadSet step creates 1 job per readSet per sample
    if ($step->getLoop() eq 'readSet') {
      foreach my $sample (@{$pipeline->getSamples()}) {
        foreach my $readSet (@{$sample->getReadSets()}) {
          debug "main: processing read set " . $readSet->getName();

          my $rO_job = &$stepName(\%cfg, $readSet);
          if ($rO_job) {
            $rO_job->setLoopTags([$sample->getName(), $readSet->getName()]);
            debug "main: readSet job " . join(".", @{$rO_job->getLoopTags()});
            $step->submitJob(\%cfg, $rO_job);
          }
        }
      }
    # Sample step creates 1 job per sample
    } elsif ($step->getLoop() eq 'sample') {
      foreach my $sample (@{$pipeline->getSamples()}) {
        debug "main: processing sample " . $sample->getName();

        my $rO_job = &$stepName(\%cfg, $sample);
        if ($rO_job) {
          $rO_job->setLoopTags([$sample->getName()]);
          debug "main: sample job " . join(".", @{$rO_job->getLoopTags()});
          $step->submitJob(\%cfg, $rO_job);
        }
      }
    # Global step creates 1 job only
    } else {
      my $rO_job = &$stepName(\%cfg);
      if ($rO_job) {
        $rO_job->setLoopTags([]);
        $step->submitJob(\%cfg, $rO_job);
      }
    }
  }
}

# Step functions
#---------------

sub bamToFastq {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $rO_job;

  if ($readSet->getBAM() and not($readSet->getFASTQ1())) {
    if ($readSet->getRunType() eq "PAIRED_END") {
      $readSet->setFASTQ1($readSet->getBAM());
      $readSet->getFASTQ1() =~ s/\.bam$/.pair1.fastq.gz/;
      $readSet->setFASTQ2($readSet->getBAM());
      $readSet->getFASTQ2() =~ s/\.bam$/.pair2.fastq.gz/;
      $rO_job = Picard::samToFastq($rH_cfg, $readSet->getBAM(), $readSet->getFASTQ1(), $readSet->getFASTQ2());
    } elsif ($readSet->getRunType() eq "SINGLE_END") {
      $readSet->setFASTQ1($readSet->getBAM());
      $readSet->getFASTQ1() =~ s/\.bam$/.single.fastq.gz/;
      $rO_job = Picard::samToFastq($rH_cfg, $readSet->getBAM(), $readSet->getFASTQ1());
    }
  }
  return $rO_job;
}

sub trim {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";

  my $pairedOutput1;
  my $unpairedOutput1;
  my $pairedOutput2;
  my $unpairedOutput2;
  my $singleOutput;

  if ($readSet->getRunType() eq "PAIRED_END") {
    $pairedOutput1 = $trimFilePrefix . "pair1.fastq.gz";
    $unpairedOutput1 = $trimFilePrefix . "single1.fastq.gz";
    $pairedOutput2 = $trimFilePrefix . "pair2.fastq.gz";
    $unpairedOutput2 = $trimFilePrefix . "single2.fastq.gz";
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    $singleOutput = $trimFilePrefix . "single.fastq.gz";
  }

  return Trimmomatic::trim(
    $rH_cfg,
    $readSet->getFASTQ1(),
    $readSet->getFASTQ2(),
    $pairedOutput1,
    $unpairedOutput1,
    $pairedOutput2,
    $unpairedOutput2,
    $singleOutput,
    $readSet->getQualityOffset(),
    $trimFilePrefix . "out",
    $trimFilePrefix . "stats.csv"
  );
}

sub normalization {
  my $rH_cfg = shift;
  my $readSet = shift;

  my $trimFilePrefix = "\$WORK_DIR/trim/" . $readSet->getSample()->getName() . "/" . $readSet->getName() . ".trim.";
  my $normalizationDirectory = "\$WORK_DIR/normalization/" . $readSet->getSample()->getName() . "/" . $readSet->getName();

  if ($readSet->getRunType() eq "PAIRED_END") {
    return Trinity::normalize_by_kmer_coverage(
      $rH_cfg,
      [$trimFilePrefix . "pair1.fastq.gz"],
      [$trimFilePrefix . "pair2.fastq.gz"],
      undef,
      $normalizationDirectory
    );
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    return Trinity::normalize_by_kmer_coverage(
      $rH_cfg,
      undef,
      undef,
      $trimFilePrefix . "single.fastq.gz",
      $normalizationDirectory
    );
  }
}

sub normalizationMergeResults {
  my $rH_cfg = shift;

  my $kmerSize = LoadConfig::getParam($rH_cfg, 'normalization', 'kmerSize', 1, 'int');
  my $maxCoverage = LoadConfig::getParam($rH_cfg, 'normalization', 'maxCoverage', 1, 'int');
  my $maxPctStdev = LoadConfig::getParam($rH_cfg, 'normalization', 'maxPctStdev', 1, 'float');

  my $normalizedFileSuffix = ".normalized_K" . $kmerSize . "_C" . $maxCoverage . "_pctSD" . $maxPctStdev . ".fq";

  my $rA_leftNormalizedFiles;
  my $rA_rightNormalizedFiles;
  my $singleNormalizedFile;

  # Retrieve single/paired end normalized files for each readSet
  foreach my $sample (@{$pipeline->getSamples()}) {
    foreach my $readSet (@{$sample->getReadSets()}) {
      if ($readSet->getRunType() eq "PAIRED_END") {
        my $normalizationDirectory = "\$WORK_DIR/normalization/" . $sample->getName() . "/" . $readSet->getName();
        push(@$rA_leftNormalizedFiles, "$normalizationDirectory/left$normalizedFileSuffix");
        push(@$rA_rightNormalizedFiles, "$normalizationDirectory/right$normalizedFileSuffix");

      } elsif ($readSet->getRunType() eq "SINGLE_END") {
        # Assume only 1 SINGLE_END read file (multiple SINGLE_END read files not supported)
        die "Error: normalizationMergeResults for multiple SINGLE_END reads is not implemented!";
      }
    }
  }

  return Trinity::normalize_by_kmer_coverage(
    $rH_cfg,
    $rA_leftNormalizedFiles,
    $rA_rightNormalizedFiles,
    $singleNormalizedFile,
    "\$WORK_DIR/normalization/global"
  );
}

sub trinity {
  my $rH_cfg = shift;

  if (scalar(@{$pipeline->getRunTypes()}) > 1) {
    die "Error in trinity: assembly of multiple run types (" . join(", ", @{$pipeline->getRunTypes()}) . ") is not supported!";
  }
  my $runType = @{$pipeline->getRunTypes()}[0];

  my $maxCoverage = LoadConfig::getParam($rH_cfg, 'normalization', 'maxCoverage', 1, 'int');
  my $kmerSize = LoadConfig::getParam($rH_cfg, 'normalization', 'kmerSize', 1, 'int');
  my $maxPctStdev = LoadConfig::getParam($rH_cfg, 'normalization', 'maxPctStdev', 1, 'float');

  my $normalizedFileSuffix = ".normalized_K" . $kmerSize . "_C" . $maxCoverage . "_pctSD" . $maxPctStdev . ".fq";

  my $normalizationDirectory = "\$WORK_DIR/normalization/global";
  my $trinityDirectory = "\$WORK_DIR/trinity_out_dir";

  if ($runType eq "PAIRED_END") {
    return Trinity::trinity(
      $rH_cfg,
      ["$normalizationDirectory/left$normalizedFileSuffix"],
      ["$normalizationDirectory/right$normalizedFileSuffix"],
      undef,
      $trinityDirectory
    );
  } elsif ($runType eq "SINGLE_END") {
    return Trinity::trinity(
      $rH_cfg,
      undef,
      undef,
      ["$normalizationDirectory/single$normalizedFileSuffix"],
      $trinityDirectory
    );
  }
}

sub blastSplitQuery {
  my $rH_cfg = shift;

  my $trinityFastaFile = "\$WORK_DIR/trinity_out_dir/Trinity.fasta";
  my $trinityIndexFile = "\$WORK_DIR/trinity_out_dir/Trinity.idx";
  my $reducedTrinityFastaFile = "\$WORK_DIR/trinity_out_dir/Trinity.longest_transcript.fasta";
  my $numJobs = LoadConfig::getParam($rH_cfg, 'blast', 'blastJobs', 1, 'int');

  my $rO_job = new Job([$trinityFastaFile], ["\$WORK_DIR/blast/chunks/Trinity.longest_transcript.fasta_chunk_" . sprintf("%07d", $numJobs - 1)]);
  if (!$rO_job->isUp2Date()) {
    $rO_job->addModules($rH_cfg, [
      ['blast', 'moduleVersion.exonerate']
    ]);

    # Remove previous Trinity assembly FASTA index if present
    my $command = "rm -f $trinityIndexFile && \\\n";
    # Create Trinity assembly FASTA index
    $command .= "fastaindex $trinityFastaFile $trinityIndexFile && \\\n";
    # Create Trinity assembly FASTA subset with longest transcript per component only
    $command .= "fastalength $trinityFastaFile | perl -pe 's/ ((\\S+)_seq\\S+)/\\t\\1\\t\\2/' | sort -k3,3 -k1,1gr | uniq -f2 | cut -f2 | fastafetch $trinityFastaFile -i $trinityIndexFile -q stdin > $reducedTrinityFastaFile && \\\n";

    # Split Trinity assembly FASTA into chunks for BLAST parallelization
    my $chunkDir = "\$WORK_DIR/blast/chunks";
    $command .= "mkdir -p $chunkDir && \\\n";
    $command .= "fastasplit -f $reducedTrinityFastaFile -o $chunkDir -c " . LoadConfig::getParam($rH_cfg, 'blast', 'blastJobs', 1, 'int');

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

sub blast {
  my $rH_cfg = shift;

  my $numJobs = LoadConfig::getParam($rH_cfg, 'blast', 'blastJobs', 1, 'int');

  for (my $jobIndex = 0; $jobIndex < $numJobs; $jobIndex++) {
    # fastasplit creates FASTA chunk files numbered with 7 digits and padded with leading 0s
    my $chunkIndex = sprintf("%07d", $jobIndex);

    my $rO_job = new Job();
    if (!$rO_job->isUp2Date()) {
      $rO_job->addModules($rH_cfg, [
        ['blast', 'moduleVersion.tools'],
        ['blast', 'moduleVersion.exonerate'],
        ['blast', 'moduleVersion.blast']
      ]);

      my $cores = LoadConfig::getParam($rH_cfg, 'blast', 'blastCPUperJob', 1, 'int');
      my $program = LoadConfig::getParam($rH_cfg, 'blast', 'blastProgram');
      my $db = LoadConfig::getParam($rH_cfg, 'blast', 'blastDb');

      # Check if BLAST db files are available
      my $blastDbHome = "\$MUGQIC_INSTALL_HOME/genomes/blast_db";
      `ls $blastDbHome/$db.*[np]hr` or die "Error: $db BLAST db files do not exist in $blastDbHome!";

      my $options = LoadConfig::getParam($rH_cfg, 'blast', 'blastOptions');
      my $chunkDir = "\$WORK_DIR/blast/chunks";
      my $chunkQuery = "$chunkDir/Trinity.longest_transcript.fasta_chunk_$chunkIndex";
      my $chunkResult = "$chunkDir/$program" . "_Trinity.longest_transcript_$db" . "_chunk_$chunkIndex.tsv";

      # Each FASTA chunk is further divided in subchunk per CPU per job as a second level of BLAST parallelization
      # The user must adjust BLAST configuration to optimize num. jobs vs num. CPUs per job, depending on the cluster
      my $command = "parallelBlast.pl -file $chunkQuery --OUT $chunkResult -n $cores --BLAST \\\"$program -db $db $options\\\"";

      $rO_job->addCommand($command);
    }
    $rO_job->setStep($H_steps{"blast"});
    $rO_job->setLoopTags(["chunk_$jobIndex"]);
    $pipeline->getStepByName("blast")->submitJob($rH_cfg, $rO_job);
  }
}

sub blastMergeResults {
  my $rH_cfg = shift;

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {
    my $program = LoadConfig::getParam($rH_cfg, 'blast', 'blastProgram');
    my $db = LoadConfig::getParam($rH_cfg, 'blast', 'blastDb');
    my $blastDir = "\$WORK_DIR/blast";
    my $chunkResults = "$blastDir/chunks/$program" . "_Trinity.longest_transcript_$db" . "_chunk_*.tsv";
    my $result = "$blastDir/$program" . "_Trinity.longest_transcript_$db.tsv";

    # All BLAST chunks are merged into one file named after BLAST program and reference database
    my $command = "cat $chunkResults > $result.tmp && \\\n";
    # Remove all comment lines except "Fields" one which is placed as first line
    $command .= "cat <(grep -m1 '^# Fields' $result.tmp) <(grep -v '^#' $result.tmp) > $result && \\\n";
    $command .= "rm $result.tmp && \\\n";

    # Create a BLAST results ZIP file for future deliverables
    $command .= "zip -j $result.zip $result";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

# The RSEM reference assembly is created once only, and then used by all RSEM read-set jobs in parallel
sub rsemPrepareReference {
  my $rH_cfg = shift;

  return Trinity::rsemPrepareReference($rH_cfg, "\$WORK_DIR/trinity_out_dir/Trinity.fasta");
}

# RSEM abundance estimation is performed by read-set (there should be 1 read-set / sample for RNA-Seq)
sub rsem {
  my $rH_cfg = shift;
  my $sample = shift;

  # Trinity wrapper around RSEM only support one readset as input as of 2013-11-10
  ($sample->getNbReadSets() == 1) or die "Error in rsem: one and only one readset per sample is supported!";

  my $readSet = @{$sample->getReadSets()}[0];
  my $trimFilePrefix = "\$WORK_DIR/trim/" . $sample->getName() . "/" . $readSet->getName() . ".trim.";

  my $rA_input1;
  my $rA_input2;

  if ($readSet->getRunType() eq "PAIRED_END") {
    $rA_input1 = [$trimFilePrefix . "pair1.fastq.gz"];
    $rA_input2 = [$trimFilePrefix . "pair2.fastq.gz"];
  } elsif ($readSet->getRunType() eq "SINGLE_END") {
    $rA_input1 = [$trimFilePrefix . "single.fastq.gz"];
  } else {
    die "Error in rsem: unknown read type\n";
  }

  return Trinity::rsem(
    $rH_cfg,
    "\$WORK_DIR/trinity_out_dir/Trinity.fasta",
    $rA_input1,
    $rA_input2,
    $sample->getName(),
    "\$WORK_DIR/rsem/" . $sample->getName()
  );
}

sub differentialGeneExpression {
  my $rH_cfg = shift;

  my $rO_job = new Job();
  if (!$rO_job->isUp2Date()) {

    # Check design file
    unless (defined $designFile) {die "Error: design file is not defined! (use -d option)\n" . getUsage()};
    unless (-f $designFile) {die "Error: design file $designFile does not exist!\n" . getUsage()};
    $designFile = abs_path($designFile);

    # Retrieve BLAST result file
    my $program = LoadConfig::getParam($rH_cfg, 'blast', 'blastProgram');
    my $db = LoadConfig::getParam($rH_cfg, 'blast', 'blastDb');
    my $blastDir = "\$WORK_DIR/blast";
    my $blastResult = "$blastDir/$program" . "_Trinity.longest_transcript_$db.tsv";

    my $dgeDir = "\$WORK_DIR/DGE";
    my $isoformsMatrix = "$dgeDir/isoforms.counts.matrix";
    my $isoformsAnnotatedMatrix = "$dgeDir/isoforms.counts.$db.matrix";
    my $genesMatrix = "$dgeDir/genes.counts.matrix";
    my $genesAnnotatedMatrix = "$dgeDir/genes.counts.$db.matrix";

    $rO_job->addModules($rH_cfg, [
      ['differentialGeneExpression', 'moduleVersion.trinity'],
      ['differentialGeneExpression', 'moduleVersion.cranR'],
      ['differentialGeneExpression', 'moduleVersion.tools']
    ]);

    my $command = "mkdir -p $dgeDir && \\\n";

    # Create isoforms and genes matrices with counts of RNA-seq fragments per feature using Trinity RSEM utility
    $command .= "merge_RSEM_frag_counts_single_table.pl " . join(" ", map("\$WORK_DIR/rsem/" . $_->getName() . "/" . $_->getName() . ".isoforms.results", @{$pipeline->getSamples()})) . " > $isoformsMatrix && \\\n";
    $command .= "merge_RSEM_frag_counts_single_table.pl " . join(" ", map("\$WORK_DIR/rsem/" . $_->getName() . "/" . $_->getName() . ".genes.results", @{$pipeline->getSamples()})) . " > $genesMatrix && \\\n";

    # Extract isoforms and genes length values from any one of RSEM's results files
    my $sampleName = @{$pipeline->getSamples()}[0]->getName();
    $command .= "cut -f 1,3,4 \$WORK_DIR/rsem/$sampleName/$sampleName.isoforms.results > \$WORK_DIR/rsem/isoforms.lengths.tsv && \\\n";
    $command .= "cut -f 1,3,4 \$WORK_DIR/rsem/$sampleName/$sampleName.genes.results > \$WORK_DIR/rsem/genes.lengths.tsv && \\\n";

    # Merge isoforms and genes matrices with BLAST annotations if any:
    # edger.R requires a matrix with gene/isoform annotation as second column 
    # Keep BLAST best hit only
    # Remove from column headers ".(genes|isoforms).results" created by RSEM
    $command .= "grep -v '^#' $blastResult | awk '!x[\\\$1]++' | awk -F\\\"\\t\\\" 'FNR==NR {a[\\\$1]=\\\$2; next}{OFS=\\\"\\t\\\"; if (a[\\\$1]) {print \\\$1, a[\\\$1]} else {print \\\$1, \\\$1}}' - $isoformsMatrix | sed '1s/^\\t/Isoform\\tSymbol/' | paste - <(cut -f 2- $isoformsMatrix) | sed '1s/\\.isoforms\\.results//g' > $isoformsAnnotatedMatrix && \\\n";
    # Remove "_seq" from isoform BLAST query name and keep BLAST isoform best hit as BLAST gene best hit
    $command .= "grep -v '^#' $blastResult | awk '!x[\\\$1]++' | awk -F\\\"\\t\\\" 'FNR==NR {sub(/_seq.*/, \\\"\\\", \\\$1); a[\\\$1]=\\\$2; next}{OFS=\\\"\\t\\\"; if (a[\\\$1]) {print \\\$1, a[\\\$1]} else {print \\\$1, \\\$1}}' - $genesMatrix | sed '1s/^\\t/Gene\\tSymbol/' | paste - <(cut -f 2- $genesMatrix) | sed '1s/\\.genes\\.results//g' > $genesAnnotatedMatrix && \\\n";

    # Perform edgeR
    $command .= "Rscript \\\$R_TOOLS/edger.R -d $designFile -c $isoformsAnnotatedMatrix -o $dgeDir/isoforms_$db && \\\n";
    $command .= "Rscript \\\$R_TOOLS/edger.R -d $designFile -c $genesAnnotatedMatrix -o $dgeDir/genes_$db && \\\n";

    # Perform DESeq
    $command .= "Rscript \\\$R_TOOLS/deseq.R -d $designFile -c $isoformsAnnotatedMatrix -o $dgeDir/isoforms_$db && \\\n";
    $command .= "Rscript \\\$R_TOOLS/deseq.R -d $designFile -c $genesAnnotatedMatrix -o $dgeDir/genes_$db && \\\n";

    # Merge edgeR results with gene/isoform length values and BLAST description
    $command .= "for gi in genes isoforms; do for f in $dgeDir/\\\${gi}_$db/*/dge_results.csv; do sed '1s/gene_symbol/$db.id/' \\\$f | awk -F\\\"\\t\\\" 'FNR==NR {a[\\\$1]=\\\$2\\\"\\t\\\"\\\$3; next}{OFS=\\\"\\t\\\"; if (a[\\\$1]) {print \\\$0, a[\\\$1]} else {print \\\$0, \\\"\\\", \\\"\\\"}}' \$WORK_DIR/rsem/\\\${gi}.lengths.tsv - | sed '1s/\\t\\\$/length\\teffective_length/' | awk -F\\\"\\t\\\" 'FNR==NR {a[\\\$2]=\\\$NF; next}{OFS=\\\"\\t\\\"; if (a[\\\$2]) {print \\\$0, a[\\\$2]} else {print \\\$0, \\\"\\\"}}' <(grep -v '^#' $blastResult) - | sed '1s/\\\$/description/' > \\\${f/.csv/_$db.csv}; done; done";

    $rO_job->addCommand($command);
  }
  return $rO_job;
}

# Merge all sample Trimmomatic results
sub metrics {
  my $rH_cfg = shift;

  if (scalar(@{$pipeline->getRunTypes()}) > 1) {
    die "Error in metrics: metrics of multiple run types (" . join(", ", @{$pipeline->getRunTypes()}) . ") is not supported!";
  }
  my $runType = @{$pipeline->getRunTypes()}[0];
  my $trimDirectory = "\$WORK_DIR/trim";
  my $pattern = "trim.stats.csv";
  my $outputFile = "\$WORK_DIR/metrics/trimming.stats";

  # Merge all sample Trimmomatic results
  return Metrics::mergeTrimmomaticStats($rH_cfg, $runType, $pattern, $trimDirectory, $outputFile);
}

sub deliverable {
  my $rH_cfg = shift;

  return GqSeqUtils::clientReport($rH_cfg, abs_path($configFile), "\$WORK_DIR", "RNAseqDeNovo");
}