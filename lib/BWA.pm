#!/usr/env/perl

=head1 NAME

I<BWA>

=head1 SYNOPSIS

BWA->aln()

=head1 DESCRIPTION

B<BWA> is a library that aligns fastqs on a reference genome

Input = file_name

Output = array


=head1 AUTHOR


=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

=cut

package BWA;

# Strict Pragmas
#--------------------------
use strict;
use warnings;

#--------------------------

# Dependencies
#-----------------------
use LoadConfig;

# SUB
#-----------------------
sub aln {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $single1         = shift;
  my $single2         = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $group           = shift;

    my $command         = "";
    if ( defined($pair1) && defined($pair2) ) {
        $command = pairCommand( $rH_cfg, $sampleName, $pair1, $pair2, $optOutputPrefix, $rgId, $rgSample, $rgLibrary, $rgPlatformUnit, $rgCenter, $group);
    }
    elsif ( defined($single1) ) {
        $command = singleCommand( $rH_cfg, $sampleName, $single1, $optOutputPrefix, $rgId, $rgSample, $rgLibrary, $rgPlatformUnit, $rgCenter, $group);
    }
    else {
        die "Unknown runType, not paired or single\n";
    }

    return $command;
}

sub pairCommand {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $pair1           = shift;
  my $pair2           = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $group           = shift;    # Variable used by deNovoAssembly pipeline

  my $bwaRefIndex = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaRefIndex' );
  if (defined $group){
    $bwaRefIndex = "assembly/" . $group .  $bwaRefIndex;
  }
  
  my $outputSai1Name = $optOutputPrefix.'.pair1.sai';
  my $outputSai2Name = $optOutputPrefix.'.pair2.sai';
  my $outputBAM = $optOutputPrefix.'.sorted.bam';
  my $bamFileDate = -M $outputBAM;

    my @commands;
    # -M gives modified date relative to now. The bigger the older.
    if ( !defined($bamFileDate) || !defined( -M $pair1 ) || !defined( -M $pair2 ) || $bamFileDate > -M $pair1 || $bamFileDate > -M $pair2 ) {
        my $sai1Command = "";
        my $sai2Command = "";
        my $bwaCommand  = "";
        $sai1Command .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' ;';
        $sai1Command .= 'bwa aln';
        $sai1Command .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
        $sai1Command .= ' ' . $bwaRefIndex;
        $sai1Command .= ' ' . $pair1;
        $sai1Command .= ' -f ' . $outputSai1Name;
        push( @commands, $sai1Command );

        $sai2Command .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' ;';
        $sai2Command .= 'bwa aln';
        $sai2Command .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
        $sai2Command .= ' ' . $bwaRefIndex;
        $sai2Command .= ' ' . $pair2;
        $sai2Command .= ' -f ' . $outputSai2Name;
        push( @commands, $sai2Command );

        my $rgTag = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSample . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";
        $bwaCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' ;';
        $bwaCommand .= ' module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.picard').' ;';
        $bwaCommand .= ' bwa sampe ';
        $bwaCommand .= ' '.LoadConfig::getParam( $rH_cfg, 'aln', 'bwaExtraSamXeFlags' );
        $bwaCommand .= ' -r ' . $rgTag;
        $bwaCommand .= ' ' . $bwaRefIndex;
        $bwaCommand .= ' ' . $outputSai1Name;
        $bwaCommand .= ' ' . $outputSai2Name;
        $bwaCommand .= ' ' . $pair1;
        $bwaCommand .= ' ' . $pair2;
        $bwaCommand .= ' | java -Xmx' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortSamRam' ) . ' -jar \${PICARD_HOME}/SortSam.jar INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
        $bwaCommand .= ' OUTPUT=' . $outputBAM;
        $bwaCommand .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortSamRecInRam' );
        push( @commands, $bwaCommand );
    }

    return \@commands;
}

sub singleCommand {
  my $rH_cfg          = shift;
  my $sampleName      = shift;
  my $single          = shift;
  my $optOutputPrefix = shift;
  my $rgId            = shift;
  my $rgSample        = shift;
  my $rgLibrary       = shift;
  my $rgPlatformUnit  = shift;
  my $rgCenter        = shift;
  my $group           = shift;
    
  my $bwaRefIndex = LoadConfig::getParam( $rH_cfg, 'aln', 'bwaRefIndex' );

  if (defined $group){
    $bwaRefIndex = "assembly/" . $group .  $bwaRefIndex;
  }
  
  my $outputSaiName = $optOutputPrefix.'.single.sai';
  my $outputBAM = $optOutputPrefix.'.sorted.bam';
  my $bamFileDate = -M $outputBAM;

    my @commands;
    # -M gives modified date relative to now. The bigger the older.
    if ( !defined($single) || !defined( -M $single ) || !defined($bamFileDate) || $bamFileDate > -M $single ) {
        my $saiCommand = "";
        my $bwaCommand = "";

        $saiCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' ;';
        $saiCommand .= 'bwa aln';
        $saiCommand .= ' -t ' . LoadConfig::getParam( $rH_cfg, 'aln', 'bwaAlnThreads' );
        $saiCommand .= ' ' . $bwaRefIndex;
        $saiCommand .= ' ' . $single;
        $saiCommand .= ' -f ' . $outputSaiName;
        push( @commands, $saiCommand );

        my $rgTag = "'" . '@RG\tID:' . $rgId . '\tSM:' . $rgSample . '\tLB:' . $rgLibrary . '\tPU:run' . $rgPlatformUnit . '\tCN:' . $rgCenter . '\tPL:Illumina' . "'";
        $bwaCommand .= 'module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.bwa').' ;';
        $bwaCommand .= ' module load '.LoadConfig::getParam($rH_cfg, 'aln', 'moduleVersion.picard').' ;';
        $bwaCommand .= ' bwa samse';
        $bwaCommand .= ' '.LoadConfig::getParam( $rH_cfg, 'aln', 'bwaExtraSamXeFlags' );
        $bwaCommand .= ' -r ' . $rgTag;
        $bwaCommand .= ' ' . $bwaRefIndex;
        $bwaCommand .= ' ' . $outputSaiName;
        $bwaCommand .= ' ' . $single;
        $bwaCommand .= ' | java -Xmx' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortSamRam' ) . ' -jar \${PICARD_HOME}/SortSam.jar INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate';
        $bwaCommand .= ' OUTPUT=' . $outputBAM;
        $bwaCommand .= ' MAX_RECORDS_IN_RAM=' . LoadConfig::getParam( $rH_cfg, 'aln', 'sortSamRecInRam' );
        push( @commands, $bwaCommand );
    }

    return \@commands;
}

sub index {
    my $rH_cfg      = shift;
    my $sampleName  = shift;
    my $rH_laneInfo = shift;

    my $laneDirectory = "assembly/" . $sampleName . "/";
    my %retVal;

    my $command = 'module add mugqic/bwa/0.6.2;';
    $command .= ' bwa index ' . $laneDirectory . $rH_cfg->{'aln.bwaRefIndex'};

    $retVal{'command'} = $command;
    return ( \%retVal );

}

1;

