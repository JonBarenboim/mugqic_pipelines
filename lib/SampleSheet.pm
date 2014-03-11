#!/usr/bin/env perl

=head1 NAME

I<SampleSheet>

=head1 SYNOPSIS

SampleSheet->parseSampleSheet(sampleSheet_file_name)

=head1 DESCRIPTION

B<SampleSheet> is a library that parses a Nanuq generated
sample sheet and populates an array with he parsed values.
Each row from the sheet is a new hash in the array.

Input = /path/samplesheet_file_name

Output = %hash


=head1 AUTHOR

B<Louis Letourneau> - I<louis.letourneau@mail.mcgill.ca>
B<Maxime Caron> - I<max.caron@mail.mcgill.ca>

=head1 DEPENDENCY

B<Pod::Usage> Usage and help output.

B<File::Basename> path parsing

B<Cwd> path parsing

=cut

package SampleSheet;

# Strict Pragmas
#---------------------
use strict;
use warnings;
#---------------------

# Dependencies
#--------------------
use File::Basename;
use Text::CSV;
use Cwd 'abs_path';
#--------------------


# SUB
#--------------------

sub parseSampleSheetAsHash {
  my $fileName = shift;
  my $rawReadFormat = shift;

  my $rA_SampleLaneInfos = parseSampleSheet($fileName, $rawReadFormat);
  my %sampleInfo;
  for my $rH_Sample (@$rA_SampleLaneInfos) {
    if (!defined $sampleInfo{$rH_Sample->{'name'}}) {
      $sampleInfo{$rH_Sample->{'name'}} = [];
    }

    push(@{$sampleInfo{$rH_Sample->{'name'}}}, $rH_Sample);
  }
  return \%sampleInfo;
}

sub parseSampleSheetAsHashByProcessingId {
  my $fileName = shift;

  my $rA_SampleLaneInfos = parseSampleSheet($fileName);
  my %sampleInfo;
  for my $rH_Sample (@$rA_SampleLaneInfos) {
    if (!defined($rH_Sample->{'processingSheetId'})) {
      die "Missing processingSheetId";
    }
    if(!defined $sampleInfo{ $rH_Sample->{'processingSheetId'} }) {
      $sampleInfo{ $rH_Sample->{'processingSheetId'} } = [];
    }

    push(@{$sampleInfo{ $rH_Sample->{'processingSheetId'} }}, $rH_Sample);
  }
  return \%sampleInfo;
}

sub parsePairedSampleSheet {
  my $fileName = shift;

  my @retVal;
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";

  while (<SAMPLE_SHEET>) {
    chomp;

    my $line = $_;
    if ($line =~ /^#/) {
      next;
    }

    $line =~ s/"//g;
    my @values = split(/,/, $line);

    my %sampleInfo;
    $sampleInfo{'sample'} = $values[0];
    $sampleInfo{'normal'} = $values[1];
    $sampleInfo{'tumor'} = $values[2];
    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

sub parseSampleSheet {
  my $fileName = shift;
  my $rawReadFormat = shift;

  my @retVal;
  my $csv = Text::CSV->new();
  open(SAMPLE_SHEET, "$fileName") or die "Can't open $fileName\n";
  my $line = <SAMPLE_SHEET>;
  $csv->parse($line);
  my @headers = $csv->fields();
  my ($nameIdx, $libraryBarcodeIdx, $runIdIdx, $laneIdx, $runTypeIdx, $statusIdx, $qualOffsetIdx, $bedFilesIdx, $processingSheetIdIdx, $libSourceIdx) = parseHeaderIndexes(\@headers);

  while($line = <SAMPLE_SHEET>) {
    $csv->parse($line);
    my @values = $csv->fields();
    if ($values[$statusIdx] =~ /invalid/) {
      warn "Invalid: $values[$nameIdx] $values[$runIdIdx] $values[$laneIdx]\n";
      next;
    }

    my %sampleInfo;
    $sampleInfo{'name'} = $values[$nameIdx];
    $sampleInfo{'libraryBarcode'} = $values[$libraryBarcodeIdx];
    $sampleInfo{'runId'} = $values[$runIdIdx];
    $sampleInfo{'lane'} = $values[$laneIdx];
    $sampleInfo{'runType'} = $values[$runTypeIdx];
    $sampleInfo{'qualOffset'} = $values[$qualOffsetIdx];
    my @bedFiles = undef;
    if ($bedFilesIdx != -1) {
        my @bedFiles = split(';', $values[$bedFilesIdx]);
    }
    $sampleInfo{'bedFiles'} = \@bedFiles;
    if ($processingSheetIdIdx > -1) {
      $sampleInfo{'processingSheetId'} = $values[$processingSheetIdIdx];
    }
    $sampleInfo{'libSource'} = $values[$libSourceIdx];

    my $rawReadPrefix = $sampleInfo{'name'} . "." . $sampleInfo{'libraryBarcode'} . "." . $sampleInfo{'qualOffset'} . ".";

    if (defined($rawReadFormat) and $rawReadFormat eq "fastq") {
      if ($values[$runTypeIdx] eq "PAIRED_END") {
        $sampleInfo{'read1File'} = $rawReadPrefix . "pair1.fastq.gz";
        $sampleInfo{'read2File'} = $rawReadPrefix . "pair2.fastq.gz";
      } elsif ($values[$runTypeIdx] eq "SINGLE_END") {
        $sampleInfo{'read1File'} = $rawReadPrefix . "single.fastq.gz";
      } else {
        die "Error in SampleSheet::parseSampleSheet: unrecognized run type $values[$runTypeIdx]";
      }
    } else {    # BAM format by default
      $sampleInfo{'read1File'} = $rawReadPrefix . "bam";
    }

    push(@retVal, \%sampleInfo);
  }

  return \@retVal;
}

sub parseHeaderIndexes {
  my $rA_headers = shift;
  my $nameIdx=-1;
  my $libraryBarcodeIdx=-1;
  my $runIdIdx=-1;
  my $laneIdx=-1;
  my $runTypeIdx=-1;
  my $statusIdx=-1;
  my $qualOffsetIdx=-1;
  my $bedFilesIdx=-1;
  my $processingSheetIdIdx=-1;
  my $libSourceIdx = -1;
	
  for(my $idx=0; $idx < @{$rA_headers}; $idx++) {
    my $header = $rA_headers->[$idx];
    $header =~ s/"//g;
    if ($header eq "Name") {
      $nameIdx=$idx;
    } elsif ($header eq "Library Barcode") {
      $libraryBarcodeIdx=$idx;
    } elsif ($header eq "Run") {
      $runIdIdx=$idx;
    } elsif ($header eq "Region") {
      $laneIdx=$idx;
    } elsif ($header eq "Run Type") {
      $runTypeIdx=$idx;
    } elsif ($header eq "Status") {
      $statusIdx=$idx;
    } elsif ($header eq "Quality Offset") {
      $qualOffsetIdx=$idx;
    } elsif ($header eq "BED Files") {
      $bedFilesIdx=$idx;
    }
    elsif($header eq "ProcessingSheetId") {
      $processingSheetIdIdx=$idx;
    }
    elsif($header eq "Library Source") {
      $libSourceIdx=$idx;
    }
  }

  my $sampleSheetErrors="";
  my $sampleSheetWarnings="";
  if ($nameIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Sample Name\n";
  }
  if ($libraryBarcodeIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Library Barcode\n";
  }
  if ($runIdIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Run ID\n";
  }
  if ($laneIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Lane\n";
  }
  if ($runTypeIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Run Type\n";
  }
  if ($statusIdx==-1) {
    $sampleSheetErrors .= "[Error] Missing Status\n";
  }
  if($qualOffsetIdx==-1) {
      $sampleSheetErrors.="Missing Quality Offset\n";
  }
  if($bedFilesIdx==-1) {
    $sampleSheetWarnings.="[Warning] Missing BED Files\n";
  }
  if($libSourceIdx==-1) {
    $sampleSheetErrors.="Missing Library Source\n";
  }
  if(length($sampleSheetWarnings) > 0) {
    warn $sampleSheetWarnings;
  }
  if (length($sampleSheetErrors) > 0) {
    die $sampleSheetErrors;
  }

  return ($nameIdx, $libraryBarcodeIdx, $runIdIdx, $laneIdx, $runTypeIdx, $statusIdx, $qualOffsetIdx, $bedFilesIdx, $processingSheetIdIdx, $libSourceIdx);
}
1;
