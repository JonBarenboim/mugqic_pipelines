#!/usr/local/bin/perl -w

#my ($read_files, $db_type, $map_type, $read_type, $cutoff_type, $cutoff0, $cutoff1, $cutoff2, $cutoff3) = @ARGV;
my ($infile, $lengthfile, $outfile1, $outfile2, $outfile3, $read_files, $db_type, $map_type, $read_type, $cutoff_type, $cutoff0, $cutoff1, $cutoff2, $cutoff3) = @ARGV;

if ($cutoff_type == 0) {
    $cutoff0 = 0;
    $cutoff1 = 0;
    $cutoff2 = 0;
    $cutoff3 = 0;
}
#
#my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
#my $Workpath = "";
#
#my $tmp = $Workpath.$read_files."_".$read_type;
#my $lengthfile = $tmp."_IDs_length.txt";
#if ($read_type eq 'singletons') {
#    my $tmp1 = $Workpath.$read_files."1_".$read_type."_IDs_length.txt";
#    my $tmp2 = $Workpath.$read_files."2_".$read_type."_IDs_length.txt";
#    system("cat $tmp1 $tmp2 > $lengthfile");
#}
#
#my $infile;
#my $outfile1;
#my $outfile2;
#if ($db_type eq 'micro_cds') {
#    $infile = $tmp."_n_micro_cds_sorted.".$map_type."out";
#    $outfile1 = $tmp."_n_micro_cds_".$map_type."_IDs.txt";
#    $outfile2 = $tmp."_n_micro_cds_".$map_type."_pairs.txt";
#    $outfile3 = $tmp."_n_micro_cds_".$map_type."_hitsID.txt";
#} elsif ($db_type eq 'nr') {
#    $infile = $tmp."_nr.".$map_type."out";
#    if ($read_type eq 'singletons') {
#        $infile = $tmp."_nr_sorted.".$map_type."out";
#    }
#    $outfile1 = $tmp."_nr_".$map_type."_IDs.txt";
#    $outfile2 = $tmp."_nr_".$map_type."_pairs.txt";
#    $outfile3 = $tmp."_nr_".$map_type."_hitsID.txt";
#} elsif ($db_type eq 'ecoli_ppi') {
#    if ($read_type eq 'genes') {
#        $infile = $Workpath."microbial_cds_sub_ecoli_ppi.diamondout";
#        $lengthfile = $Workpath."microbial_cds_sub_IDs_length.txt";
#    } elsif ($read_type eq 'proteins') {
#        $infile = $Workpath."nr_all_sub_ecoli_ppi.diamondout";
#        $lengthfile = $Workpath."nr_all_sub_IDs_length.txt";
#    }
#    my ($root, $suff_in) = $infile =~ m/(.*)\.(.*)/;
#    $outfile1 = $root."_IDs.txt";
#    $outfile2 = $root."_pairs.txt";
#    $outfile3 = $root."_hitsID.txt";
#}
#print "Outputfile: $outfile1\n\n";
#print "infile: $infile\n";
#print "length: $lengthfile\n";
#print "outfile1: $outfile1\n";
#print "outfile2: $outfile2\n";
#print "outfile3: $outfile3\n";
#die;

my %reads;
open(INPUT0, $lengthfile) or die "Error opening $lengthfile : $!\n";
while( my $line = <INPUT0> ) {
    chomp($line);
    my @values = split(/\t/, $line);
    if ($read_type eq 'proteins') {
        $reads{$values[1]} = $values[2];
    } elsif ($read_type eq 'genes') {
        $reads{$values[0]} = $values[1];
    } else {
        $reads{$values[0]} = $values[2];
    }
}
close INPUT0;

my %baits;
my %hits;
my $max_score = 0;
my $line1 = 'na';
my $k = 0;
unlink($outfile1);
open(OUTPUT1, '>>', $outfile1) or die "Error opening $outfile1 : $!\n";
unlink($outfile2);
open(OUTPUT2, '>>', $outfile2) or die "Error opening $outfile2 : $!\n";
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while  (my $myline = <INPUT>) {
    chomp($myline);
    my @line = split(/\t/, $myline);

    my $mydbID = $line[0];
    if (($read_files eq 'nr_sub') and (not exists $reads{$mydbID})) {
        $mydbID = join('', $line[0], ',');
    }
    if (($read_files eq 'nr_sub') and (not exists $reads{$mydbID})) {
        $mydbID = join('', $line[0], ';');
    }

    if ($line1 eq $line[0]) {
        $k++;
    }
    else {
        $k = 0;
    }

    if ($k == 0) {
        $line1 = $mydbID;
        $max_score = $line[11];

        my $hitgiID = $line[1];

        my $pident = $line[2];
        $pident = sprintf("%.2f", $pident);

        my $poverlap;
        if (not exists $reads{$mydbID}) {
            $reads{$mydbID} = 100;
        }
        if ($map_type eq 'diamond') {
            $poverlap = 100 * (3 * $line[3]) / $reads{$mydbID};
        } else {
            $poverlap = 100 * $line[3] / $reads{$mydbID};
        }
        $poverlap = sprintf("%.2f", $poverlap);

        my $evalue = $line[10];

        if ($reads{$mydbID} >= $cutoff0) {
            if ($max_score >= $cutoff3) {
                $baits{$mydbID} = 1;
                print OUTPUT1 "$mydbID\n";
                print OUTPUT2 "$mydbID\t$hitgiID\t$pident\t$poverlap\t$evalue\t$max_score\n";
                $hits{$hitgiID} = 1;
            } else {
                $max_score++;
            }
        } else {
            if (($pident >= $cutoff1) and ($poverlap >= $cutoff2)) {
                $baits{$mydbID} = 1;
                print OUTPUT1 "$mydbID\n";
                print OUTPUT2 "$mydbID\t$hitgiID\t$pident\t$poverlap\t$evalue\t$max_score\n";
                $hits{$hitgiID} = 1;
            } else {
                $max_score++;
            }
        }
    } else {
        if ($line[11] == $max_score) {
            my $hitgiID = $line[1];

            my $pident = $line[2];
            $pident = sprintf("%.2f", $pident);

            my $poverlap;
            if (not exists $reads{$mydbID}) {
                $reads{$mydbID} = 100;
            }
            if ($map_type eq 'diamond') {
                $poverlap = 100 * (3 * $line[3]) / $reads{$mydbID};
            } else {
                $poverlap = 100 * $line[3] / $reads{$mydbID};
            }
            $poverlap = sprintf("%.2f", $poverlap);

            my $evalue = $line[10];

            if ($reads{$mydbID} >= $cutoff0) {
                if ($max_score >= $cutoff3) {
                    print OUTPUT2 "$mydbID\t$hitgiID\t$pident\t$poverlap\t$evalue\t$max_score\n";
                    $hits{$hitgiID} = 1;
                }
            } else {
                if (($pident >= $cutoff1) and ($poverlap >= $cutoff2)) {
                    print OUTPUT2 "$mydbID\t$hitgiID\t$pident\t$poverlap\t$evalue\t$max_score\n";
                    $hits{$hitgiID} = 1;
                }
            }
        }
    }
}
close INPUT;
close OUTPUT1;
close OUTPUT2;

print "number of mapped reads:  ".keys( %baits ).".\n";

print "number of hits:  ".keys( %hits ).".\n\n\n";
unlink($outfile3);
open(OUTPUT3, '>>', $outfile3) or die "Error opening $outfile3 : $!\n";
foreach $key (sort keys %hits) {
    print OUTPUT3 $key."\n";
}
close(OUTPUT3);

exit 0;

