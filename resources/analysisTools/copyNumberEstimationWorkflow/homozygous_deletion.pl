#!/usr/bin/env perl

use strict;
use warnings;
use v5.10;
use List::Util;

#This script adds crest information to segments which have start or end in commmon with defined crest segment
#and adds a mappablilty, homoyyzgDel, unmappable tag depending on the mappability value


my %opts;

use Getopt::Long;

#GetOptions( 'a=s' => \$opts{file},
#            'b=s' => \$opts{crest_file},
#            'c=s' => \$opts{outfile},
#            'd=s' => \$opts{min_seg_map},
#          );
GetOptions( 'a=s' => \$opts{crest_file},
            'b=s' => \$opts{min_seg_map},
          );
          
my %fields;
my @cols;
my $header;
my $crests;
my @crestCols = qw(chr start end length type chr1 chr2 id);
my %bp;
my $start;
my $end;

#open (OUT, " | bgzip -c > $outfile") or die "could not open $outfile for writing: $!\n";

say join "\t", "chromosome", "start", "end", "crest", "length", "tcnNbrOfLoci", "tcnMean", "tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci", "map";

my $crest = $opts{crest_file};
open (CREST, "$crest") or die $!;

while (<CREST>){
  chomp;
  $crests = {};
  @{$crests}{@crestCols} = split(/\t/);
  $crests->{type} =~ s/([IC]TX).+/$1/g;
  $bp{"$crests->{chr}".":"."$crests->{start}"}=$crests;
  $bp{"$crests->{chr}".":"."$crests->{end}"}=$crests;
}

#$header = <IN>;
$header = <STDIN>;
chomp $header;
@cols = qw(chromosome tcnID dhId start end tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhMean c1Mean c2Mean map);

#while (<IN>) {
while (<STDIN>) {
  chomp;
#  print $_."\n";
  @fields{@cols} = split(/\t/);
  
  if ($fields{chromosome} eq "NA") {
    next;
  }
  
  $start = "";
  $end   = "";
  
  $fields{start} =~ s/\.0$//;
  $fields{end} =~ s/\.0$//;
 
  $start = "$fields{chromosome}".":"."$fields{start}";
  $end   = "$fields{chromosome}".":"."$fields{end}";

  if (defined($bp{"$start"})){
    $fields{crest} = $bp{"$start"}{type};
  }elsif (defined($bp{"$end"})){
    $fields{crest} = $bp{"$end"}{type};
  }else{
    $fields{crest} = "NA";
  }

  $fields{length} = $fields{end}-$fields{start}+1;

#unmappable
  if ($fields{map} < $opts{min_seg_map}*$fields{length}){
    say join "\t", $fields{chromosome}, $fields{start}, $fields{end}, $fields{crest}, $fields{length}, $fields{tcnNbrOfLoci}, $fields{tcnMean}, $fields{tcnNbrOfSNPs}, $fields{tcnNbrOfHets}, $fields{dhNbrOfLoci}, "unmappable";

#homozygous deletion
  }elsif(($fields{tcnMean} eq "NA") && ($fields{map} > $opts{min_seg_map}*$fields{length})){
    say join "\t", $fields{chromosome}, $fields{start}, $fields{end}, $fields{crest}, $fields{length}, $fields{tcnNbrOfLoci}, 0, $fields{tcnNbrOfSNPs}, $fields{tcnNbrOfHets}, $fields{dhNbrOfLoci}, "homozygousDel";

#mappable
  }elsif(($fields{tcnMean} ne "NA") && ($fields{map} > $opts{min_seg_map}*$fields{length})){
    say join "\t", $fields{chromosome}, $fields{start}, $fields{end}, $fields{crest}, $fields{length}, $fields{tcnNbrOfLoci}, $fields{tcnMean}, $fields{tcnNbrOfSNPs}, $fields{tcnNbrOfHets}, $fields{dhNbrOfLoci}, "mappable";
  }
}
#close OUT or die "could not close $outfile: $!\n";
#system("tabix -s 1 -b 2 -e 3 $outfile") && die "could not tabix $outfile: $!\n";
