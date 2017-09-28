#!/usr/bin/env perl

# Copyright (c) 2017 The ACEseq workflow developers.
# Distributed under the GNU GENERAL PUBLIC LICENSE (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE_GNU.txt).

use strict;
use warnings;
use v5.10;

use List::Util qw(sum);

use constant DEBUG => 0; # set to 1 and it will also print b-file lines

my %opts;
my @standardColnames;

BEGIN {
  use Getopt::Long;
  %opts = ('tabix_bin' => 'tabix',	  
 #'bgzip_bin' => 'bgzip',
	   'bFileType' => 'vcf',
	   'padding' => 0,
	   'maxBorderDistanceSum' => -1, # negative values in one of these disable this filter
	   'minOverlapFraction' => -1,   # CAVE: it is one filter --- both conditions are connected by OR#	   'reportLevel' => 2 # 1: exact or position or all; 2: exact or all; 3: all
	  );

 GetOptions( 'a=s' => \$opts{afile},
	      'b=s' => \$opts{bfile},
####	      'columnName=s' => \$opts{columnname},
	      'padding:i' => \$opts{padding},
	      'tabix_bin:s' => \$opts{tabix_bin},
	      #'bgzip_bin:s' => \$opts{bgzip_bin},
	      'bAdditionalColumns:s' => \$opts{bcolumns},
	      'minOverlapFraction:f' => \$opts{min_overlap_frac},
	      'maxBorderDistanceSum:i' => \$opts{max_border_distance_sum},
	      'bFileType:s' => \$opts{bFileType},
	      'reportLevel:i' => \$opts{reportLevel}
	    );


#  $opts{afile} = "/ibios/co02/bludau/pscbs/workflow/results/163298/prune/fit_prune.txt";
#  $opts{bfile} = "/ibios/co02/bludau/pscbs/workflow/wgEncodeCrgMapabilityAlign100mer.bedGraph.gz";
#  $opts{bFileType} = "bed";
#  $opts{tabix_bin} = "tabix";

  ### Process user input and check validity
  
  $opts{bFileType} = lc($opts{bFileType});

  my %valid_bFileTypes = ( 'vcf' => 1,
			   'bed' => 1,
			 );
  $valid_bFileTypes{$opts{bFileType}} || die "Only b-files of type ". join(" ", sort keys %valid_bFileTypes). " are supported";



  ### Now make constants
  require constant;
  for (keys %opts) {
    constant->import(uc($_), $opts{$_});
  }
  constant->import('COLUMNNAME', 'MAP');
}

if (BFILETYPE eq 'vcf') {
  @standardColnames = qw(CHROM POS ID REF ALT QUAL FILTER INFO);
} else {
  @standardColnames = qw(chrom chromStart chromEnd name);
}





#my $newcol = $opts{columnname};
#my $pad = $opts{padding};
 
open(A, "$opts{afile}") || die "Could not open a-file $opts{afile}\n";

# guess b file chr format ## TODO: this does not work in every case
my ($b_chr_prefix, $b_chr_suffix);
open(GUESS, TABIX_BIN . " -l $opts{bfile} | ");
while (<GUESS>) {
  chomp;
  if (/([^\d]*)\d+(.*)/) {
    $b_chr_prefix = $1;
    $b_chr_suffix = $2;
    last;
  }
}
close GUESS;
$b_chr_prefix = '' if (! defined($b_chr_prefix));
$b_chr_suffix = '' if (! defined($b_chr_suffix));


# get b-file column names
my @b_colnames;
my $b_header_cmd = TABIX_BIN() . " -h $opts{bfile} $b_chr_prefix" . '1' . $b_chr_suffix . ":0-0 | tail -n1";
my $b_header = `$b_header_cmd`;
if ($b_header =~ s/^\#CHROM/CHROM/ || $b_header =~ s/^\#chrom/chrom/) {
  @b_colnames = split(/\t/);
} else {
  # no header; give standard col names
  @b_colnames =  @standardColnames;
}

#open(HEAD, BGZIP_BIN . " -c -d $opts{bfile} |");
#while (<HEAD>) {
#  if (s/^\#CHROM/CHROM/) {
#    @b_colnames = split(/\t/);
#      last;
#  }
#  if (/^[^\#]/) {
#      # no header; give standard col names
#    @b_colnames =  @standardColnames;
#    warn "b-file $opts{bfile} has no header; standard column names used";
#      last;
#  }
#}

my $newcolidx;
my %a_fields;
my @b_fields;
my @a_columns;
my @b_columns;
@b_columns = split(',',$opts{bcolumns}) if (defined($opts{bcolumns}));
my $chr='';
my $chr_raw = '';

my $header = <A>;
#while ($header = <A>) {
#  last if ($header =~ /^\#CHROM/); # that is the line with the column names
#  print $header; # print out every preceeding line
#}
chomp $header;
#$header =~ s/^\#CHROM/CHROM/;
#@a_columns = split(/\s/,$header);
@a_columns = qw(CHROM tcnId dhId START END tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhMean c1Mean c2Mean);
# my $t_label = $a_columns[9]; # this column position is fixed in the vcf spec so nobody should ever change it
#$newcolidx = 0;
#foreach (@a_columns) {
#  if ($_ eq COLUMNNAME()) {
#    last;
#  }
#  $newcolidx++;
#}
#if ($newcolidx == @a_columns) {
#  say '#', $header, "\t", COLUMNNAME();
push @a_columns, COLUMNNAME();
#} else {
#  say '#', $header;
#}

my @matches;
my $a_line;
my $b_line;
my %b_lines;
my $b_fh;
my $next_b_line;

my ($a_left, $a_right);
my $chr_changed;
my $b_linectr = 0; # need this as hash key
my $b_linenr;
my ($ol_left, $ol_right, $ol_length, $ol_frac);
my ($dist_left, $dist_right);
my $alt;

AFILE_LOOP: while ($a_line=<A>) {
  @matches = ();
  chomp($a_line);
  @a_fields{@a_columns} = split(/\s/, $a_line);
  #say join "\t", map { "$_:$a_fields{$_}" } @a_columns;
  if ($a_fields{CHROM} ne $chr_raw) {
    $chr_changed = 1;
    $chr_raw = $a_fields{CHROM};
    $chr = $chr_raw =~ s/[^\dXY]*([\dXY]+).*/$1/r;
  } else {
    $chr_changed = 0;
  }
  $a_left = $a_fields{START};
  $a_right = $a_fields{END};
 
  #if (defined($a_fields{INFO}) && $a_fields{INFO} =~ /END=(\d+)/) {
  #  $a_right = $1 - 1;
  #} else {
  #  $a_right = $a_left + length($a_fields{REF}) - 1;
  #}
  
  if ($chr_changed) {
    %b_lines = ();
    $next_b_line = {};
    close $b_fh if (ref($b_fh));
    open $b_fh, TABIX_BIN . " $opts{bfile} ${b_chr_prefix}${chr}${b_chr_suffix} |" or die "opening b file $opts{bfile} with tabix failed";
    $b_linectr = 0;
  }
  if (! defined($next_b_line->{left}) || $next_b_line->{left}-PADDING() <= $a_right) {
    # read new b_lines until we have one where the left coordinate is higher than a_right + pad
    while ($b_line=<$b_fh>) {
      if (defined($next_b_line->{left})) {
	$b_lines{$b_linectr} = $next_b_line;
	$next_b_line = {};
      }
      chomp($b_line);
      DEBUG && say $b_line;
      $b_linectr++;
      @b_fields = split(/\t/, $b_line);
      ($b_fields[0] eq $b_chr_prefix.$chr.$b_chr_suffix) || die "Chromosomes between tumor and germline file do not match: $b_fields[0] ne ${b_chr_prefix}${chr}${b_chr_suffix}";
      if (BFILETYPE eq 'vcf') { ### b-file is vcf file
	$next_b_line->{left} =  $b_fields[1];
        $next_b_line->{ref}=$b_fields[3];
        $next_b_line->{alt}=$b_fields[4];
        $next_b_line->{report} = 'POS='.$next_b_line->{left}.';';
        if ($b_fields[7] =~ /END=(\d+)/) {
          $next_b_line->{right} = $1;
        } else {
          $next_b_line->{right} = $next_b_line->{left} + length($next_b_line->{ref});
          $next_b_line->{report} .= 'END='.$next_b_line->{right}.';';
        }
        $next_b_line->{report} .= join ';', (map { $b_colnames[$_] . '=' .$b_fields[$_] } @b_columns), $b_fields[7];
      } else { ### b-file is bed file
	$next_b_line->{left} =  $b_fields[1] + 1;
        $next_b_line->{ref} = 'NA';
        $next_b_line->{alt} = 'NA';
        #$next_b_line->{report} = 'POS='.$next_b_line->{left}.';';
	$next_b_line->{right} = $b_fields[2];
	#$next_b_line->{report} .= 'END=' . $next_b_line->{right} . ';';
        #$next_b_line->{report} .= join ';', (map { $b_colnames[$_] . '=' .$b_fields[$_] } @b_columns);
	$next_b_line->{map} = $b_fields[3];
      }
      last if ($next_b_line->{left}-PADDING() > $a_right);
    }
  }
  # Now compare a_line to all b_lines in b_lines hash
 B_LINE_LOOP: foreach $b_linenr (sort keys %b_lines) {
    if ( $b_lines{$b_linenr}{right}+PADDING() < $a_left ) {
      delete($b_lines{$b_linenr});
      next;
    }
    if ( $b_lines{$b_linenr}{left}-PADDING() > $a_right ) {
      next;
    }
    #if (BFILETYPE eq 'vcf') {
    #  foreach $alt (split(',', $b_lines{$b_linenr}{alt})) { # to check for exact match we must go through all alt alleles in b files
#	if ( $b_lines{$b_linenr}{left} == $a_left && $b_lines{$b_linenr}{right} == $a_right && $b_lines{$b_linenr}{ref} eq $a_fields{REF} && $alt eq $a_fields{ALT} ) {#
#	  push(@matches,'MATCH=exact;' . $b_lines{$b_linenr}{report});
#	  next B_LINE_LOOP ;
#	}
 #     }
  #  }
   # if ( $b_lines{$b_linenr}{left} == $a_left && $b_lines{$b_linenr}{right} == $a_right ) {
    #  push(@matches,'MATCH=position;' . $b_lines{$b_linenr}{report});
    #  next;
    #}
    if ($b_lines{$b_linenr}{left} <= $a_right  && $b_lines{$b_linenr}{right} >= $a_left) {
      $ol_left = ($b_lines{$b_linenr}{left} > $a_left) ? $b_lines{$b_linenr}{left} : $a_left;
      $ol_right = ($b_lines{$b_linenr}{right} < $a_right) ? $b_lines{$b_linenr}{right} : $a_right;
      $ol_length = $ol_right - $ol_left + 1;
      #$ol_frac = $ol_length / (($b_lines{$b_linenr}{right}-$b_lines{$b_linenr}{left} > $a_right-$a_left) ? $b_lines{$b_linenr}{right}-$b_lines{$b_linenr}{left} + 1 : $a_right-$a_left + 1);
      #$ol_frac = sprintf("%.2f", $ol_frac);
      push(@matches, $ol_length * $b_lines{$b_linenr}{map});
    } #else {
      #$ol_frac = 0;
    #}
    #$dist_left = $b_lines{$b_linenr}{left} - $a_left;
    #$dist_right = $b_lines{$b_linenr}{right} - $a_right;

    #if (MIN_OVERLAP_FRAC < 0 || $ol_frac >= MIN_OVERLAP_FRAC || MAX_BORDER_DISTANCE_SUM < 0 || abs($dist_left) + abs($dist_right) <= MAX_BORDER_DISTANCE_SUM) {
    #  push(@matches, "MATCH=inexact(OF=$ol_frac;DL=$dist_left;DR=$dist_right);" . $b_lines{$b_linenr}{report});
    #  next;
    #}
  }
  $a_fields{COLUMNNAME()} = sum(@matches);
  $a_fields{COLUMNNAME()} = 0 if (! defined($a_fields{COLUMNNAME()}));

#  $a_fields{COLUMNNAME()} = '';
#  if (REPORTLEVEL < 3) {
#    foreach (@matches) {
#      if (/^MATCH=exact/) { # if we have an exact match report only this
#	$a_fields{COLUMNNAME()} = ($a_fields{COLUMNNAME()}) ? $a_fields{COLUMNNAME()} . '&' . $_ : $_;
#      }
#      if (REPORTLEVEL == 1 && /^MATCH=position/) { # if we have a position match report only this 
#	$a_fields{COLUMNNAME()} = ($a_fields{COLUMNNAME()}) ? $a_fields{COLUMNNAME()} . '&' . $_ : $_;
#      }
#    }
#  }
#  if (! $a_fields{COLUMNNAME()}) {
#    $a_fields{COLUMNNAME()} = @matches ? join('&', @matches) : '.'; # else report all matches
#  }
  #say join "\t", map { "$_:$a_fields{$_}" } @a_columns;
  say join "\t", @a_fields{@a_columns};
} # AFILE_LOOP
close A;
close $b_fh if (ref($b_fh));
