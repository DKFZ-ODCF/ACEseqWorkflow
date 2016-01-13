#!/usr/bin/env perl

use strict;
use warnings;
use v5.10;

my %opts;
my @standardColnames;

BEGIN {
  use Getopt::Long;
  use Pod::Usage;
  
  if (! @ARGV) {
    pod2usage({-verbose => 1, -exitval => 2, -output => \*STDERR});
  }
  
  %opts = ('tabix_bin' => 'tabix',
	   #'bgzip_bin' => 'bgzip',
	   'aFileType' => 'vcf',
	   'bFileType' => 'vcf',
	   'padding' => 0,
	   'maxBorderDistanceSum' => -1, # negative values in one of these disable this filter
	   'minOverlapFraction' => -1,   # CAVE: it is one filter --- both conditions are connected by OR
	   'reportLevel' => 2, # 1: exact or position or all; 2: exact or all; 3: all,
	   'aPosOffset' => 0,
	   'aEndOffset' => 0,
	   'maxNrOfMatches' => 0,
	   'aColNameLineStart' => '#CHR',
	   'noMatchString' => '.',
	   'emptyBFields' => 0,
	  );



  GetOptions( 'a=s' => \$opts{afile},
	      'b=s' => \$opts{bfile},
	      'columnName=s' => \$opts{columnname},
	      'padding:i' => \$opts{padding},
	      'tabix_bin:s' => \$opts{tabix_bin},
	      #'bgzip_bin:s' => \$opts{bgzip_bin},
	      'bReportColumn:i' => \$opts{bReportColumn},
	      'bAdditionalColumns:s' => \$opts{bcolumns},
	      'minOverlapFraction:f' => \$opts{minOverlapFraction},
	      'maxBorderDistanceSum:i' => \$opts{maxBorderDistanceSum},
	      'bFileType:s' => \$opts{bFileType},
	      'aFileType:s' => \$opts{aFileType},
	      'aChromColumn:s' => \$opts{aChromColumn},
	      'aPosColumn:s' => \$opts{aPosColumn},
	      'aEndColumn:s' => \$opts{aEndColumn},
	      'aPosOffset:i' => \$opts{aPosOffset},
	      'aEndOffset:i' => \$opts{aEndOffset},
	      'aColNameLineStart:s' => \$opts{aColNameLineStart},
	      'reportLevel:i' => \$opts{reportLevel},
	      'reportMatchType!' => \$opts{reportMatchType},
	      'reportBFeatCoord!' => \$opts{reportBFeatCoord},
	      'reportOnlyMatches!' => \$opts{reportOnlyMatches},
	      'maxNrOfMatches:i' => \$opts{maxNrOfMatches},
	      'noMatchString:s' => \$opts{noMatchString},
	      'breakPointMode!' => \$opts{breakPointMode},
	      'chromXtr:s' => \$opts{chromXtr},
	      'chromYtr:s' => \$opts{chromYtr},
	      'emptyBFields' => \$opts{emptyBFields},
	      'debug!' => \$opts{debug},
	      'help' => \$opts{help},
	    ) || pod2usage({-verbose => 1, -message => "$!", -exitval => 2, -output => \*STDERR});

  if ($opts{help}) {
    pod2usage({-verbose => 2, -exitval => 1});
  }

  ### Process user input and check validity
  if (! $opts{afile} || ! $opts{bfile} || ! $opts{columnname}) {
    pod2usage({-verbose => 1, -message => "Required parameter missing (a, b and columnname are required)", -exitval => 2, -output => \*STDERR});
  }

  $opts{bFileType} = lc($opts{bFileType});

  my %valid_bFileTypes = ( 'vcf' => 1,
			   'bed' => 1,
			   'gff3' => 1,
			   'bed_refvar' => 1,
			 );
  $valid_bFileTypes{$opts{bFileType}} || die "Unknown b-file type $opts{bFileType}. Only b-files of type ". join(" ", sort keys %valid_bFileTypes). " are supported";

  $opts{aFileType} = lc($opts{aFileType});
  my %valid_aFileTypes = ( 'vcf' => 1,
			   'custom' => 1,
			 );
  $valid_aFileTypes{$opts{aFileType}} || die "Only a-files of type ". join(" ", sort keys %valid_aFileTypes). " are supported";

  if (! defined($opts{bReportColumn})) {
    $opts{bReportColumn} = ($opts{bFileType} eq 'vcf') ? 7 : ($opts{bFileType} eq 'gff3') ? 8 : 3;
  }

  ### Now make constants
  require constant;
  for (keys %opts) {
    constant->import(uc($_), $opts{$_});
  }
  if ($opts{aFileType} eq 'vcf') {
    constant->import(CHROMCOL => 'CHROM');
    constant->import(POSCOL => 'POS');
  } elsif ($opts{aFileType} eq 'custom') {
    constant->import(CHROMCOL => $opts{aChromColumn});
    constant->import(POSCOL => $opts{aPosColumn});
    constant->import(ENDCOL => $opts{aEndColumn});
  }

}
if (BFILETYPE eq 'vcf') {
  @standardColnames = qw(CHROM POS ID REF ALT QUAL FILTER INFO);
} elsif (BFILETYPE eq 'gff3') {
  @standardColnames = qw(seqid source type start end score strand phase attributes);
} elsif (BFILETYPE eq 'bed_refvar') {
  @standardColnames = (qw(chrom chromStart chromEnd ), 'ref>var');
} else {
  @standardColnames = qw(chrom chromStart chromEnd name);
}


DEBUG && say "MINOVERLAPFRACTION: ", MINOVERLAPFRACTION(), "MAXBORDERDISTANCESUM: ", MAXBORDERDISTANCESUM();


#my $newcol = $opts{columnname};
#my $pad = $opts{padding};
 
open(A, "$opts{afile}") || die "Could not open a-file $opts{afile}\n";

if (! -e $opts{bfile})
{
	die "The b-file $opts{bfile} does not exist, exiting!\n";
}
if (! -e "$opts{bfile}.tbi")
{
	die "There is no .tbi index for b-file $opts{bfile}, exiting!\n";
}

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
if (BFILETYPE ne 'gff3') { # gff3 files have no column names in header
  my $b_header_cmd = TABIX_BIN() . " -h $opts{bfile} $b_chr_prefix" . '1' . $b_chr_suffix . ":0-0 |";
  open(HEAD, $b_header_cmd);
  my @b_header = <HEAD>;
  
  #### if I have a multi-line header print out all lines but the last
  #for (my $i=0; $i < @b_header-1; $i++) {
  #  print $b_header[$i];
  #}
  
  ### the last line should contain the column names
  if (@b_header) {
    my $b_header = $b_header[-1];
    chomp($b_header);
    $b_header =~ s/^\#//;
    @b_colnames = split(/\t/, $b_header);
  }
}
if (! @b_colnames) {
  @b_colnames =  @standardColnames;
}


DEBUG && say "b_colnames: @b_colnames";


my $newcolidx;
my %a_fields;
my @b_fields;
my @a_columns;
my %a_columns;
my @b_columns = split(',',$opts{bcolumns}) if (defined($opts{bcolumns}));
my $chr='';
my $chr_raw = '';
my $mh;
my $rs;

my $header;
while ($header = <A>) {
  last if ($header =~ /^$opts{aColNameLineStart}/i); # that is the line with the column names
  print $header; # print out every preceeding line
  die "Invalid a-file header" if ($header =~ /^[^\#]/);
}

chomp $header;
my $poundkey_rm = 1 if ($header =~ s/^\#//);
$header =~ s/(?:(?<=\t)\s+)||(?:\s+(?=\t|$))//g; ### remove leading and trailing whitespace around column names
@a_columns = split(/\t/, $header);
$newcolidx = 0;
foreach (@a_columns) {
  if ($_ eq COLUMNNAME()) {
    last;
  }
  $newcolidx++;
}
if ($newcolidx == @a_columns) {
  say $poundkey_rm ? '#' : '', $header, "\t", COLUMNNAME();
  push @a_columns, COLUMNNAME();
} else {
  say $poundkey_rm ? '#' : '', $header;
}

%a_columns = map { $_ => 1 } @a_columns;

if (! exists $a_columns{CHROMCOL()}) {
  die "CHROMCOL not found in a-file ";
}

if (! exists $a_columns{POSCOL()}) {
  die "POSCOL not found in a-file ";
}

if ($opts{aFileType} eq  'custom' && ! exists $a_columns{ENDCOL()}) {
  die "ENDCOL not found in a-file ";
}


my @matches;
my $a_line;
my $b_line;
my @b_lines;
my $b_fh;
my $next_b_line;

my ($a_left, $a_right);
my $chr_changed;
# my $b_linectr = 0; # need this as hash key
# my $b_linenr;
my ($ol_left, $ol_right, $ol_length, $ol_frac);
my ($dist_left, $dist_right);
my ($chrXa, $chrXb, $chrYa, $chrYb);
($chrXa, $chrXb) = split(':', CHROMXTR()) if (CHROMXTR());
($chrYa, $chrYb) = split(':', CHROMYTR()) if (CHROMYTR());
my $alt;
my %a_alts;

AFILE_LOOP: while ($a_line=<A>) {
  @matches = ();
  chomp($a_line);
  @a_fields{@a_columns} = split(/\t/, $a_line);
  # DEBUG && say join '; ', map {"$_ => $a_fields{$_}"} @a_columns;
  if ($a_fields{CHROMCOL()} ne $chr_raw) {
    $chr_changed = 1;
    $chr_raw = $a_fields{CHROMCOL()};
    $chr = $chr_raw =~ s/[^\dXY]*([\dXY]+).*/$1/r;
    $chr = $chrXb if (CHROMXTR() && $chr eq $chrXa);
    $chr = $chrYb if (CHROMYTR() && $chr eq $chrYa);
  } else {
    $chr_changed = 0;
  }
  #DEBUG && say "CHR: $chr CHR_RAW: $chr_raw ($a_fields{CHR})";

  $a_left = $a_fields{POSCOL()};
  if (AFILETYPE eq 'vcf' ) {
  #  if (defined($a_fields{INFO}) && $a_fields{INFO} =~ /END=(\d+)/) {
  #    $a_right = $1 - 1;
  #  } else {
      $a_right = $a_left + length($a_fields{REF}) - 1;
  #  }
  } else {
    $a_right = $a_fields{ENDCOL()};
  }

  if (APOSOFFSET()) {
    $a_left += APOSOFFSET();
  }
  if (AENDOFFSET()) {
    $a_right += AENDOFFSET();
  }
  

  if ($chr_changed) {
    @b_lines = ();
    $next_b_line = {};
    close $b_fh if (ref($b_fh));
    open $b_fh, TABIX_BIN . " $opts{bfile} ${b_chr_prefix}${chr}${b_chr_suffix} |" or die "opening b file $opts{bfile} with tabix failed";
    warn "Tabix returned no b-features for chromosome ${b_chr_prefix}${chr}${b_chr_suffix}" if (! ref($b_fh));
    # $b_linectr = 0;
  }
  if ((! defined($next_b_line->{left}) || $next_b_line->{left}-PADDING() <= $a_right && ref($b_fh))) {
    # read new b_lines until we have one where the left coordinate is higher than a_right + pad
    while ($b_line=<$b_fh>) {
      if (defined($next_b_line->{left})) {
	push(@b_lines,$next_b_line);
	$next_b_line = {};
      }
      chomp($b_line);
      DEBUG && say $b_line;
      # $b_linectr++;
      @b_fields = split(/\t/, $b_line);
      # ($b_fields[0] eq $b_chr_prefix.$chr.$b_chr_suffix) || die "Chromosomes between tumor and germline file do not match: $b_fields[0] ne ${b_chr_prefix}${chr}${b_chr_suffix}";
      if (BFILETYPE eq 'vcf') { ### b-file is vcf file
	$next_b_line->{left} =  $b_fields[1];
        $next_b_line->{ref} = $b_fields[3];
        $next_b_line->{alt} = $b_fields[4];
        $next_b_line->{report} = 'POS='.$next_b_line->{left}.';' if (REPORTBFEATCOORD);
	$next_b_line->{right} = $next_b_line->{left} + length($next_b_line->{ref}) - 1;
        if ($b_fields[7] !~ /END=(\d+)/) {
          $next_b_line->{report} .= 'END='.$next_b_line->{right}.';' if (REPORTBFEATCOORD);
        }
      } elsif (BFILETYPE eq 'gff3') { ### b-file is gff3 file
	$next_b_line->{left} =  $b_fields[3];
	$next_b_line->{right} =  $b_fields[4];  # here we read END as the last included base
      } elsif (BFILETYPE eq 'bed_refvar') {
	$next_b_line->{left} =  $b_fields[1] + 1;
	$next_b_line->{right} = $b_fields[2]; # here we read END as the last included base
	($next_b_line->{ref},$next_b_line->{alt}) = split(/>/,$b_fields[3]);
      } else { ### b-file is bed file
	$next_b_line->{left} =  $b_fields[1] + 1;
	$next_b_line->{right} = $b_fields[2]; # here we read END as the last included base
      }

      if ($next_b_line->{right}+PADDING() < $a_left) { # b-feature is left of a-feature; skip...
	$next_b_line = {};
	next;
      }
      if (REPORTBFEATCOORD && (BFILETYPE eq 'gff3' || BFILETYPE eq 'bed' || BFILETYPE eq 'bed_refvar')) {
	$next_b_line->{report} = 'POS='.$next_b_line->{left}.';END='.$next_b_line->{right}.';';
      }
      if (EMPTYBFIELDS()) {
	for (@b_columns) {
	  $b_fields[$_] = '.' if (! defined($b_fields[$_]));
	}
      }
      if (BREPORTCOLUMN() >= 0 ) {
	$next_b_line->{report} .= join ';',  $b_fields[BREPORTCOLUMN()], (map { $b_colnames[$_] . '=' .$b_fields[$_] } @b_columns);
      } else {
	$next_b_line->{report} .= join ';',  (map { $b_colnames[$_] . '=' .$b_fields[$_] } @b_columns);
      }
      last if ($next_b_line->{left}-PADDING() > $a_right);
    }
    if (defined($next_b_line->{left}) && $next_b_line->{left}-PADDING() <= $a_right) { # the while loop terminated because eof at b_fh
      push(@b_lines, $next_b_line);
      $next_b_line = {};
    }
  }
  # Now compare a_line to all b_lines in b_lines hash
 B_LINE_LOOP: foreach $b_line (@b_lines) {
    if ( $b_line->{right}+PADDING() < $a_left ) {
      $b_line = undef; # delete($b_lines{$b_linenr}); #####
      next;
    }
    if ( $b_line->{left}-PADDING() > $a_right ) {
      last;
    }
    if ((BFILETYPE eq 'vcf' || BFILETYPE eq 'bed_refvar') && AFILETYPE() eq 'vcf') {
      %a_alts = ();
      map { $a_alts{uc($_)} = 1 } split(',',$a_fields{ALT});
      foreach $alt (split(',', $b_line->{alt})) { # to check for exact match we must go through all alt alleles in b files
	if ( $b_line->{left} == $a_left && $b_line->{right} == $a_right && uc($b_line->{ref}) eq uc($a_fields{REF}) && $a_alts{uc($alt)} ) {
	  $mh = {mt => 'exact',
		 mr => $b_line->{report},
		 ofsd => 3};
	  push(@matches, $mh);
	  next B_LINE_LOOP ;
	}
      }
    }
    next B_LINE_LOOP if (REPORTLEVEL == 4);
    if ( $b_line->{left} == $a_left && $b_line->{right} == $a_right ) {
      $mh = {mt => 'position',
	     mr => $b_line->{report},
	     ofsd => 2};
      push(@matches, $mh);
      next;
    }
    if ($b_line->{left} <= $a_right  && $b_line->{right} >= $a_left) {
      if (BREAKPOINTMODE()) {
	next if ($b_line->{left}-PADDING() > $a_left && $b_line->{right}+PADDING() < $a_right);
      }
      $ol_left = ($b_line->{left} > $a_left) ? $b_line->{left} : $a_left;
      $ol_right = ($b_line->{right} < $a_right) ? $b_line->{right} : $a_right;
      $ol_length = $ol_right - $ol_left + 1;
      $ol_frac = $ol_length / (($b_line->{right}-$b_line->{left} > $a_right-$a_left) ? $b_line->{right}-$b_line->{left} + 1 : $a_right-$a_left + 1);
      $ol_frac = sprintf("%.2f", $ol_frac);
    } else {
      $ol_frac = 0;
    }
    $dist_left = $b_line->{left} - $a_left;
    $dist_right = $b_line->{right} - $a_right;

    if (MINOVERLAPFRACTION() < 0 || $ol_frac >= MINOVERLAPFRACTION() || MAXBORDERDISTANCESUM() < 0 || abs($dist_left) + abs($dist_right) <= MAXBORDERDISTANCESUM()) {
      $mh = {mt => "inexact(OF=$ol_frac;DL=$dist_left;DR=$dist_right)",
	     mr => $b_line->{report},
	     ofsd => $ol_frac};
      push(@matches, $mh);
      next;
    }
  }

  @b_lines = grep { defined } @b_lines; # this should be improved!


  $a_fields{COLUMNNAME()} = '';
  if (REPORTLEVEL < 3) {
    foreach (@matches) {
      if ($_->{mt} eq 'exact') { # if we have an exact match report only this
	$rs = REPORTMATCHTYPE() ? 'MATCH=' . $_->{mt} . ';' . $_->{mr} : $_->{mr};
	$a_fields{COLUMNNAME()} = ($a_fields{COLUMNNAME()}) ? $a_fields{COLUMNNAME()} . '&' . $rs : $rs;
	next;
      }
      if (REPORTLEVEL == 1 && $_->{mt} eq 'position') { # if we have a position match report only this (and exact matches if we have)
	$rs = REPORTMATCHTYPE() ? 'MATCH=' . $_->{mt} . ';' . $_->{mr} : $_->{mr};
	$a_fields{COLUMNNAME()} = ($a_fields{COLUMNNAME()}) ? $a_fields{COLUMNNAME()} . '&' . $rs : $rs;
      }
    }
  }

  if (MAXNROFMATCHES() && @matches) {
    @matches = (sort {$b->{ofsd} <=> $a->{ofsd}} @matches)[0 .. ((MAXNROFMATCHES()-1 < $#matches) ? MAXNROFMATCHES()-1 : $#matches)];
  }

  if (! $a_fields{COLUMNNAME()}) { # else report all matches
    if (REPORTMATCHTYPE()) {
      $a_fields{COLUMNNAME()} = join('&', (map {'MATCH=' . $_->{mt} .';' . $_->{mr}} @matches));
    } else {
      $a_fields{COLUMNNAME()} =  join('&', (map {$_->{mr}} @matches));
    }
  }
  if (! @matches) {
    next if (REPORTONLYMATCHES());
    $a_fields{COLUMNNAME()} = NOMATCHSTRING();
  }
  say join "\t", @a_fields{@a_columns};
} # AFILE_LOOP
close A;
close $b_fh if (ref($b_fh));


__END__

=head1 NAME

annotate_vcf.pl

=head1 SYNOPSIS

annotate_vcf.pl -a afile -b bfile --columnName=newcolname [options]

=head1 OPTIONS

=over 8

=item B<-a>

a-file (must be coordinate sorted; must be either vcf format (at least the 1st eight columns; additional columns are okay) or custom format where CHROM, POS and END columns are set by the respective options; type - for stdin (required)

=item B<-b>

b-file (must be vcf, gff3 or bed format; must be sorted and indexed with tabix (and thus compressed with bgzip)) (required)

=item B<--columnName>

Name of column to which overlaps with b-file should be reported. If not yet present it is added to the a-file; if already present in input a-file, the values are replaced  (required)

=item B<--bFileType>

'bed' or 'vcf' or 'bed_refvar'[defaul: vcf]

=item B<--aFileType>

'vcf' or 'custom'. For custom, aChromColumn, aPosColumn and aEndColumn must be set[defaul: vcf]

=over 4

=item B<--aChromColumn>

Name of chromosome column in a-file (starting # do not belong to the column name, e.g. in a vcf file the column name is CHROM and not #CHROM. All other starting characters are not removed and thus belong to the column name) (required for aFileType=custom)

=item B<--aPosColumn>

Name of position (start) column in a-file (required for aFileType=custom)

=item B<--aEndColumn>

Name of feature end column in a-file (required for aFileType=custom)

=back

=item B<--aColNameLineStart>

First characters of last header (i.e. line with column names); no line before the last header line must start with this identifier [default: #CHR]

=item B<--aPosOffset>

Offset of POS read from file relative to 1-based including coordinates (e.g. 1 when a-file is bed-like, but 0 is a-file is gff3-like) [default: 0]

=item B<--aEndOffset>

Offset of END read from file relative to 1-based including coordinates (e.g. 0 when a-file is gff3-like, and 0 if a-file is bed-like) [default: 0]

=item B<--chromXtr>

Translation of a-file chromosome X identifier to b-file; e.g. X:23 if a-file has X or chrX and b-file has 23 or chr23 [default: '']

=item B<--chromYtr>

see --chromXtr

=item B<--padding>

the number of bases added to both sides of each b-feature when looking for overlaps (useful for finding non-exact overlaps) [default: 0]

=item B<--tabix_bin>

the tabix binary [default: tabix]

=item B<--bReportColumn>

zero-based index of column from b-file to be reported. Only the value is reported. Set to negative value to disable [default: 7 for vcf / 3 for bed]

=item B<--bAdditionalColumns>

zero-based indices of columns from b-file which should be reported if b-feature matches; they are reported as COLUMNNAME=VALUE [default: '']

=item B<--minOverlapFraction>

minimum fraction of overlap (overlap_length / max(length(a_feature), length(b_feature))) required for a match to be reported. Set to negative value to disable this filter. THIS FILTER IS CONNECTED WITH LOGICAL OR TO maxBorderDistanceSum [default: -1]

=item B<--maxBorderDistanceSum>

maximum sum of distances between left borders of a- and b-feature and right borders of a- and b-feature.  Set to negative value to disable this filter. THIS FILTER IS CONNECTED WITH LOGICAL OR TO maxBorderDistanceSum [default: -1]

=item B<--reportLevel>

Which matches should be reported? (1: Report only exact or position if present and otherwise all; 2: Report only exact if present and otherwise all; 3: report all, 4: report only exact matches) [default: 2]

=item B<--reportMatchType>

Start report for each match with MATCH=exact or similar [default: noreportMatchType]

=item B<--reportBFeatCoord>

Report POS and END coordinates of each matching b-feature [default: noreportBFeatCoord]

=item B<--reportOnlyMatches>

Report only a-file lines which matched a b-file feature [default: noreportOnlyMatches]

=item B<--maxNrOfMatches>

Maxmimum number of matches to be reported for one a-feature. When activated, matches are sorted as follows: exact > position > inexact with decreasing overlap fraction. Set to 0 to report all matches. [default: 0]

=item B<--noMatchString>

String to print into report column if a-feature does not overlap a b-feature [default: .]

=item B<--breakPointMode>

Report only b-features which overlap a POS or END [default: nobreakpointmode]

=item B<--emptyBFields>

B-file contains lines with empty fields in bAdditionalColumns [default: noemptyBFields]

=item B<--debug>

Prints out also b-file lines (for debugging purpose only) [default: nodebug]

=back

=head1 DESCRIPTION

B<annotate_vcf.pl> overlaps features from a vcf file (or vcf-like file where the first eight columns confirm to vcf standard, or a custon tab-separated file; see below) (a-file) with a vcf-, bed- or gff3-file (b-file) allowing also nearby matches by "padding" the b-file features. Comparing REF and ALT sequences for "exact matching" is currently only supported when both a- and b-file are vcf.

Matches are reported to one column of the a-file; in case of multiple matches, these are joined by '&', so the line count in the a-file is not changed.

The a-file must be either vcf-like, or it must be a custom tab-separated file when CHROM, POS, and END column are set by the respective options. The a-file must have a header, and the last header line must contain the column names (and match /^\#chr/i or the string specified in --aColNameLineStart).

The b-file must be compressed with bgzip and indexed with tabix.

Both files must be coordinate sorted; the order of chromosomes does not matter and need's not to be consistent between both files. Also the chromosome identifier (chr1, 1, chr1.fa, ...) is not required to be consistent as long as the chromosome numbers / letters are the same ('X' and 'Y' in a-file and '23' and '24' in b-file requires to set --chrXtr and --chrYtr).

Reported positions are 1-based. POS gives the first bp within the feature, END the last bp of the feature. This is according to vcf specifications. In former versions END was the first base not included in the feature anymore, but this was changed.




=cut

