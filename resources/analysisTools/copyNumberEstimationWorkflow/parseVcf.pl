#!/usr/bin/perl
use strict;
use warnings;

my %filehandle; my $filename; my $filenametmp;

for (my $i=1; $i < 25; $i++){
	my $id = $i;
	$id =~ s/23/X/;
	$id =~ s/24/Y/;
	$filename = "$ARGV[0]/$ARGV[1]$id.$ARGV[2]";
	$filenametmp = "$ARGV[0]/$ARGV[1]$id.$ARGV[2].tmp";
	open( $filehandle{$id}, ">", $filenametmp) or die "Could not open $filenametmp!\n" ;
	print {$filehandle{$id}} "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcontrol\n";
}

while (my $line = <STDIN> ){
	my @fields = split("\t", $line);
	next if $fields[8] eq ".\n";
	next if $fields[8] eq "0/0\n";
	$fields[9] = $fields[8];
	$fields[8] = "GT";
	if ($filehandle{$fields[0]}){
		print {$filehandle{$fields[0]}} join("\t", @fields);
	}
}

for (my $i=1; $i < 25; $i++){
	my $id = $i;
	$id =~ s/23/X/;
	$id =~ s/24/Y/;
	$filename = "$ARGV[0]/$ARGV[1]$id.$ARGV[2]";
	$filenametmp = "$ARGV[0]/$ARGV[1]$id.$ARGV[2].tmp";
	close( $filehandle{$id});
	my $move=system("mv $filenametmp $filename 2>&1");
	die "Could not move $filenametmp!\n" if $move != 0;
}
