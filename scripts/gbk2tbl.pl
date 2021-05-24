#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "gbk2tbl v. 0.1\n";
my $help    = "Use as:
  gbk2tbl < genbank_file.gbk > feature_table.tbl
";

our($opt_v, $opt_h);

die $help if not getopts('hv') or $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	print { shift } $version;
}

sub HELP_MESSAGE {
	print { shift } $help;
}

my $genbank = Bio::SeqIO->new(-fh => \*STDIN, -format => 'genbank') or die "$!\n";

while (my $seq = $genbank->next_seq) {
	printf ">Features %s\n", $seq->id;
	for my $feat ($seq->get_SeqFeatures) {
		my $rev  = $feat->location->strand < 0;
		my @locs = $feat->location->each_Location;
		@locs = reverse @locs if $rev;

		my $fst = 1;
		for my $loc (@locs) {
			my $start = sprintf "%s%d", pos_type($loc->start_pos_type, $rev), $loc->start;
			my $end   = sprintf "%s%d", pos_type($loc->end_pos_type,   $rev),   $loc->end;
			my $line  = sprintf "%s	%s", $rev ? $end : $start, $rev ? $start : $end;
			$line .= sprintf "	%s", $feat->primary_tag if $fst;
			$fst = 0;
			print $line . "\n";
		}
		for my $tag ($feat->get_all_tags) {
			next if $tag eq 'translation' or $tag eq 'score';
			for my $val ($feat->get_tag_values($tag)) {
				printf "			%s	%s\n", $tag, $val;
			}
		}
	}
}

sub pos_type {
	my $type = shift;
	my $rev = shift;
	return '<' if $type eq 'BEFORE' and not $rev or $type eq 'AFTER' and $rev;
	return '>' if $type eq 'BEFORE' and $rev or $type eq 'AFTER' and not $rev;
	return '';
}
