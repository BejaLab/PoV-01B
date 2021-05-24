#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "biotags.pl v. 1.2\n";
my $help    = "Use as:
  biotags.pl -i {input} [-p {primary}] -t {tags} [-H] [-d {delimiter}] > stdout

 -i - input file to parse
 -p - comma-separated list of primary tags
 -a - comma-separated list of accessions
 -H - whether to include header in the output
 -d - field delimiter (tab is the default)
 -T - comma-separated list of sequence tags (call -T? to show the list)
 -t - comma-separated list of feature tags. Any named tag can be used
      (call -t? to show the list of predefined tags)

from https://github.com/har-wradim/miscngs/
";

our($opt_v, $opt_h, $opt_i, $opt_T, $opt_t, $opt_p, $opt_H, $opt_d, $opt_a);

die $help    if not getopts('hvHi:a:T:t:p:d:') or $opt_h;
die $version if $opt_v;
if (defined $opt_T and $opt_T eq '?') {
	print "accession (=accession_number)
accession_number
alphabet (dna, rna or protein)
annotation (defaults to comment)
    comment
    date_changed
    keyword
    reference
    # etc.
authority
can_call_new (useless)
definition (=desc)
desc
description
display_id
feature_count
get_num_of_annotations
id
is_circular
keywords
length
locus (=display_id)
namespace
object_id
organism (=species)
primary_id
revcom
seq
source (=species)
species (defaults to binomial)
    binomial
    classification
    common_name
    division
    genus
    ncbi_taxid
    organelle
    species
    sub_species
    taxon (not implemented)
    tree (not implemented)
    variant
translate
version
";
	exit;
}

	if (defined $opt_t and $opt_t eq '?') {
		print "display_name
end
entire_seq
get_all_tags
get_tagset_values
get_tag_values
gff_string
location (formatted)
phase
primary_id
primary_tag
seq
seq_id
source_tag
spliced_seq
start
_static_gff_formatter (not implemented)
strand
";
	exit;
}

sub VERSION_MESSAGE {
	print { shift } $version;
}

sub HELP_MESSAGE {
	print { shift } $help;
}

# check input args
die "No input file specified\n" if not    $opt_i;
print STDERR "No tags specified\n" if not $opt_t and not $opt_T;

$opt_d = "\t" if not defined $opt_d or $opt_d eq '\t';

my %primaries;
if (defined $opt_p) { $primaries{$_}  = 1 for split /,/, $opt_p }

my %accessions;
if (defined $opt_a) { $accessions{$_} = 1 for split /,/, $opt_a }
my $accn = -1;
$accn = scalar keys %accessions if defined $opt_a;

my @stags = ();
my @ftags = ();

@stags = split /,/, $opt_T if $opt_T;
@ftags = split /,/, $opt_t if $opt_t;

# print the header if requested
printf "%s\n", join $opt_d, (@stags, @ftags) if defined $opt_H;

# open the input and iterate over individual seqs/features
my $in = Bio::SeqIO->new(-file => $opt_i) or die "$!\n";
while (my $seq = $in->next_seq) {
	next if $opt_a and not $accessions{$seq->id};
	if ($opt_t) {
		for my $feature ($seq->get_SeqFeatures) {
			next if $opt_p and not $primaries{$feature->primary_tag};
			my @vals = ();
			push @vals, get_seq_values($seq, $_)     for @stags;
			push @vals, get_tag_values($feature, $_) for @ftags;
			printf "%s\n", join $opt_d, @vals if @vals;
		}
	}
	else {
		my @vals = ();
		push @vals, get_seq_values($seq, $_) for @stags;
		printf "%s\n", join $opt_d, @vals if @vals;
	}
	last if --$accn == 0;
}

# get value for a given feature/tag pair
sub get_tag_values {
	my $feature = shift;
	my $tag = shift;

	# first check explicit tags
	if ($feature->has_tag($tag)) {
		return join ';', $feature->get_tag_values($tag);
	}

	return '' if $tag =~ /^(new|can)$/;

	# then object members
	if ($feature->can($tag)) {
		my @vals;
		eval { @vals = $feature->$tag };
		return '' if not defined $vals[0];
		my $val = $vals[0];
		my $ref = ref($val);
		return $val->seq if $ref eq 'Bio::PrimarySeq';
		return feature_location($feature) if $ref =~ /Bio::Location/;
		return join ";", @vals;
	}
	return '';
}

# get value for a given seq/tag pair
sub get_seq_values {
	my $seq = shift;
	my $tag = shift;

	return '' if $tag =~ /^(new|can)$/;

	my %aliases = (
		locus      => 'display_id',
		definition => 'desc',
		accession  => 'accession_number',
		organism   => 'species',
		source     => 'species'
	);

	$tag = $aliases{$tag} if defined $aliases{$tag};
	if ($seq->can($tag)) {
		my @vals;
		eval { @vals = $seq->$tag };
		return '' if not defined $vals[0];
		my $val  = $vals[0];
		return '' if not defined $val;
		my $ref  = ref($val);
		return $val->seq                     if $ref eq 'Bio::Seq::RichSeq';
		return $val->binomial                if $ref eq 'Bio::Species';
		return $val->get_all_annotation_keys if $ref eq 'Bio::Annotation::Collection';
		return join ";", @vals;
	}
	my @ac = $seq->get_Annotations($tag);
	if (@ac) {
		my @vals;
		push @vals, $_->display_text for (@ac);
		return join ';', @vals;
	}
	my $species = $seq->species;
	if ($species->can($tag)) {
		my @vals;
		eval { @vals = $species->$tag };
		return '' if not defined $vals[0];
		return join ';', @vals if defined $vals[0];
	}
	return '';
}

# extract join'ed segment coordinates
sub feature_location {
	my $feature = shift;
	my $loc = $feature->location;
	return location_range($loc) if not $loc->isa('Bio::Location::SplitLocationI');
	my @f = ();
	push @f, location_range($_) for ($loc->sub_Location);
	return join ",", @f;
}

sub location_range {
	my $loc = shift;
	return sprintf "%s%s..%s%s", pos_type($loc->start_pos_type), $loc->start, pos_type($loc->end_pos_type), $loc->end;
}

sub pos_type {
	my $pos = shift;
	my %vals = ( BEFORE => '<', AFTER => '>', EXACT => '' ); # WITHIN => '', BETWEEN => ''
	return $vals{$pos} if defined $vals{$pos};
	return '';
}
