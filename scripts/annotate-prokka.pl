#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Std;

$Getopt::Std::STANDARD_HELP_VERSION = 1;
my $version = "gb.pl v. 0.1\n";
my $help    = "Use as:
  annotate.pl -p {prokka gbk} -g {gmsn gff3} -l {locus tag} -t {tRNAscan-SE txt} -c {codon table} -M {protein min length, ignored} [-m {manual corrections in gbk}] > output.gb
";

our($opt_v, $opt_h, $opt_p, $opt_g, $opt_l, $opt_t, $opt_c, $opt_m, $opt_M);

die $help if not getopts('hvHi:p:g:l:t:c:m:M:') or $opt_h;
die $version if $opt_v;

sub VERSION_MESSAGE {
	print { shift } $version;
}

sub HELP_MESSAGE {
	print { shift } $help;
}

# check input args
die "Input prokka file not specified\n"   if not $opt_p;
die "Input gmsn gff file not specified\n" if not $opt_g;
die "No locus tag specified\n"            if not $opt_l;
$opt_c = 1 if not $opt_c;
$opt_l =~ s/[^_\w]/_/g;

my $prokka = Bio::SeqIO->new(-file => $opt_p)     or die "$!\n";
my $gmsn   = Bio::Tools::GFF->new(-file => $opt_g)   or die "$!\n";
my $out    = Bio::SeqIO->newFh(-format => 'genbank') or die "$!\n";

my %features;

while (my $feature = $gmsn->next_feature()) {
	next if $feature->primary_tag ne "CDS";
	my ($name) = split / /, $feature->seq_id;
	$features{$name} = () if not $features{$name};

	$feature->add_tag_value("inference", "COORDINATES:profile:" . $feature->source_tag);
	push @{$features{$name}}, $feature;
}

if ($opt_m) {
	my $manual = Bio::SeqIO->new(-file => $opt_m) or die "$!\n";
	while (my $seq = $manual->next_seq) {
		my $name = $seq->id;
		$features{$name} = () if not $features{$name};
		for my $feature ($seq->get_SeqFeatures) {
			next if $feature->primary_tag ne "CDS";
			my $i = 0;
			my $add = 1;
			my @remove = ();
			for my $feat0 (@{$features{$name}}) {
				my $startDiff = $feature->start - $feat0->start;
				my $endDiff   = $feature->end   - $feat0->end;
				if ($startDiff <= 0 and $endDiff >= 0) {
					#my $note = sprintf "Replaces an annotation made with %s (%s..%s)", $feat0->source_tag, $feat0->start, $feat0->end;
					#$feature->add_tag_value("note", $note);
					if ($feature->has_tag("note") and $feat0->start == $feature->start and $feat0->end == $feature->end) {
						my ($note) = $feature->get_tag_values("note");
						push @remove, $i if $note eq "remove";
					}
					$features{$name}[$i] = $feature;
					$add = 0;
					last;
				}
			} continue {
				$i += 1;
			}
			push @{$features{$name}}, $feature if $add;
			for my $i (@remove) {
				splice @{$features{$name}}, $i;
			}
		}
	}
}

if (defined $opt_t) {
	open TRNA, $opt_t or die "$!\n";
	while (<TRNA>) {
		chomp;
		my ($name, $num, $begin, $end, $type, $codon, $intronBegin, $intronEnd, $score) = split /\s+/;
		next if $num !~ /^\d+$/;

		my $strand = $begin > $end ? -1 : 1;
		my $gene = sprintf "tRNA-%s(%s)", $type, $codon;
		my $product = sprintf "tRNA-%s(%s)", $type, $codon;

		my $tag = { inference => "COORDINATES:profile:tRNAscan-SE", score => $score, gene => $gene, product => $product };

		my $loc;
		if ($intronBegin > 0) {
			$loc = Bio::Location::Split->new();
			$loc->add_sub_Location(Bio::Location::Simple->new(-start => $begin, -end => $intronBegin - 1, -strand => $strand));
			$loc->add_sub_Location(Bio::Location::Simple->new(-start => $intronEnd + 1, -end => $end, -strand => $strand));
		} else {
			$loc = Bio::Location::Simple->new(-start => $begin, -end => $end, -strand => $strand);
		}

		my $feature = Bio::SeqFeature::Generic->new(-location => $loc, -strand => $strand, -primary_tag => 'tRNA', -tag => $tag);

		$features{$name} = () if not $features{$name}; 
		push @{$features{$name}}, $feature;
	}
	close TRNA;
}

my $geneNum = 0;

while (my $seq = $prokka->next_seq) {
	my $len = length $seq->seq;
	my $name = $seq->id;
	my @new_features = ();
	for my $feature ($seq->get_SeqFeatures) {
		next if $feature->primary_tag ne "CDS";
		my ($translation) = $feature->get_tag_values("translation");
		next if $translation =~ /[*]/;
		my ($locus_tag) = $feature->get_tag_values("locus_tag"); 
		my $add = 1;
		my $i = 0;
		for my $feat0 (@{$features{$name}}) {
			my $startDiff = $feature->start - $feat0->start;
			my $endDiff   = $feature->end   - $feat0->end;
			if ($startDiff <= 5 and $endDiff >= -5 or $startDiff >= -5 and $endDiff <= 5) {
				if ($feature->strand == $feat0->strand and ($feature->start == $feat0->start or $feature->end == $feat0->end)) {
					if ($feat0->end - $feat0->start < $feature->end - $feature->start) {
						$features{$name}[$i] = $feature;
					} elsif ($feature->has_tag("product")) {
						my ($product) = $feature->get_tag_values("product");
						$feat0->add_tag_value("product", $product) if $product ne "hypothetical protein";
					}
				}
				$add = 0;
			}
		} continue {
			$i += 1;
		}
		push @new_features, $feature if $add;
	}
	$seq->remove_SeqFeatures();
	$seq->desc(".");
	$seq->annotation(Bio::Annotation::Collection->new());

	for my $feature (sort { $a->start <=> $b->start } (@{$features{$name}}, @new_features)) {

		if ($feature->primary_tag eq "CDS") {
			my $cds = Bio::Seq->new(-seq => $seq->subseq($feature->start, $feature->end));
			$cds = $cds->revcom if $feature->strand < 0;
			my $translation;
			if ($feature->has_tag("translation")) {
				($translation) = $feature->get_tag_values("translation");
			} else {
				$translation = $cds->translate(-codontable_id => $opt_c)->seq;
			}
			# prokka corrects non-ATG starts to hyphen
			$translation =~ s/^-/M/;
			my $ini = substr($translation, 0, 1);
			my $stop = substr($translation, -1);
			my $strand = $feature->strand;

			my $location;
			if ($feature->start < 4 and ($strand < 0 and $stop ne "*" or $strand > 0 and $ini ne "M")) {
				$location = Bio::Location::Fuzzy->new(-start => $feature->start, -end => $feature->end, -strand => $feature->strand, -start_fuz => 'BEFORE');
			}
			if ($len - $feature->end < 4 and ($strand > 0 and $stop ne "*" or $strand < 0 and $ini ne "M")) {
				$location = Bio::Location::Fuzzy->new(-start => $feature->start, -end => $feature->end, -strand => $feature->strand, -end_fuz => 'AFTER');
			}
			if (defined $location) {
				my $newFeature = Bio::SeqFeature::Generic->new(-location => $location, -primary_tag => $feature->primary_tag);
				for my $tag (qw(note inference score product)) {
					next if not $feature->has_tag($tag);
					my ($val) = $feature->get_tag_values($tag);
					$newFeature->add_tag_value($tag, $val);
				}
				$feature = $newFeature;
			}
			$translation =~ s/[*]$//;

			$feature->remove_tag("translation") if $feature->has_tag("translation");
			$feature->add_tag_value("translation", $translation);

			$feature->add_tag_value("product", "hypothetical protein") if not $feature->has_tag("product");
			$feature->remove_tag("locus_tag")  if $feature->has_tag("locus_tag");
			$feature->remove_tag("protein_id") if $feature->has_tag("protein_id");

			$geneNum += 10;
			my $locus_tag = sprintf("%s_%04d", $opt_l, $geneNum);
			$feature->add_tag_value("locus_tag", $locus_tag);
			$feature->add_tag_value("protein_id", $locus_tag);
			my $gene = Bio::SeqFeature::Generic->new(-location => $feature->location, -primary_tag => "gene");
			$gene->add_tag_value("locus_tag", $locus_tag);
			$gene->add_tag_value("protein_id", $locus_tag);
			$seq->add_SeqFeature($gene);
		}
		$seq->add_SeqFeature($feature);
	}
	print $out $seq;
}
