#!/usr/bin/perl
use warnings;
use strict;


# Single-sample exomes are not usually encoded with positions in which no variants are found
# However, these positions may be useful for detecting UPD signatures, in some cases, depending on genotypes of the other family members
# How do we decide when to add homozygous reference positions to our exome files?
# We assume that positions in highly covered regions (exome-targetted regions with high coverage) that are 
#  polymorphic at that position (polymorphic='not rare', 5%-95% in allele frequency in the population)
#  are absent from the exome because no variant was present in sequencing data

# Using the Agilent SureSelect v3 captured regions, and using captured regions of reasonably high depth
#  we made a list of positions for which we are confident we should have called a variant if it was in the sequencing data.
# The file containing this list of positions is included, and called 'common_variants_within_well_covered_target_regions.txt'

# If your exome design is not the Agilent SureSelect v3, your captured regions will differ the ones we used
# So, we provide this script for you to generate your own list of positions that you can add as homozygous reference to your exomes
# if your exomes are not encoded with homozygous reference positions

# With this file you can design your own file of polymorhpic-highly covered positions
# You'll need to supply a file holding minor-allele-frequencies (1000G is publically available)
# This file isn't uploaded due to size constraints
# But, it looks like this:
#CHROM	POS	rsID	REF	ALT	AF_AFR	AF_AMR	AF_ASN	AF_EUR	AF_MAX
#1	10583	rs58108140	G	A	0.040650	0.171271	0.131119	0.207124	0.207124
#1	10611	rs189107123	C	G	0.014228	0.027624	0.013986	0.021108	0.027624

# The MedianRead depth file looks like this:
# #chr	start	end	Median_100Samples
#1	14367	15164	75.44
#1	15571	16090	8.22
#1	16491	18548	98.84
#1	18699	19270	49.17
#1	20503	20823	101.72
#1	24348	25015	28.69

# Note; you'll need tabix in your bin for this to run.

my $tab_file_holding_mafs = "annots-rsIDs-AFs.2012-07-19.tab.gz";
my $tsv_file_holding_target_regions = "MedianReadDepthAcross100Samples.txt";

open (REGION_COVERAGE, $tsv_file_holding_target_regions) or die "$!: Can't open region coverage file\n";
while (<REGION_COVERAGE>) {
	chomp;
	next if ($_ =~ /^#/);

	#	1	14367	15164	75.44

	my ($chr,$start, $end, $avg_coverage) = split /\t/, $_;
	next if ($avg_coverage < 30);

	my $decile = int(($end - $start)/10);
	# avoid the outer 10%s of the region which generally have lower coverage
	my $new_start = $start + $decile;
	my $new_end = $end - $decile;

	# Get MAFs for positions in this good-coverage region.
	# Can improve this with an iterator
	my $query = `tabix $tab_file_holding_mafs $chr:$new_start-$new_end`;
	chomp $query;
	my @query = split /\n/, $query;
	my %maf = ();
	for (@query) {
		my $line_ref = [ split /\t/, $_ ];
		my ($chr, $pos, $ref, $alt, $maf) = ($line_ref->[0], $line_ref->[1], $line_ref->[3], $line_ref->[4], $line_ref->[8]);
		next if (length($ref) > 1 || length($alt) > 1);
		$maf{$chr}{$pos} = [ $maf, $ref, $alt ];
	}

	my $value = $new_start;
	while ($value <= $new_end) {
		if (exists($maf{$chr}{$value})) {
			my ($maf, $ref, $alt) = @{$maf{$chr}{$value}};
			# if the variant is common we should be able to see it if it's really there, so we'll guess it's probably homozygous reference
			if ($maf >= 0.05 && $maf < 0.95) {
				print "$chr\t$value\t$maf\t$ref\t$alt\n"	
			}
			# reduce memory footprint
			delete($maf{$chr}{$value});
		}
		$value++
	}
}
