#!/usr/bin/env perl

# DK6 19_4_2013
# Synapsis:
# Reads through pre-processed VCF and if a polymorphic_site is absent from VCF, it gets added as a homREF site.
# We only add hom-ref sites that are in generally highly covered regions (Agilent SureSelect Version3)
# If you have your own exome capture design you can supply your own list of highly covered polymorphic sites
#  using a script we provide: get_common_variants_in_target_region_file.pl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Iterator::Simple qw(iter);
use Vcf;

sub process_options {
	my %opt;
	GetOptions(\%opt, qw(polymorphic_sites=s no_homREF_vcf=s input));
	if (! $opt{polymorphic_sites})	{ die "please include --polymorphic_sites \n"};
	if (! $opt{no_homREF_vcf} ) 	{ die "please include --no_homREF_vcf\n"};
	return \%opt
}

# MAIN
{
	my $opt_href = process_options();

	my $common_site_iter = generate_polysites_iter($opt_href->{polymorphic_sites});
	my ($sample_id, $vcf_iter) = generate_vcf_iter($opt_href->{no_homREF_vcf});	
	my $cs_inst = $common_site_iter->next;
	my $vcf_inst = $vcf_iter->next;

	print "##fileformat=VCFv4.1\n";
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_id\n";

	my ($chr, $pos, $ref, $alt, $id, $filter, $gt);
	while ($cs_inst && $vcf_inst) {
		last if (! $cs_inst || ! $vcf_inst);
		my $cs_vcf_cmp = cmp_pos($cs_inst, $vcf_inst);
		if ($cs_vcf_cmp < 0) {	# CS position is before VCF position
#			print "vcf has passed common site so presumably the common site was hom-ref and is not annotated in this VCF file\n";
			$chr = $cs_inst->{chr};
			$pos = $cs_inst->{pos};
			$ref = $cs_inst->{ref};
			$alt = $cs_inst->{alt};
			$id  = "common";
			$filter = "PASS";
			$gt = "0/0";
			$cs_inst = $common_site_iter->next;
		}
		elsif ( $cs_vcf_cmp == 0) {		# some position both files; use the vcf genotype at this position
#			print "USE VCF genotype\n";
			$chr = $vcf_inst->{chr};
			$pos = $vcf_inst->{pos};
			$ref = $vcf_inst->{ref};
			$alt = $vcf_inst->{alt};
			$id  = $vcf_inst->{id};
			$gt = $vcf_inst->{gt};
			$filter = join (",", @{$vcf_inst->{filter}});
			$cs_inst = $common_site_iter->next;
			$vcf_inst = $vcf_iter->next;
		}
		else {						# VCF position needs to catch up, it is before the CS position; use its genotype
#			print "vcf needs to catch up to cs file\n";
			$chr = $vcf_inst->{chr};
			$pos = $vcf_inst->{pos};
			$ref = $vcf_inst->{ref};
			$alt = $vcf_inst->{alt};
			$id  = $vcf_inst->{id};
			$gt = $vcf_inst->{gt};
			$filter = join (",", @{$vcf_inst->{filter}});
			$vcf_inst=$vcf_iter->next;
					
		}
#		print Dumper "vcf", $vcf_inst, "common", $cs_inst;
	#	print Dumper $chr, $pos, $id, $ref, $alt, $filter, 'FORMAT', $gt;
#		print "I will print\t$chr\t$pos\t$id\t$ref\t$alt\t$filter\tFORMAT\t$gt\n";
		print "$chr\t$pos\t$id\t$ref\t$alt\t.\t$filter\t.\tGT\t$gt\n";

	}
	# Purge last entry
	if ($vcf_inst) {
		$chr = $vcf_inst->{chr};
		$pos = $vcf_inst->{pos};
		$ref = $vcf_inst->{ref};
		$alt = $vcf_inst->{alt};
		$id  = $vcf_inst->{id};
		$gt = $vcf_inst->{gt};
		$filter = join (",", @{$vcf_inst->{filter}});
		print "$chr\t$pos\t$id\t$ref\t$alt\tQUAL\t$filter\tINFO\tGT\t$gt\n";
	}
}

sub generate_polysites_iter {
	my $polysite_file = shift;
	open (my $polysite_fh, "<", "$polysite_file") or die "$! can't open polymorphic_in_TRs_file\n";
	# start iterator subroutine
	# an iterator is just like reading through a file line by line
	# http://www.perl.com/pub/2005/06/16/iterators.html
	return iter sub {
		while (my $line = <$polysite_fh>) {
			chomp $line;
			my ($chr,$pos,$maf,$ref,$alt) = split /\t/, $line;
			return {
				chr => chr2num($chr),
				pos => $pos,
				maf => $maf,
				ref => $ref,
				alt => $alt
			}
		}
	}
}	

sub generate_vcf_iter {
	my $vcf_file = shift;
	my $vcf = Vcf->new( file => $vcf_file );
	$vcf->parse_header();
	my $sample_id = $vcf->get_samples();
	return $sample_id, iter sub {
		my $x_href = $vcf->next_data_hash() or return;
		if (keys %{ $x_href->{gtypes} } != 1) {
        	die "This version does not handle multi-sample VCF\n" 
		};
		die if (! chr2num($x_href->{CHROM})); # Fix this if this is broken
        my $gt = ( values %{ $x_href->{gtypes} } )[0]->{GT};
        if (! $gt) { die }
#        $gt =~ s{[/|]}{};
        return { 
            chr    => chr2num( $x_href->{CHROM} ),
            pos    => $x_href->{POS},
            ref    => $x_href->{REF},
            alt    => join (",", @{$x_href->{ALT}}),
            gt     => $gt,
            filter => $x_href->{FILTER},
            id     => $x_href->{ID}
        }
    }
}

sub chr2num {
	my $in = shift;
    (my $chr = $in) =~ s/^[Cc]hr//i;
	return $chr if $chr =~ m/^\d(\d)?$/;
	return 23   if $chr =~ /^X$/;
	return 24   if $chr =~ /^Y$/;
	return 25   if $chr =~ /^XY$/;
	return 26   if $chr =~ /^M$/;
	die "Check chromosome name: $chr \n";
}


sub cmp_pos {
    my ( $x, $y ) = @_;
	# If the chromosomes are the same, then compare the positions.  
	# <=> comparison operator; returns 0 if x == y, 1 if x > y, -1 if x < y
	# ( ex: 14 <=> 14 = 0.  The left side is 0, so the || then compares the right side )
    return $x->{chr} <=> $y->{chr} || $x->{pos} <=> $y->{pos};
}

