#!/usr/bin/env perl

# for help, type "perldoc UPDio.pl" at the command line.

# Load Modules
use warnings;
use strict;
use Statistics::R;
use Data::Dumper; 
use Path::Class;
use Getopt::Long;
use Vcf;
use Iterator::Simple qw( iter );
use List::MoreUtils qw( any uniq all );
use Cwd;
use Math::Round qw(nearest_floor);

# Global variables
my $log_fh;
my $genotype_events_fh;
my $table_fh;
my $upd_fh;

my %opt = ( 
	verbose                     =>	1, 
	plots                       =>	1, 
	increase_cnv_filtering      =>	0,
	significance_level          =>	0,
	include_MI                  =>	0,
	include_X                   =>	0,
	common_cnv_file             => "sample_data/common_cnvs_1percent.tsv",
	output_path                 => cwd() . "/output_dir",
	R_scripts_dir               => "scripts"
);
	

#MAIN BLOCK
{ 
	my ($child_vcf_path, $mom_vcf_path, $dad_vcf_path, $cnv_status) = process_options();

	my $child_name = get_child_name_from_vcf_path($child_vcf_path);

	my $output_dir = open_output_filehandles($child_name);
    
	print $log_fh "Analysing UPD in proband $child_name trio\n";

	my $position_events_href = parse_vcfs_to_gather_genotype_events( $child_vcf_path, $mom_vcf_path, $dad_vcf_path );
	
	my $tabulate_event_href = tabulate_events($position_events_href);

	my $binomial_tests_result_aref = run_binomial_test_on_tabulated_events($tabulate_event_href);

	my $significant_events_aref = is_significant($binomial_tests_result_aref);

	output_upd($significant_events_aref,$position_events_href, $cnv_status);
	
	plot_upd($significant_events_aref, $child_vcf_path, $child_name, $output_dir);

	close_output_filehandles();

	print $log_fh "Done!  Finished analyzing $child_name.\n";
}

# Subroutines
sub process_options {
	my @options = qw(	child_vcf=s mom_vcf=s dad_vcf=s 
						child_cnv_data=s R_scripts_dir=s
						output_path=s common_cnv_file=s significance_level=f 
						plots increase_cnv_filtering include_MI include_X
	);
	
    GetOptions(\%opt, @options);
	
	for my $person ( qw(child_vcf mom_vcf dad_vcf) ) {
		if (! $opt{$person}) {
			die "FATAL: Please supply $person file\n";
		}
		if (! -r $opt{$person} ) {
			die "FATAL: Cannot read supplied $person file\n. Check it exists and read permissions are enabled.\n" 
		}
	}
	return ( $opt{child_vcf}, $opt{mom_vcf}, $opt{dad_vcf} );
}

sub get_child_name_from_vcf_path {
	my $child_vcf_path = shift;
	my $child_name = $child_vcf_path;
	$child_name =~ s/\.gz$//;
	$child_name =~ s/\.vcf$//;
	$child_name =~ s/^.*\///;
	
 	# my ($child_name) = $child_vcf_path =~ m{(.*/)?(.+)\.vcf(?:\.gz)?$};
	if (! $child_name ) { print Dumper $child_vcf_path, $child_name; die "FATAL: problem formatting child_name\n" };
    return $child_name;
}

sub open_output_filehandles {
	my $child_name = shift;
	my $output_dir = $opt{output_path};
	my $dir = dir($output_dir);
	if ( ! -r $dir) {
		$dir->mkpath() or die "$! can't create directory $opt{output_path}; do you have write permissions there?\n"
	}
    $log_fh 	= file("$output_dir" . "/" . "$child_name.log")->open('w') || die "Problem making log file handle\n";
    $table_fh 	= file( "$output_dir" . "/" . "$child_name.table" )->open('w') || die "Problem making table file handle\n";
    $upd_fh 	= file( "$output_dir" . "/" . "$child_name.upd" )->open('w') || die "Problem making upd file handle\n";
	if ($opt{verbose}) {
    	$genotype_events_fh = file( "$output_dir" . "/" . "$child_name.events_list" )->open('w');
	}
	return($output_dir);
}

sub close_output_filehandles {
	if ($genotype_events_fh) {
    	$genotype_events_fh->close or die "FATAL: can't close $opt{'output_dir'}: $!";
	}
}

sub parse_vcfs_to_gather_genotype_events {
	my ( $child_vcf_path, $mom_vcf_path, $dad_vcf_path) = @_;
	print $log_fh "Preparing vcfs for reading...\n";	
	my $child_vcf  = parse_vcf( $child_vcf_path );
	my $mom_vcf = parse_vcf( $mom_vcf_path );
	my $dad_vcf = parse_vcf( $dad_vcf_path );
	my $cnv_filter_href = load_cnv_data();
	my $child_gt  = $child_vcf->next;
	my $mom_gt = $mom_vcf->next;
	my $dad_gt = $dad_vcf->next;
	my %position_events;	
	while ( $child_gt && $mom_gt && $dad_gt ) {
	    my $child_mom = cmp_pos( $child_gt, $mom_gt );
	    if ( $child_mom < 0 ) {
			$child_gt = $child_vcf->next;
			next;
	    }
	    elsif ( $child_mom > 0 ) {
			$mom_gt = $mom_vcf->next;
			next;
	    }
	    else { # child_gt == mom_gt
			my $child_dad = cmp_pos( $child_gt, $dad_gt );
			if ( $child_dad < 0 ) {
			    $child_gt = $child_vcf->next;
			    next;
			}
			elsif ( $child_dad > 0 ) {
			    $dad_gt = $dad_vcf->next;
			    next;
			}
			else { # child_gt == dad_gt
				my $upd_status = compare_genotypes( $child_gt, $mom_gt, $dad_gt, $cnv_filter_href );
				if ($upd_status) {
					my ($chr, $pos) = ($child_gt->{chr}, $child_gt->{pos});
					$position_events{$chr}{$pos} = $upd_status;
				}
			    $child_gt = $child_vcf->next;
			    $mom_gt = $mom_vcf->next;
			    $dad_gt = $dad_vcf->next;
			}
	    }
	}
	return(\%position_events);
}

sub parse_vcf {
    my $path = shift;
	print $log_fh "Reading genotypes for $path...\n";
    my $vcf = Vcf->new( file => $path );
    $vcf->parse_header();
    return iter sub {
		my $x_href = $vcf->next_data_hash() or return;
		if (keys %{ $x_href->{gtypes} } != 1) {
			die "This version does not handle multi-sample VCF\n" 
		};
		if (! chr2num($x_href->{CHROM})) { print Dumper $x_href; die "FATAL: I do not undersatnd chromosome identity in this line.\n"};
		if (! $x_href->{ID} ) { print Dumper $x_href; die "FATAL: I do not understand ID in this line.\n" };
		if ($x_href->{ID} eq "common") { 
			return {
				chr    => chr2num($x_href->{CHROM} ),
			    pos    => $x_href->{POS},
			    ref    => $x_href->{REF},
			    alt    => $x_href->{ALT},
			    gt     => "00",
			    filter => $x_href->{FILTER},
				id     => "common"
			}
		}
		else {	
			my $gt = ( values %{ $x_href->{gtypes} } )[0]->{GT};
			if (! $gt) { print "no genotype found!\n", Dumper $path, $x_href;}
			$gt =~ s{[/|]}{};
			return { 
			    chr    => chr2num( $x_href->{CHROM} ),
			    pos    => $x_href->{POS},
			    ref    => $x_href->{REF},
			    alt    => $x_href->{ALT},
			    gt     => $gt,
			    filter => $x_href->{FILTER},
			    id     => $x_href->{ID}
			};
		}
    }
}

sub chr2num {
    my $chr = shift;
    $chr =~ s/^[Cc]hr//i;

    if ( $chr =~ m/\d+/ ) {
		return $chr;
    }
    elsif ( $chr =~ /^X$/) {
			return 23
	}
    elsif ( $chr =~ /^[YM]$/) {
			return
	}
	else {
  		die "Check chromosome name: $chr \n";
	}
}

sub cmp_pos {
    my ( $x, $y ) = @_;
	# If the chromosomes are the same, then compare the positions.  
	# <=> comparison operator; returns 0 if x == y, 1 if x > y, -1 if x < y
	# ( ex: 14 <=> 14 = 0.  The left side is 0, so the || then compares the right side )
	if (! $x || ! $y) { print Dumper $x, $y; die "FATAL: I do not understand comparison\n" };
	if ($x->{chr} - $y->{chr} > 1) { die "$x->{chr} and $y->{chr} !!  Files not sorted correctly\n" };
	if ($y->{chr} - $x->{chr} > 1) { die "$x->{chr} and $y->{chr} !!  Files not sorted correctly\n" };

    return $x->{chr} <=> $y->{chr} || $x->{pos} <=> $y->{pos};
}


sub is_wanted_genotype {
    my ( $child, $mom, $dad, $cnv_filter_href ) = @_;

	# if any uses List::MoreUtils module.  it's like the 'map' function
    return if all { $_->{id} eq 'common' } $child, $mom, $dad;

    my %zygosity;
	return if position_lies_within_known_cnv($child->{chr},$child->{pos},$cnv_filter_href);
	
	for my $gt ( $child, $mom, $dad ) {
		return unless ( ${$gt->{filter}}[0] =~ m/(^[.]$)|(PASS)/ );
		return unless scalar @{$gt->{alt}} == 1;
		return unless $gt->{alt}[0] =~ m/^[ACGT]$/;
		return unless $gt->{gt} =~ m{^[01][01]$};
		return unless $gt->{ref} =~ m/^[ACGT]$/;
		$zygosity{$gt->{alt}[0]}++;	
    }
    return if keys %zygosity > 1;
    return 1;
}

sub position_lies_within_known_cnv {
	my ($chr,$pos,$cnvs_href) = @_;
	my $nearest_hashkey = nearest_floor(1E6,$pos);
	for my $cnv_region (@ { $cnvs_href->{$chr}{$nearest_hashkey} }) {
		my @interval = split /[,]/, $cnv_region;
		my $start = $interval[0];
		my $end = $interval[1];
		if ($pos >= $start && $pos <= $end) {
			return 1
		}
	}
	return 0
}

sub compare_genotypes {
    my ( $child_gt, $mom_gt, $dad_gt, $cnv_href ) = @_;
	my %position_events;

	return unless is_wanted_genotype( $child_gt, $mom_gt, $dad_gt, $cnv_href );

	my $upd_status = calculate_upd_status( $child_gt->{gt}, $mom_gt->{gt}, $dad_gt->{gt} );

    my ($chr, $pos, $ref, $alt) = ( $child_gt->{chr}, $child_gt->{pos}, $child_gt->{ref}, ${$child_gt->{alt}}[0]);

	my ($childs_geno, $moms_geno, $dads_geno) = ( $child_gt->{gt}, $mom_gt->{gt}, $dad_gt->{gt} );
	if ($opt{verbose}) {
		my $translated_status = translate_nomenclature($upd_status);
		print $genotype_events_fh "$chr\t$pos\t$ref\t$alt\t$childs_geno\t$moms_geno\t$dads_geno\t$translated_status\n";
	}
    return if ( $upd_status eq 'uninformative' );

	return($upd_status);
}

sub tabulate_events {
	my $position_events_href = shift;
	print $log_fh "Tabulating events...\n";
	my %tabulate_event;
	for my $chr (keys %{$position_events_href}) {
		for my $pos (keys %{$position_events_href->{$chr}}) {	
			my $event = $position_events_href->{$chr}{$pos};
			my @events = split /[,]/,$event;
			for (@events) {$tabulate_event{$chr}{$_}++};
			$tabulate_event{$chr}{informative}++;
		}
		for my $informative_event_type qw(iUPI_P hUPI_P BPI iUPI_M hUPI_M MI_S MI_D informative) {
			if (!$tabulate_event{$chr}{$informative_event_type}) {
				$tabulate_event{$chr}{$informative_event_type} = 0
			}
			$tabulate_event{all_chrs}{$informative_event_type} += $tabulate_event{$chr}{$informative_event_type};
		}
		if ($tabulate_event{all_chrs}{informative} == 0) { die "it is VERY unlikely no informative events on any chr; i bet there is a problem with how the genotypes are recorded, or there are no genotypes\n"};
	}
	my @events = qw(iUPI_P hUPI_P iUPI_M hUPI_M BPI informative);
	if ($opt{include_MI}) {
		@events = qw(iUPI_P hUPI_P iUPI_M hUPI_M MI_S MI_D BPI informative);
	}

#	print $table_fh "chr\t", join ("\t", @events),"\n";
	print $table_fh "chr\t";
	for (@events) {
		print $table_fh translate_nomenclature($_), "\t";
	}
	print $table_fh "\n";
	my @chrs;
	if ( $opt{include_X} ) { 
		print $log_fh "including X chromosome\n";
		@chrs = 1..23 
	}	
	else { #Default
		print $log_fh "analysing autosomes\n";
		@chrs = 1..22
	};
	for my $chr (@chrs,"all_chrs") {
		print $table_fh "$chr\t";
		for my $informative_event_type (@events) {
			if (!$tabulate_event{$chr}{$informative_event_type}) {$tabulate_event{$chr}{$informative_event_type}=0};
			print $table_fh translate_nomenclature($tabulate_event{$chr}{$informative_event_type}),"\t"
		}
		print $table_fh "\n";
	}
	return(\%tabulate_event);
}

sub translate_nomenclature {
	my $before = shift;
	my %trans = (
					'iUPI_P' => "UI_P",
					'hUPI_P' => "AI_P",
					'iUPI_M' => "UI_M",
					'hUPI_M' => "UI_P",
					'hUPI_M,iUPI_M' => "AI_M",
					'hUPI_P,iUPI_P' => "AI_P",
					BPI => "B"
	);
	if (exists $trans{$before}) {
		return $trans{$before}
	}
	else {
		return $before
	}
}

sub calculate_concordance {
    my ( $x, $y ) = @_;
	if (! $x || ! $y ) { print Dumper $x, $y; die "FATAL: I do not understand comparison being made in calculate_concordance\n"};
    my $x1 = substr($x,0,1);
    my $x2 = substr($x,1,1);
    my $y1 = substr($y,0,1);
    my $y2 = substr($y,1,1);

    if ( ( $x1 eq $y1 && $x2 eq $y2 ) || ( $x1 eq $y2 ) && ( $x2 eq $y1 ) ) {
		return 1;
    }
    elsif ( ( $x1 eq $y1 ) || ( $x1 eq $y2 ) || ( $x2 eq $y1 ) || ( $x2 eq $y2 ) ) {
		return 0.5;
    }
    else {
		return 0;
    }
}

sub calculate_zygosity {
    my ( $gt ) = @_;
    # Return 1 for homozygous, 2 for heterozygous genotype
	# uniq is List::MoreUtils to take unique items from a list
    scalar uniq split '', $gt
}

sub is_lookup_match {
    my ( $params, $status_href ) = @_;
    for my $col ( grep { $_ ne 'value' } keys %{$status_href} ) {
        return unless defined $params->{$col} and $params->{$col} == $status_href->{$col};
    }
    return 1;
}

sub calculate_upd_status {
    my ( $child_gt, $mom_gt, $dad_gt ) = @_;
    my %params = (
		child_zygos => calculate_zygosity( $child_gt ),
		mom_zygos   => calculate_zygosity( $mom_gt ),
		dad_zygos   => calculate_zygosity( $dad_gt ),
		mom_dad     => calculate_concordance( $mom_gt, $dad_gt ),
		child_mom   => calculate_concordance( $child_gt, $mom_gt ),
		child_dad   => calculate_concordance( $child_gt, $dad_gt )
    );

	### UPD Testing & Cataloging Informative SNPS
	# inspired by Figure 1 of TrioUPD paper published in Human Mutation 2007 with the following title:
	# "Visualization of uniparental inheritance, Mendelian inconsistencies, deletions,
	# and parent of origin effects in single nucleotide polymorphism trio data with SNPtrio."
	# Informative sites:  Mendelian Incompatibility singly or doubly (MI-S, MI-D) (not included by default)
	# Informative sites:  Uniparental isodisomy ("Absolute Isodiomsy") or heterodisomy ("Undistinguished Disomy") maternal -M or paternal -P
	# Informative sites:  Bi-parental inheritance (BPI)
	my @event_lookup_table = (
	    { child_zygos => 2, mom_zygos => 1, dad_zygos => 1,	mom_dad => 0,	child_mom => 0.5,	child_dad => 0.5,	value => 'BPI'},
	    { child_zygos => 1, mom_zygos => 2, dad_zygos => 1, mom_dad => 0.5,	child_mom => 0.5,	child_dad => 0,		value => 'iUPI_M'},
	    { child_zygos => 1, mom_zygos => 1, dad_zygos => 1, mom_dad => 0,	child_mom => 0,		child_dad => 1,		value => 'hUPI_P,iUPI_P'},
	    { child_zygos => 1, mom_zygos => 1, dad_zygos => 1, mom_dad => 0,	child_mom => 1,		child_dad => 0,		value => 'hUPI_M,iUPI_M'},
	    { child_zygos => 1, mom_zygos => 1, dad_zygos => 2, mom_dad => 0.5,	child_mom => 0,		child_dad => 0.5,	value => 'iUPI_P'}
	);
	if ($opt{include_MI}) {
		push @event_lookup_table, (
	    { child_zygos => 2, mom_zygos => 1, dad_zygos => 1, mom_dad => 1,	child_mom => 0.5,	child_dad => 0.5,	value  => 'MI_S'},
	    { child_zygos => 1, mom_zygos => 1,	dad_zygos => 1,	mom_dad => 1,	child_mom => 0,		child_dad => 0,		value => 'MI_D'},
	) };

    for my $row ( @event_lookup_table ) {
		if ( is_lookup_match( \%params, $row ) ) {
		    return $row->{value};
		}
    }
    return "uninformative";
}

sub load_cnv_data {
	# Load CNVs
	my $common_cnvs_href = load_common_cnvs();
	my $cnv_positions_in_child_href = (load_sample_specific_cnvs($common_cnvs_href));
	if (! $cnv_positions_in_child_href) { die "Problem loading cnv data\n" };
	my $wobble = 10000;
	my %wobble_cnv_positions_in_child;
	for my $chr (keys %{$cnv_positions_in_child_href}) {
		while (my $range = shift @{$cnv_positions_in_child_href->{$chr}}) {
			my $interval = [ split /[,]/, $range ];
			my ($start,$end) = ($interval->[0],$interval->[1]);
			my $new_start = $start - $wobble;
			my $new_end = $end + $wobble;
			my $ceiling = nearest_floor(1E6, $new_end);
			my $floor = nearest_floor(1E6, $new_start);
			my $current_position = $ceiling;
			while ($current_position >= $floor) {
				my $nearest_mil = nearest_floor(1E6, $current_position);
				push @{$wobble_cnv_positions_in_child{$chr}{$nearest_mil}}, "$new_start,$new_end";
				$current_position -= 1E6;
			}
		}
	}
	return(\%wobble_cnv_positions_in_child);
}

sub load_common_cnvs {
	my ( %common_cnvs, $common_cnv_fh);
	my $common_cnv_file = $opt{common_cnv_file};
	open ($common_cnv_fh, "<", $common_cnv_file) or die "$!: Can't open common cnv file common_cnv_file: did you specify a proper location with the appropriate option?\n";
	while (<$common_cnv_fh>) {
		chomp;
		next if ($_ =~ /^#/);
		my ($cnv_chr,$cnv_start,$cnv_stop) = split /\t/, $_;
		push @{$common_cnvs{$cnv_chr}}, "$cnv_start,$cnv_stop";
	}
	return(\%common_cnvs);
}

sub load_sample_specific_cnvs {
	my $common_cnvs_href = shift;
	my $sample_cnvs_href = $common_cnvs_href;
	my $sample_cnv_path = $opt{child_cnv_data} || 0;
	my %sample_cnvs;
	my $cnv_fh;
	if ( -r $sample_cnv_path) {
		open ($cnv_fh, "<", $sample_cnv_path) or die "YIkes, can't read common cnv file $sample_cnv_path!\n";
		print $log_fh "Sample-specific CNV data found for proband\n";	
		while (my $cnv_line = <$cnv_fh>){
			next if ($cnv_line =~ /^#/);
			my @fields = split /\t/, $cnv_line;
			my ($cnv_chr, $cnv_start, $cnv_stop, $cnv_copies) = @fields;
			$cnv_chr =~ s/[Cc]hr//;
			push @{$sample_cnvs_href->{$cnv_chr}}, "$cnv_start,$cnv_stop";
		}
	}
	else {
		print $log_fh "WARNING: no proband cnv data available; deletions can mimic upd events\n";
	}
	return($sample_cnvs_href);
}

sub	run_binomial_test_on_tabulated_events {
	my $tabulated_events_href = shift;
	print $log_fh "Running binomial test...\n";
	my %tabulated_events = %{$tabulated_events_href};
	my @binom_test_output;
	my @chrs_to_run_tests_on = grep  { $_ ne "all_chrs" } (keys %tabulated_events);
	my @events_to_test;
	if ($opt{include_MI}) {
		@events_to_test = qw(iUPI_P hUPI_P iUPI_M hUPI_M MI_S MI_D);
	}
	else {
		@events_to_test = qw(iUPI_P hUPI_P iUPI_M hUPI_M);
	}
	my $R = Statistics::R->new();
	for my $chr (@chrs_to_run_tests_on) {
		next if ($tabulated_events{$chr}{informative} == 0);
		for my $event_type (@events_to_test) {
			next if ($event_type eq "B");
			next if ($tabulated_events{$chr}{$event_type} == 0);
			my $event_avg = $tabulated_events{all_chrs}{$event_type} / $tabulated_events{all_chrs}{informative};
			my $success = $tabulated_events{$chr}{$event_type};
			my $failure = $tabulated_events{$chr}{informative} - $tabulated_events{$chr}{$event_type};
			$R->set('success', $success); $R->set('failure',$failure); $R->set('prob_of_success', $event_avg);
			$R->send('binom.test(x=c(success,failure),p=prob_of_success,alternative="g")$p.value');
			(my $binom_test_pval = $R->read) =~ s/.*( )(.*)/$2/;
			push (@binom_test_output , [ $binom_test_pval, $event_type , $chr ]);
		}
	}
	if (scalar @binom_test_output == 0) {
		print $log_fh "None\n"
	}
	return(\@binom_test_output);
}

sub is_significant {
	my ($binomial_tests_results_aref) = @_;
	my @binomial_tests_results = @{$binomial_tests_results_aref};
	my $p_val;
	if ($opt{significance_level}) {
		$p_val = $opt{significance_level};
		print $log_fh "User defined p_value cut off of $p_val will be used\n"
	}
	else {
		my $num_chr = 22;
		$num_chr = 23 if ($opt{include_X});
		my $number_of_tests = $num_chr * 4;
		my $alpha = 0.05;
		$p_val = $alpha/$number_of_tests;
		print $log_fh "Testing significance at p-value: $alpha/$number_of_tests = ", $alpha/$number_of_tests, ".\n";
	}
	my @significant_results;
	for (@binomial_tests_results) {
		my ($sig,$event, $chr) = @{$_};
		if ($sig <= $p_val) {
			push @significant_results, \@{$_}
		}
	}
	if (scalar @significant_results == 0) {
		print $upd_fh "None\n";
	}
	return (\@significant_results);
}

sub output_upd {
	my ($significant_results_aref,$position_events_href, $cnv_status) = @_;
	print $log_fh "Printing results...\n";
	my @significant_results = @{$significant_results_aref};
	my %position_events = %{$position_events_href};
	my %sig_results;
	for (@significant_results) {
		my ($sig, $event, $chr) = @{$_};
		my @positions_for_this_significant_result;
		for my $pos (sort { $a <=> $b } keys %{$position_events{$chr}}) {
			if ($position_events{$chr}{$pos} =~ /$event/) {
				push @positions_for_this_significant_result, $pos;
			}
		}
		push @{$sig_results{$sig}},  [ $event , $chr , \@positions_for_this_significant_result ];
	}
	
	for my $sig (sort {$a<=>$b} keys %sig_results) {
		for my $check_for_multiple_results_with_same_pval (@{$sig_results{$sig}}) {
			my ($event, $chr, $positions) = @{$check_for_multiple_results_with_same_pval};
			if (! $opt{child_cnv_data} && $opt{increase_cnv_filtering} && cluster_test($positions) ) {
				print $upd_fh "small_cluster_of_positions_consider_CNV\t$sig\t$event\t$chr\t", join ("\t", @{$positions}), "\n";
			}
			else {
				print $upd_fh "$sig\t$event\t$chr\t", join ("\t", @{$positions}), "\n";
			}
		}
	}
	return 1;
}

sub cluster_test {
	my $positions_aref = shift;
	my $R = Statistics::R->new();
	$R->set('positions',$positions_aref);
	$R->send('mad(positions)');
	(my $pos_mad = $R->read) =~ s/(.*[ ]+)(.*)/$2/;
	if ($pos_mad < 5000000) {
		return(1);
	}
}

sub plot_upd {
	my ($significant_events_aref, $child_vcf_path, $child_name, $output_dir) = @_;
	if (! $opt{plots}) {
		print $log_fh "No plot option selected. Will not plot\n";	
		return
	}
	if (! $opt{verbose}) {
		print $log_fh "will plot zygosity\n";

		if (scalar @{$significant_events_aref} == 0) { 
			print $log_fh "No events to plot\n"
		}
		plot_zygosity($child_vcf_path, $output_dir);
		return
	}	
	plot_zygosity_and_position_events($child_name, $child_vcf_path, $output_dir);
	return 1;
}
	
sub plot_zygosity {
	my ($child_vcf_path, $output_dir) = @_;
	if (! $output_dir || ! $child_vcf_path ) { print Dumper $child_vcf_path, $output_dir; die "Problem with preparing data for plotting\n"};
	print $log_fh "Plotting zygosity\n";
	`R-2.14.1 --slave -f $opt{R_scripts_dir}/vcf_to_zygosity_plot.R --args $child_vcf_path $output_dir 2> /dev/null`;
	return 1;
}

sub plot_zygosity_and_position_events {
	my ($child_name, $child_vcf_path, $output_dir) = @_;
	if (! $output_dir || ! $child_vcf_path ) { print Dumper $child_vcf_path, $output_dir; die "Problem with preparing data for plotting\n"};
	print $log_fh "Plotting positions and their events...\n";
	$output_dir =~ s/\/$//;
	my $events_list_path = "$output_dir" . "/" . "$child_name.events_list"; 
	my $plot_sex = "FALSE";
	$plot_sex = "TRUE" if $opt{include_X};
	`R --slave -f $opt{R_scripts_dir}/plot_zygosity_and_events.R --args $events_list_path $child_vcf_path $output_dir $plot_sex 2> /dev/null`;	
	return 1;
}	 

## END SUB

__END__

=over 25

=head1 NAME

 UPDio.pl

=head1 VERSION

 Version 0.9 (Pre-Publication Beta)

=head1 DESCRIPTION 

 UPDio is designed to detect UPD events in probands from vcf trio data.

=head1 SYNOPSIS

 
 perl UPDio.pl --child_vcf $child_file --mom_vcf $mom_file --dad_vcf $dad_file --child_cnv_data $child_cnv_data_file


=head2 Required Arguments


=item C<child_vcf>

the child vcf file (or vcf.gz file)

=item C<mom_vcf> 

the mom vcf file (or vcf.gz file)

=item C<dad_vcf> 

the dad vcf file (or vcf.gz file)

=head2 Recommended Arguments

=item C<child_cnv_data>

CNV data for the proband; useful to avoid mis-calling CNV events as UPD events

=head2 Options

=item C<output_path>

for specifying your own output path; use absolute path

=item C<include_X>

for processing the X chromomsome; only makes sense for female probands

=item C<significance_level>

for specifying your own significance value, accepts floating point number, such as "0.0001"

=item C<common_cnv_file>

for specifying your own tab-separated list of cnv (or other) regions to avoid


=item C<increase_cnv_filtering>

adds annotation meant to flag regions that are likely to be smaller cnv regions (unlikely to be UPD)

=item C<noplots>

disables plotting feature

=item C<noverbose>

quiets logging output

=item C<include_MI>

for searching for mendelian inconsistent sites; occassionally useful for finding loss of transmitted allele

=back

=head1 EXAMPLE

 perl UPDio.pl --child_vcf Family101_child.vcf.gz --mom_vcf Family101_mom.vcf.gz --dad_vcf Family101_dad.vcf.gz --child_cnv_data Family101_child_cnvs.txt

=head1 AUTHORS

 Dan King  < dk6 at sanger.ac.uk >
 Coding oversight by Ray Miller
 Project oversight by Matt Hurles

=head1 INSTITUTION

 Wellcome Trust Sanger Institute
 Wellcome Trust Genome Campus,
 Hinxton, Cambridge, CB10 1SA, UK
 
=head1 LICENSE

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  For more information, please go to <http://www.gnu.org/licenses/>.

=cut
