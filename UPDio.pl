#!/usr/bin/env perl
# for help, type "perldoc UPDio.pl" at the command line.

# Version 1.0
# Added multisampling (trio) calling
# Fixed table output header

# find bugs or have suggestions? email me dk6 at sanger dot ac dot uk

# Load Modules
use warnings;
use strict;
use Statistics::R;
use Data::Dumper;
use Path::Class;
use Getopt::Long;
use Vcf;
use Iterator::Simple qw( iter );
use List::Util qw(sum);
use List::MoreUtils qw( any uniq all );
use Cwd;
use Math::Round qw(nearest_floor);
use File::Basename;
use File::Spec;
use FindBin;

# Global variables
my $fh_href;

my %opt = (
    increase_cnv_filtering   => 0,
    significance_level       => 0,
    include_MI               => 0,
    include_X                => 0,
    common_cnv_file          => get_updio_dir() . "/sample_data/common_dels_1percent.tsv",
    output_path              => cwd() . "/output_dir",
    R_scripts_dir            => get_updio_dir() . "/scripts",
    path_to_R                => "R",
    high_qual                => 0,
    testing                  => 0
);


#MAIN BLOCK
{
	my $paths = process_options();

    my $child_name = get_child_name_from_vcf_path($paths);

    $fh_href = open_output_filehandles($child_name);

    print {$fh_href->{log}} "Analysing UPD in proband $child_name trio\n";

    my $events_href = gather_events( $paths );

    my $tabulated_events_href = tabulate_events($events_href);

    my $significant_results_aref = run_binomial_test_on_tabulated_events_to_get_significant_results($tabulated_events_href);

    exit_routine() if ( ! $significant_results_aref);

    print_upd_results($significant_results_aref);

    plot_upd_results($child_name,$paths);

    exit_routine();
}

# Subroutines
sub get_updio_dir {
    my $script_link = readlink __FILE__; 
    if (! File::Spec->file_name_is_absolute($script_link) ) {
      return dirname(File::Spec->rel2abs($script_link, $FindBin::Bin ));
    }
    return dirname($script_link);
}

sub process_options {
    my @options = qw(	multisample_vcf=s childID=s momID=s dadID=s
						child_vcf=s mom_vcf=s dad_vcf=s
                        child_cnv_data=s R_scripts_dir=s path_to_R=s
                        output_path=s common_cnv_file=s significance_level=f
                        testing highqual include_MI include_X no_plot no_event_list);

    GetOptions(\%opt, @options);

	# See if a multisample_vcf has been loaded
	if ($opt{multisample_vcf}) {
        if (! -r $opt{multisample_vcf} ) {
            die "FATAL: Cannot read supplied multisample_vcf file: $opt{multisample_vcf}\n. Check it exists and read permissions are enabled.\n"
        }
		if (! $opt{childID} || ! $opt{momID} || ! $opt{dadID}) {
			die "Must supply childID momID and dadID, must match VCF header sample ids\n"
		}
    	return { multisample_path => $opt{multisample_vcf} };
	}

	# If no multisample file is given, then all 3 family member vcfs must be provided
    for my $person ( qw(child_vcf mom_vcf dad_vcf) ) {
        if (! -r $opt{$person} ) {
            die "FATAL: Cannot read supplied $person file\n. Check it exists and read permissions are enabled.\n"
        }
    }
    return { child_path => $opt{child_vcf}, 
			 mom_path => $opt{mom_vcf}, 
			 dad_path => $opt{dad_vcf} }
}

sub get_child_name_from_vcf_path {
    my $paths = shift;
	# Getting childID from multisample_vcf is easy. 
	return $opt{childID} if exists $opt{childID};

	# Otherwise, calculate it from the vcf path
    my $child_name = $paths->{child_path};
    $child_name =~ s/\.gz$//;
    $child_name =~ s/\.vcf$//;
    $child_name =~ s/^.*\///;
    # my ($child_namNe) = $child_vcf_path =~ m{(.*/)?(.+)\.vcf(?:\.gz)?$};
    if (! $child_name ) {
        print Dumper $paths, $child_name; die "FATAL: problem formatting child_name\n";
    };
    return $child_name;
}

sub open_output_filehandles {
    my $child_name = shift; 
    my $output_dir = $opt{output_path};
    $output_dir =~ s/\/$//;
    my $dir = dir($output_dir);
    if ( ! -r $dir) {
        $dir->mkpath() or die "$! can't create directory $opt{output_path}; do you have write permissions there?\n"
    }
	my ($log, $genotype_events_fh, $genotype_events_path, $table, $upd_call, $fig_fh);
    $log 				= file( "$output_dir" . "/" . "$child_name.log")->open('w') || die "Problem making log file handle $output_dir, $child_name\n";
    $table				= file( "$output_dir" . "/" . "$child_name.table" )->open('w') || die "Problem making table file handle\n";
    $upd_call 			= file( "$output_dir" . "/" . "$child_name.upd" )->open('w') || die "Problem making upd file handle\n";
    $genotype_events_path =     "$output_dir" . "/" . "$child_name.events_list";
    $genotype_events_fh = file( $genotype_events_path )->open('w') || die "Problem making events_list file hande\n";
	
	return { 	log => $log,
				table => $table,
				upd_results => $upd_call,
				figure => $fig_fh,
				events => $genotype_events_fh,
                events_path => $genotype_events_path,                    
	 			output_dir => $output_dir };
}

sub gather_events {
	my $paths = shift;
	my $events_href;
	# If multisample_vcf is the input, then just read the file top to bottom
	if (exists $paths->{multisample_path}) {
		$events_href = parse_multisample_vcf_to_gather_genotype_events($paths);
	}
	else {	# otherwise, we'll stagger-read the 3 vcfs and gather the genotypes
		$events_href = parse_vcfs_to_gather_genotype_events($paths);	
	}
	return $events_href;

}

sub parse_multisample_vcf_to_gather_genotype_events {
    my $paths = shift;
	my $multisample_path = $paths->{multisample_path};
    my $multi_iter  = parse_vcf($multisample_path);
    my $multi_inst = $multi_iter->next;
	my %events;
	my $cnv_filter_href = load_cnv_data();
    while ($multi_inst) {
 	   last if ($multi_inst->{chr} >= 24 || $multi_inst->{chr} >= 24 || $multi_inst->{chr} >= 24);
 	   last if ( ( ! $opt{include_X} ) && ( $multi_inst->{chr} >= 23 || $multi_inst->{chr} >= 23 || $multi_inst->{chr} >= 23 ) );
		my $is_wanted = is_wanted( { multisample => $multi_inst, cnv_filter_href => $cnv_filter_href } );
		if (! $is_wanted) {
			$multi_inst = $multi_iter->next;
			next;
		}
		my ($chr, $pos, $event_status) = ($is_wanted->{chr}, $is_wanted->{pos}, $is_wanted->{event_status});
		$events{$chr}{$event_status}++;
		$multi_inst = $multi_iter->next;
    }
	return \%events;
}

sub parse_vcfs_to_gather_genotype_events {
	my $paths = shift;
    my ( $child_vcf_path, $mom_vcf_path, $dad_vcf_path) = ( $paths->{child_path}, $paths->{mom_path}, $paths->{dad_path} );
    print {$fh_href->{log}} "Preparing vcfs for reading...\n";
    print {$fh_href->{log}} "Reading child vcf...\n";
    my $child_vcf  = parse_vcf( $child_vcf_path );
    print {$fh_href->{log}} "Reading mom vcf...\n";
    my $mom_vcf = parse_vcf( $mom_vcf_path );
    print {$fh_href->{log}} "Reading dad vcf...\n";
    my $dad_vcf = parse_vcf( $dad_vcf_path );
    my $child_gt  = $child_vcf->next;
    my $mom_gt = $mom_vcf->next;
    my $dad_gt = $dad_vcf->next;
    my %events;
	my $cnv_filter_href = load_cnv_data();
	while ( $child_gt && $mom_gt && $dad_gt ) {
 	   last if ($child_gt->{chr} >= 24 || $mom_gt->{chr} >= 24 || $dad_gt->{chr} >= 24);
 	   last if ( ( ! $opt{include_X} ) && ( $child_gt->{chr} >= 23 || $mom_gt->{chr} >= 23 || $dad_gt->{chr} >= 23 ) );
 	   my $child_mom = cmp_pos( $child_gt, $mom_gt );
 	   if ( $child_mom < 0 ) {
 	       $child_gt = $child_vcf->next;
 	       next;
 	   } elsif ( $child_mom > 0 ) {
 	       $mom_gt = $mom_vcf->next;
 	       next;
 	   } else {                # child_gt == mom_gt
 	       my $child_dad = cmp_pos( $child_gt, $dad_gt );
 	       if ( $child_dad < 0 ) {
 	           $child_gt = $child_vcf->next;
 	           next;
 	       } elsif ( $child_dad > 0 ) {
 	           $dad_gt = $dad_vcf->next;
 	           next;
 	       } else {            # child_gt == dad_gt
 	   		my $is_wanted = is_wanted({ child_gt => $child_gt, mom_gt => $mom_gt, dad_gt => $dad_gt, cnv_filter_href => $cnv_filter_href });
 	   		
 	   		if ( ! $is_wanted) {
 	   			$child_gt = $child_vcf->next;
 	   			$mom_gt = $mom_vcf->next;
 	   			$dad_gt = $dad_vcf->next;
 	   			next;
 	   		}					
			my ($chr, $pos, $event_status) = ($is_wanted->{chr}, $is_wanted->{pos}, $is_wanted->{event_status});
			$events{$chr}{$event_status}++;

			$child_gt = $child_vcf->next;
			$mom_gt = $mom_vcf->next;
			$dad_gt = $dad_vcf->next;
 	       }
 	   }
    }
    return \%events;
}

sub is_wanted {
	# Check if genotype for child mom and dad is desirable, if so, check if the position is desirable, and if event is desirable
	my $in_hrefs = shift;
	my $cnv_filter_href = $in_hrefs->{cnv_filter_href};
	my %gt;
	my %same_alt;
	if (exists $in_hrefs->{multisample}) {	#This is a multisample VCF
		my $multi = $in_hrefs->{multisample};
		return unless is_wanted_filter( ${$multi->{filter}}[0] );
		return unless is_wanted_alt( \@{$multi->{alt}} );
		return unless is_wanted_ref($multi->{ref});
		return unless is_wanted_gt($multi->{child_gt});
		return unless is_wanted_gt($multi->{mom_gt});
		return unless is_wanted_gt($multi->{dad_gt});
        return if position_lies_within_known_cnv($multi->{chr}, $multi->{pos}, $in_hrefs->{cnv_filter_href});
		$gt { child_gt } = $multi->{child_gt};
		$gt { mom_gt } = $multi->{mom_gt};
		$gt { dad_gt } = $multi->{dad_gt}; 
		$gt { chr }    = $multi->{chr};
		$gt { pos }    = $multi->{pos};
		$gt { ref }    = $multi->{ref};
		$gt { alt }    = $multi->{alt}->[0];
	}
	else { # We are seeing the 3 vcfs
		my ($child, $mom, $dad) = ($in_hrefs->{child_gt}, $in_hrefs->{mom_gt}, $in_hrefs->{dad_gt});
		my %samples = ( child => $child, mom => $mom, dad => $dad );
		my @samples = keys %samples;
		for my $sample ( @samples ) {
			return unless ( is_wanted_filter(${$samples{$sample}->{filter}}[0]) );
			return unless ( is_wanted_alt(\@{$samples{$sample}->{alt}}) );
			return unless ( is_wanted_ref($samples{$sample}->{ref}) );
			return unless ( is_wanted_gt($samples{$sample}->{gt}) );
            return if position_lies_within_known_cnv($child->{chr}, $child->{pos}, $in_hrefs->{cnv_filter_href});
			$gt { $sample . "_gt" } = $samples{$sample}->{gt};
		}
		$gt { chr }    = $child->{chr};
		$gt { pos }    = $child->{pos};
		$gt { ref }    = $child->{ref};
		$gt { alt }    = $child->{alt}->[0];
		my $num_of_uniq_alt_positions = map {$_=()} split "", join ("", map { @{$_->{alt}} } $child, $mom, $dad);
		return if $num_of_uniq_alt_positions != 3;
		return if all { $_->{id} eq 'common' } $child, $mom, $dad;
	}
	return if ( ($gt{child_gt} . $gt{mom_gt} . $gt{dad_gt}) =~ /^0+$/  ); # exclude monomorphic sites or missing data

    my $event_status = get_event_status_from_genotypes( $gt{child_gt}, $gt{mom_gt}, $gt{dad_gt} );

	my $event_list_output_string = "$gt{chr}\t$gt{pos}\t$gt{ref}\t$gt{alt}\t$gt{child_gt}\t$gt{mom_gt}\t$gt{dad_gt}\t$event_status";
	print {$fh_href->{events}} $event_list_output_string,"\n" unless ($opt{no_event_list});

	return { event_status => $event_status, chr => $gt{chr} }

}

sub tabulate_events {
    my $events_href = shift;
    my %tabulate_event;
	
	my @informative_event_types = qw(UA_P UI_P UA_M UI_M BPI informative);
    # We need to tally results as a set of Bernulli Trials
    # Where success and failure are based on events that are informative for inheritance
    # We have four test cases:
        # uniparental maternal isodisomy
        # uniparental maternal heterdisomy
        # uniparental paternal isodisomy
        # uniparental paternal heterodisomy
    # And success is one of these conditions, failure is an informative site that is not this one condition
    # Proportion is (probability of success) is cases/(cases+failures)  = cases/(informative_sites)

	my %all_chrs = map { $_ => 0 } @informative_event_types;
	my @chrs = sort { $a <=> $b } keys %{$events_href};
    print Dumper $events_href if $opt{testing};
    for my $chr (@chrs)    {
        # Preparing binomial tests; convert event types to tallies of heterodisomy and isodisomy, maternal and paternal
		$events_href->{$chr}->{informative} = sum map { $events_href->{$chr}->{$_} }  ( grep { $_ ne "uninformative" } keys %{$events_href->{$chr}} ) ;
        $events_href->{$chr}->{mat_hUPD} += $events_href->{$chr}->{UA_M} if exists $events_href->{$chr}->{UA_M};
        $events_href->{$chr}->{mat_iUPD} += $events_href->{$chr}->{UA_M} if exists $events_href->{$chr}->{UA_M};
        $events_href->{$chr}->{mat_iUPD} += $events_href->{$chr}->{UI_M} if exists $events_href->{$chr}->{UI_M};
        $events_href->{$chr}->{pat_hUPD} += $events_href->{$chr}->{UA_P} if exists $events_href->{$chr}->{UA_P};
        $events_href->{$chr}->{pat_iUPD} += $events_href->{$chr}->{UA_P} if exists $events_href->{$chr}->{UA_P};
        $events_href->{$chr}->{pat_iUPD} += $events_href->{$chr}->{UI_P} if exists $events_href->{$chr}->{UI_P};

		$events_href->{all_chrs}->{informative}     += $events_href->{$chr}->{informative}  if exists $events_href->{$chr}->{informative};      
		$events_href->{all_chrs}->{mat_hUPD}        += $events_href->{$chr}->{mat_hUPD} if exists $events_href->{$chr}->{mat_hUPD};
		$events_href->{all_chrs}->{mat_iUPD}        += $events_href->{$chr}->{mat_iUPD} if exists $events_href->{$chr}->{mat_iUPD};
		$events_href->{all_chrs}->{pat_hUPD}        += $events_href->{$chr}->{pat_hUPD} if exists $events_href->{$chr}->{pat_hUPD};
		$events_href->{all_chrs}->{pat_iUPD}        += $events_href->{$chr}->{pat_iUPD} if exists $events_href->{$chr}->{pat_iUPD};
		$events_href->{all_chrs}->{BPI}             += $events_href->{$chr}->{BPI}  if exists $events_href->{$chr}->{BPI};      
	}
    my @binomial_test_types = qw(pat_hUPD pat_iUPD mat_hUPD mat_iUPD BPI informative);
    print {$fh_href->{table}} "Tally of Signatures for Binomial Testing\n";
	print {$fh_href->{table}} "chr\t", join ("\t", @binomial_test_types),"\n";		
	for my $chr ($chrs[0]..$chrs[$#chrs]) {
		print {$fh_href->{table}} "$chr\t";
		for (@binomial_test_types) {
			if (! $events_href->{$chr}->{$_} ) { print {$fh_href->{table}} "0\t" } 
			else {
				print {$fh_href->{table}} $events_href->{$chr}->{$_},"\t"
			}
		}
		print {$fh_href->{table}} "\n";
    }
	print {$fh_href->{table}} "all_chrs\t";
	print {$fh_href->{table}} join ("\t", map {$events_href->{all_chrs}->{$_}} @binomial_test_types) ,"\n";

    print {$fh_href->{table}} "\n\n";

    print {$fh_href->{table}} "Tally of Event Types\n";
	print {$fh_href->{table}} "chr\t", join ("\t", @informative_event_types),"\n";		
	for my $chr ($chrs[0]..$chrs[$#chrs]) {
		print {$fh_href->{table}} "$chr\t";
		for (@informative_event_types) {
			if (! $events_href->{$chr}->{$_} ) { print {$fh_href->{table}} "0\t" } 
			else {
				print {$fh_href->{table}} $events_href->{$chr}->{$_},"\t"
			}
		}
		print {$fh_href->{table}} "\n"
	}
	return $events_href;
}

sub run_binomial_test_on_tabulated_events_to_get_significant_results {
    my $tab_events_href = shift;
    print {$fh_href->{log}} "Running binomial tests...\n";
    my @binom_test_output;
    my @chrs_to_test = grep { $_ ne "all_chrs" } (keys %$tab_events_href);
    my @events_to_test = qw(pat_hUPD pat_iUPD mat_hUPD mat_iUPD);
    my $signif_level = get_significance_level();
    my $R = Statistics::R->new(r_bin => $opt{path_to_R});
    for my $chr (@chrs_to_test) {
        for my $event_type (@events_to_test) {
            next if ! $tab_events_href->{$chr}->{$event_type};
            my $event_avg_per_chr = $tab_events_href->{all_chrs}->{$event_type} / $tab_events_href->{all_chrs}->{informative};
            my $success = $tab_events_href->{$chr}->{$event_type};
            my $failure = $tab_events_href->{$chr}->{informative} - $tab_events_href->{$chr}->{$event_type};
            $R->set('success', $success); 
            $R->set('failure',$failure); 
            $R->set('prob_of_success', $event_avg_per_chr);
            $R->send('binom.test(x=c(success,failure),p=prob_of_success,alternative="g")$p.value');
            (my $binom_test_pval = $R->read) =~ s/.*( )(.*)/$2/;
            next if (  ($binom_test_pval > $signif_level) &&  (! $opt{testing}) );
            push (@binom_test_output , [ $binom_test_pval, $event_type , $chr ]);
        }
    }
    if (scalar @binom_test_output == 0) {
        print {$fh_href->{log}} "No significant results found\n";
        return
    }
    return(\@binom_test_output);
}

sub print_upd_results {
    my $results_aref = shift;
    my @results = sort { $a->[0] <=> $b->[0] } @$results_aref;
    for (@results) {
        my ($sig, $type, $chr) = ( @{$_} );
        print {$fh_href->{upd_results}} join ("\t", $sig, $type, $chr),"\n"
    }
    return 1
}
       
sub plot_upd_results {
   my ($child_name,$paths) = @_;
   return if $opt{no_plot};
   print {$fh_href->{log}} "Plotting positions and their events...\n";
   my $events_list_path = $fh_href->{events_path};
   my $output_dir = $opt{output_path};
   my $plot_sex = "FALSE";
   $plot_sex = "TRUE" if $opt{include_X};
   my $multi = "FALSE";
   my $vcf = $paths->{child_path};
   if ($opt{multisample_vcf}) { 
        print {$fh_href->{log}} "Plotting Multisample VCF in R\n";
        $multi = "TRUE";
		$vcf = $paths->{multisample_path}
   }
    print {$fh_href->{log}} "Running R\n";
    print {$fh_href->{log}} Dumper $opt{path_to_R}, "--slave -f", $opt{R_scripts_dir},
         "/plot_zygosity_and_events.R --args", $events_list_path, $vcf, $output_dir, $plot_sex, $multi;	
    `$opt{path_to_R} --slave -f $opt{R_scripts_dir}/plot_zygosity_and_events.R --args $events_list_path $vcf $output_dir $plot_sex $multi`;	
    return 1;
}


# Helper Subs

sub is_wanted_ref {
	my $ref = shift;
	return 1 if $ref =~ m/^[ACGT]$/;
}

sub is_wanted_alt {
	my $alt = shift;
	return 1 if (scalar @{$alt} == 1 && $alt->[0] =~ m/^[ACGT]$/);
}

sub is_wanted_gt {
	my $gt = shift;
	return 1 if $gt =~ m{^[01][01]$};
}

sub is_wanted_filter {
	my $filter = shift;
	return 1 if ( $filter =~ m/(^[.]$)|(PASS)/ );
}

sub get_significance_level {
    my $p_val;
    if ($opt{significance_level}) {
        $p_val = $opt{significance_level};
        print {$fh_href->{log}} "User defined p_value cut off of $p_val will be used\n";
    } else {
        my $num_chr = 22;
        $num_chr = 23 if ($opt{include_X});
        my $number_of_tests = $num_chr * 4;
        my $alpha = 0.05;
        $p_val = $alpha/$number_of_tests;
        print {$fh_href->{log}} "Testing significance at p-value: $alpha/$number_of_tests = ", $alpha/$number_of_tests, ".\n";
    }
    return $p_val
}

sub exit_routine {
    close_filehandles();   
    print {$fh_href->{log}} "UPDio finished successfully\n";
    exit();
}

sub close_filehandles {
    for (keys %{$fh_href}) {
        if (ref $_) {
            $fh_href->{$_}->close || die "Problem closing $_\n";
        }
    }
}
        
sub parse_vcf {
    my $path = shift;
    print {$fh_href->{log}} "Reading genotypes for $path...\n";
    my $vcf = Vcf->new( file => $path );
    $vcf->parse_header();
    return iter sub {
        my $x_href = $vcf->next_data_hash() or return;
		my $chr = chr2num($x_href->{CHROM});
		if (! $chr) {
            print Dumper $x_href; die "FATAL: I do not undersatnd chromosome identity in this line.$x_href->{CHROM}\n";
		}
        if (! $x_href->{ID} ) {
            print Dumper $x_href; die "FATAL: I do not understand ID in this line $x_href->{ID}.\n";
        }
		# Handle multisample vcf
        if (scalar (keys %{ $x_href->{gtypes} }) == 3) {
            return { 
                chr    => $chr,
                pos    => $x_href->{POS},
                ref    => $x_href->{REF},
                alt    => $x_href->{ALT},
                child_gt	=> format_gt($x_href->{gtypes}->{$opt{childID}}->{GT}),
                mom_gt		=> format_gt($x_href->{gtypes}->{$opt{momID}}->{GT}),
                dad_gt		=> format_gt($x_href->{gtypes}->{$opt{dadID}}->{GT}),
                filter => $x_href->{FILTER},
                id     => $x_href->{ID}
			}
        };
		# Else handle 3 single-sample vcfs
        if ($x_href->{ID} eq "common") { 
            return {
                chr    => $chr,
                pos    => $x_href->{POS},
                ref    => $x_href->{REF},
                alt    => $x_href->{ALT},
                gt     => "00",
                filter => $x_href->{FILTER},
                id     => "common"
            }
        } else {
            my $gt = ( values %{ $x_href->{gtypes} } )[0]->{GT}; # SNV genotypes should be encoded by GT
            if (! $gt) { 
                	#print "no genotype found!\n", Dumper $path, $x_href;
					#die Dumper "problem with genotype\n", $x_href, $gt;
					$gt = "./.";
					$x_href->{FILTER} = [ "FAIL" ];
                	#print "changing x_href!\n", Dumper $gt, $x_href;
            }
            return { 
                chr    => $chr,
                pos    => $x_href->{POS},
                ref    => $x_href->{REF},
                alt    => $x_href->{ALT},
                gt     => format_gt($gt),
                filter => $x_href->{FILTER},
                id     => $x_href->{ID}
            };
        }
    }
}

sub format_gt {
	my $gt = shift;
	$gt =~ s{[/|]}{};
	return $gt;
}

sub chr2num {
    my $chr = shift;
    $chr =~ s/^[Cc]hr//;
    return $chr if $chr =~ m/\d+/;
    return 23 if $chr =~ /^X$/;
    return 24 if $chr =~ /^Y$/;
    return 25 if $chr =~ /^M$/;
    die "Check chromosome name: $chr \n";
}

sub cmp_pos {
    my ( $x, $y ) = @_;
    # If the chromosomes are the same, then compare the positions.  
    # <=> comparison operator; returns 0 if x == y, 1 if x > y, -1 if x < y
    # ( ex: 14 <=> 14 = 0.  The left side is 0, so the || then compares the right side )
    if (! $x || ! $y) {
        print Dumper $x, $y; die "FATAL: I do not understand comparison\n";
    }
    ;
    if ($x->{chr} - $y->{chr} > 1) {
        die "$x->{chr} and $y->{chr} !!  Files not sorted correctly\n";
    }
    ;
    if ($y->{chr} - $x->{chr} > 1) {
        die "$x->{chr} and $y->{chr} !!  Files not sorted correctly\n";
    }
    ;
    return $x->{chr} <=> $y->{chr} || $x->{pos} <=> $y->{pos};
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

sub translate_nomenclature {
    my $before = shift;
    my %trans = (
        'iUPI_P' => "UI_P",
        'hUPI_P' => "AI_P",
        'iUPI_M' => "UI_M",
        'hUPI_M' => "UI_M",
        'hUPI_M,iUPI_M' => "AI_M",
        'hUPI_P,iUPI_P' => "AI_P",
        BPI => "B"
    );
    if (exists $trans{$before}) {
        return $trans{$before}
    } else {
        return $before
    }
}

sub calculate_concordance {
    my ( $x, $y ) = @_;
    if (! $x || ! $y ) {
        print Dumper $x, $y; die "FATAL: I do not understand comparison being made in calculate_concordance\n";
    }
    ;
    my $x1 = substr($x,0,1);
    my $x2 = substr($x,1,1);
    my $y1 = substr($y,0,1);
    my $y2 = substr($y,1,1);

    if ( ( $x1 eq $y1 && $x2 eq $y2 ) || ( $x1 eq $y2 ) && ( $x2 eq $y1 ) ) {
        return 1;
    } elsif ( ( $x1 eq $y1 ) || ( $x1 eq $y2 ) || ( $x2 eq $y1 ) || ( $x2 eq $y2 ) ) {
        return 0.5;
    } else {
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

sub get_event_status_from_genotypes {
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
    # Informative sites:  Uniparental Isodisomy (UI) and Uniparental Ambiguous (UA)
    # Informative sites:  Bi-parental inheritance (BPI)
    # homozygous = zygosity of 1, heterozygosity = zygosity of 2
    my @event_lookup_table = (
        {child_zygos => 2, mom_zygos => 1, dad_zygos => 1, mom_dad => 0,   child_mom => 0.5, child_dad => 0.5, value => 'BPI'},
        {child_zygos => 1, mom_zygos => 2, dad_zygos => 1, mom_dad => 0.5, child_mom => 0.5, child_dad => 0,   value => 'UI_M'},
        {child_zygos => 1, mom_zygos => 1, dad_zygos => 1, mom_dad => 0,   child_mom => 0,   child_dad => 1,   value => 'UA_P'},
        {child_zygos => 1, mom_zygos => 1, dad_zygos => 1, mom_dad => 0,   child_mom => 1,   child_dad => 0,   value => 'UA_M'},
        {child_zygos => 1, mom_zygos => 1, dad_zygos => 2, mom_dad => 0.5, child_mom => 0,   child_dad => 0.5, value => 'UI_P'}
    );

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
    if (! $cnv_positions_in_child_href) {
        die "Problem loading cnv data\n";
    }
    ;
    my $wobble = 10000; # Add 5kb padding to all CNVs; now trying 20kb padding because I saw increased error rates (~15% between exome samples and simulations showed decreased accuracy); yes 10kb on both sides was found to work well, and match previous simulations
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
			# This is a way to create a nested hash structure so I don't have to search through too many CNVs for each position.
			# Another way to do this is to create an iterator and read it along the file, but this is quite fast.
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
    open ($common_cnv_fh, "<", $common_cnv_file) or 
        die "$!: Can't open common cnv file common_cnv_file: $opt{common_cnv_file}; did you specify a proper location with the appropriate option?\n";
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
        open ($cnv_fh, "<", $sample_cnv_path) or die "Yikes! can't read common cnv file $sample_cnv_path!\n";
        print {$fh_href->{log}} "Sample-specific CNV data found for proband\n";	
        if ($sample_cnv_path =~ /vcf$/ || $sample_cnv_path =~ /vcf\.gz$/) {
            # parse CNV data as VCF file
            my $vcf = Vcf->new( file => $sample_cnv_path) ;
            $vcf->parse_header();
            while ( my $x_href = $vcf->next_data_hash() ) {	
				# This is how our CNV vcf files are encoded.
                next if $x_href->{INFO}->{NUMBERALGORITHMSCNSOLIDATE} < 4;
                my ($cnv_chr, $cnv_start, $cnv_stop)  =  ($x_href->{CHROM},$x_href->{POS}, $x_href->{INFO}->{END});
                $cnv_chr =~ s/[Cc]hr//;
                my $cnv_copies = $x_href->{INFO}->{SVTYPE};
                $cnv_copies =~ s/DUP/3/;
                $cnv_copies =~ s/DEL/1/;

                push @{$sample_cnvs_href->{$cnv_chr}}, "$cnv_start,$cnv_stop";
            }
        } else {
            # parse CNV data as tab separated file (chr start stop cn)
            while (my $cnv_line = <$cnv_fh>) {
                next if ($cnv_line =~ /^#/);
                my @fields = split /\t/, $cnv_line;
                my ($cnv_chr, $cnv_start, $cnv_stop, $cnv_copies) = @fields;
                $cnv_chr =~ s/[Cc]hr//;
                push @{$sample_cnvs_href->{$cnv_chr}}, "$cnv_start,$cnv_stop";
            }
        }
    } else {
        print {$fh_href->{log}} "WARNING: no proband cnv data available; deletions can mimic upd events\n";
    }
    return($sample_cnvs_href);
}


## END SUB

__END__

=over 25

=head1 NAME

 UPDio.pl

=head1 VERSION

 Version 1.0

=head1 DESCRIPTION

 UPDio is designed to detect UPD events in probands from vcf trio data.

=head1 SYNOPSIS

 perl UPDio.pl --multisample_vcf <3-sample_vcf_file> --childID <childID> --momID <momID> --dadID <dadID>

 OR

 perl UPDio.pl --child_vcf <child_file> --mom_vcf <mom_file> --dad_vcf <dad_file> [ --child_cnv_data <child_cnv_data_file>  ]

=head2 Required Arguments

=item C<multisample_vcf>

a multi-sample vcf file with child, mother, and father genotypes, must include samples IDs, IDs must match samples in vcf header

 OR

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

=item C<include_MI>

for searching for mendelian inconsistent sites; occassionally useful for finding loss of transmitted allele

=back

=head1 EXAMPLE

 perl UPDio.pl --child_vcf Family101_child.vcf.gz --mom_vcf Family101_mom.vcf.gz --dad_vcf Family101_dad.vcf.gz --child_cnv_data Family101_child_cnvs.txt

 perl UPDio.pl --multisample_vcf multisample_trios/SNP_merged_example.vcf.gz --childID MAIN528 --dadID MAIN525 --momID MAIN526


=head1 AUTHORS

 Dan King  < dk6 at sanger.ac.uk >
 Coding oversight by Ray Miller
 Project oversight by Matt Hurles

=head1 INSTITUTION

 Wellcome Trust Sanger Institute
 Wellcome Trust Genome Campus,
 Hinxton, Cambridge, CB10 1SA, UK

=head1 LICENSE

  Copyright (c) 2014 Genome Research Ltd.
  Author: Dan King
  
  Any redistribution or derivation in whole or in part including any
  substantial portion of this code must include this copyright and
  permission notice.
  
  THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  
  This code is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at
  your option) any later version (http://www.gnu.org/copyleft/gpl.txt).

=cut
