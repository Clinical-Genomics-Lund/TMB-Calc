#!/usr/bin/perl -w
use strict;
use CMD::vcf2;
use Data::Dumper;
use Getopt::Std;
use List::UtilsBy qw(max_by);

my %options= ();
getopts("hc:st:b:v:e:o:", \%options);

if ($options{h}) { my $send = options_manager();}
my $count_canon = 0;  ##Legacy
print STDOUT "./filter.pl ";
# # HANDLE VCF ##############################################################
# my $filename; #REQUIRED
# if ($options{v}) { 
# 	$filename = $options{v}; 
# 	print STDOUT "$filename ";
# }
# else { 
# 	print STDERR "MISSING VCF\n";
# 	exit;
# }
# my @variants;
# open (my $fh, '<', $filename)
# 	or die "Could not open file $filename";
# while (my $line = <$fh>) {
# 	chomp $line;
# 	push @variants, $line;
# }
# close $fh;
# ############################################################################

# HANDLE BAM ###############################################################
my $bam;
my $depth_cutoff;
if (defined $options{b}) { 
	if (defined $options{c}) {
		$bam = $options{b};
		$depth_cutoff = $options{c};
		print STDOUT "-b $bam -c $depth_cutoff ";
	}  
	else {
		print STDERR "missing coverage -c value\n"; 
		exit;
	}
}
############################################################################

# HANDLE EXAC ##############################################################
my $exac_cutoff = 0.01; #DEFAULT VALUE
if (defined $options{e}) {
	$exac_cutoff = $options{e};
	print STDOUT "-e $exac_cutoff "
} #Lägg till kontroll för numeric value
else {
	print STDOUT "-(ExAC Default) $exac_cutoff \n"
}
############################################################################

# TUMOR CONC ###############################################################
my $tumor_conc;
if (defined $options{t}) {
	$tumor_conc = $options{t};
} #Lägg till kontroll för numeric value
else {
	print STDERR "Need to define tumor concentration\n";
	exit;
}
############################################################################

# SYNONYMOUS MUT? ##########################################################
my $syn = 0;
if (defined $options{s}) {
	$syn = 1;
}
############################################################################

# OUTPUT FILE IF CHOSEN ####################################################
my $outputvcf;
my $out;
if (defined $options{o}) {
	$outputvcf = $options{o};
	open(  $out, '>>', $outputvcf) or die "failure failure!";
}
############################################################################



# Some variables coding size needs to be incorporated as a real variable.
my $pass_filter = 0;
my $coding_size = 1.421418;
my $notfilt = 0;
my $filt = 0;


# TUMOR SUPPRESSOR GENES FILTER HASH #######################################
my $tumor_sup = "supporting_files/tumor_suppressor_genes.txt";
open (my $fh4, '<', $tumor_sup) or die "Could not open $tumor_sup\n";
my %tumor_SUPP;
print STDOUT "LOADING TUMOR SUPPRESSOR GENES $tumor_sup\n";
while (my $row = <$fh4>) {
	if ($row =~ /^(\d+)\t(\w+)\t/) {
		$tumor_SUPP{$1} = $2;
	}
}
close $fh4;
############################################################################
# COSMIC SOMATIC VARIANTS FILTER HASH ######################################
my $cosmic = "supporting_files/b37_cosmic_v54_120711.vcf";
open (my $fh5, '<', $cosmic) or  die "Could not open $cosmic";
my %COSMIC;
print STDOUT "LOADING COSMIC VARIANTS $cosmic\n";
while (my $row = <$fh5>) {
	if ( $row =~ /^([0-9XY]+)\t(\d+)\t\S+\t(\w+)\t(\w+)/) {
		$COSMIC{"$1:$2"} = "$4";
	}
}
close $fh5;
############################################################################

my $vcf_file;
if ($options{v}) { 
	$vcf_file = $options{v}; 
	print STDOUT "$vcf_file \n";
}
else { 
	print STDERR "MISSING VCF\n";
	exit;
}
my $vcf = CMD::vcf2->new('file'=> $vcf_file );

### MAINSCRIPT ENSUES ###
my $c = 0;
 while ( my $a = $vcf->next_var() ) {
     $c++;
     #print Dumper($a);
     #print $a->{INFO}->{CSQ}->[0]->{Consequence},"\n";

     #if( grep /^synonymous_variant$/, $a->{INFO}->{CSQ}->[0]->{Consequence} ) {

   


 	my $ref = $a->{REF};
 	my $alt = $a->{ALT};
    my $pos = $a->{POS};
    my $chrom =  $a->{CHROM};
    print "$c _______________________ $chrom\t$pos  __________________________\n";
    if ($c > 100 ) {exit;}
    ## VAF ###########################################################
    ## id of sample is saved as only key for $a->{GT} need to retrieve
    my @keys = keys %{$a -> {GT}};
    my $VAF = $a->{GT} -> {$keys[0]} -> {AF};
    if(! defined $VAF) {next;}
    print $VAF/$tumor_conc;
    my $check_vaf = vaf($VAF,$tumor_conc);
    print "\t",$check_vaf,"\n";
    if ($check_vaf == 0) {print "VAF!\n"; next;}
    ##################################################################

    ## REMOVE COSMIC VARIANTS ###########
	if (defined $COSMIC{"$chrom:$pos"}) { 
    	if ($alt eq $COSMIC{"$chrom:$pos"} ) {
            print "COSMIC!\n";
            next;
        }
    }
    ######################################
    ## Check if tumor suppressor gene ###################
    my $hgnc_id = $a->{INFO}->{CSQ}->[0]->{HGNC_ID};
    my $is_tumor = 0;
    if (defined $tumor_SUPP{$hgnc_id}) { $is_tumor = 1;}
    #####################################################
    my $consq = $a->{INFO}->{CSQ}->[0]->{Consequence};
    #foreach my $try (@$consq) {
    #    print $try,"\n";
    #}
    my $csq = $a->{INFO}->{CSQ};
    my $vep = CSQ($alt,$csq,$syn);
# 	my $depth_filt = 1;


    ## ExAC cutoff 


 	# if (defined $options{o}) {
 	# 	if ($ /^#/ ) {
 	# 		print $out "$_\n";
 	# 	}
 	# }


	#####################################
	my $variant = $a->{CHROM}.":".$a->{POS}.":".$alt;
	
	
}
#print $filt,"\t",$notfilt,"\n";


sub CSQ {
	my ($alt,$csq,$syn) = @_;
    
    my %score = ("transcript_ablation"=>1,
    "splice_acceptor_variant"=>1,
    "splice_donor_variant"=>1,
    "stop_gained"=>1,
    "frameshift_variant"=>1,
    "stop_lost"=>1,
    "start_lost"=>1,
    "transcript_amplification"=>1,
    "inframe_insertion"=>2,
    "inframe_deletion"=>2,
    "missense_variant"=>2,
    "protein_altering_variant"=>2,
    "splice_region_variant"=>3,
    "incomplete_terminal_codon_variant"=>3,
    "start_retained_variant"=>3,
    "stop_retained_variant"=>3,
    "synonymous_variant"=>3,
    "coding_sequence_variant"=>4,
    "mature_miRNA_variant"=>4,
    "5_prime_UTR_variant"=>4,
    "3_prime_UTR_variant"=>4,
    "non_coding_transcript_exon_variant"=>4,
    "intron_variant"=>4,
    "NMD_transcript_variant"=>4,
    "non_coding_transcript_variant"=>4,
    "upstream_gene_variant"=>4,
    "downstream_gene_variant"=>4,
    "TFBS_ablation"=>4,
    "TFBS_amplification"=>4,
    "TF_binding_site_variant"=>4,
    "regulatory_region_ablation"=>2,
    "regulatory_region_amplification"=>4,
    "feature_elongation"=>4,
    "regulatory_region_variant"=>4,
    "feature_truncation"=>4,
    "intergenic_variant"=>4 );
    my %conq_tally_canon;
    my %conq_tally_nocanon;
    my %conq_tally;
    my $has_canon = 0;
                #if (grep /synonymous_variant/, @$tmp) {
                #print $b -> {Allele}, $alt, "\n";
            #}
    foreach my $b (@$csq) {
        
        #if ($alt eq $b -> {Allele}) {
            ## if canon pick consensus conseqeunce or if equal most severe ^
            if ( $b -> {CANONICAL} eq "YES" ) {
                $has_canon = 1;
                
                my $tmp = $b->{Consequence};
                foreach my $conq (@$tmp) {
                    if ($conq_tally_canon{$conq}) {
                        $conq_tally_canon{$conq}++;
                    }
                    else {
                     $conq_tally_canon{$conq} = 1;
                    }
                }
            }
            ## if no canon pick consensus conseqeunce or if equal most severe ^
            else {
                
                my $tmp = $b->{Consequence};
                foreach my $conq (@$tmp) {
                    if ($conq_tally_nocanon{$conq}) {
                        $conq_tally_nocanon{$conq}++;
                    }
                    else {
                     $conq_tally_nocanon{$conq} = 1;
                    }
                }
            }
        #}
    }
    print scalar(keys %conq_tally_nocanon);
    if ($has_canon == 1) {%conq_tally = %conq_tally_canon; print "HAS CANON\n";}
    else {%conq_tally = %conq_tally_nocanon; print "NO CANON\n";}
    
    foreach my $key (keys %conq_tally) {

        print $key,"=", $conq_tally{$key},"\n";
    }

			# if ($in_hash{$key}[1] =~ /missense_variant/) {
			# }
			# elsif ($in_hash{$key}[1] =~ /stop_gained/) {
			# }
			# elsif ($syn == 1) { 
			# 	if ($in_hash{$key}[1] =~ /synonymous_variant/ ) {

	#return $check, $somatic;
}

sub vaf {
	my ($VAF, $tumorc) = @_;
	my $check = 0;
	if ($VAF/$tumorc > 0.05 && $VAF/$tumorc < 0.4 ) {
		#if ($array[7] > 0 ) {
			$check++; 
		#}
	}
	return $check;

}


sub coverage {
	my @in = @_;
	my@IN = split(":",$in[0]);
	my $depth_cutoff = $in[2];
	my $depth_filt = 1;
	my $output = "samtools depth ";
	$output = $output.$in[1]." -r ".$IN[0].":".$IN[1]."-".$IN[1];
	my $cmd = `$output`;
	chomp $cmd;	
	$cmd =~ s/^[0-9XY]+\s+\d+\s+//;
	if ($cmd < $depth_cutoff) {
		$depth_filt = 0;
		print "$cmd\t$depth_cutoff\n";
	}
	
	return $depth_filt;

}


sub control_exac {
	my ($exac_maf,$alt,$ref,$exac_cutoff,$csq) = @_;
	my $exac_filt = 1;
	my $strict = 1;  ##Allele does not have to match exac allele, not implemented in %options

	if ($exac_maf) {
		if ($exac_maf =~ /[ACGT-]:/) {
			
			my @multi_allele = split('&',$exac_maf);
		
			if (@multi_allele > 1) {
						
				foreach my $var (@multi_allele) {
					my @exac_allele = split(':',$var);
					if ($strict == 1) { 
						if ( $exac_allele[1] > $exac_cutoff) {
							$exac_filt = 0;
						}
					}
					elsif ($alt eq $csq && $alt eq $exac_allele[0]) {
						if ( $exac_allele[1] > $exac_cutoff) {
							$exac_filt = 0;
						}
					}
				}						
			}
			else {
				my @exac_allele = split(':',$multi_allele[0]);
				if ($strict == 1) {
					if ( $exac_allele[1] > $exac_cutoff) {
							$exac_filt = 0;
						}
				}
				elsif ($alt eq $csq && $alt eq $exac_allele[0]) {
					if ( $exac_allele[1] > $exac_cutoff) {
						$exac_filt = 0;
					}
				}
			}
		}		
	}



	
	return $exac_filt;




}

sub options_manager {
	#my $in = @_;

	print "CMD\tDescription\n";
	print "-v\tinput VCF REQUIRED\n";
	print "-b\tinput BAM\n  -c\tcoverage cutoff REQUIRED if -b in use\n";
	print "-t\ttumour concentration REQUIRED\n";
	print "-s\tinclude Synonymous mutations\n";
	print "-e\tExAC MAF cutoff. Default 0.01\n";
	print "-n\tPrint new VCF, OUTPUTPATH REQ\n";
	print "-h\tthis help message\n";
	exit;
	


}

print $pass_filter/$coding_size,"\n";








# # TALLY ALL SAMPLES IN POOL (LEGACY?) ######################################
# my @all_samples;
# my $list_samples = "samples.txt";
# open (my $fh2, '<', $list_samples)
# 	or die "Could not open file $list_samples";
# while (my $line = <$fh2>) {
# 	chomp $line;
# 	push @all_samples, $line;
# }
# close $fh2;
# my %all;
# foreach my $file (@all_samples) {
# 	chomp $file;
# 	open (my $fh3, '<', $file)
# 		or die "Could not open file $file";
# 	if ($file =~ /\S+\/(\S+)\.vcf/) { $file = $1; }
# 	while (my $line = <$fh3>) {
# 		chomp $line;
# 		if ($line =~ /^([0-9XY]{1,2})\t(\d+)\t(\S+|\.)\t(\S+)\t(\S+)\t[\.0-9]+\t/) {
		
# 			my $variant = $1.":".$2.":".$5;
# 			if (defined $all{$variant}) {
# 				$all{$variant}++;
# 			}
# 			else {
# 				$all{$variant} = 0;
# 			}	
# 		}
# 	}
# close $fh3;
# }
# ############################################################################









