#!/usr/bin/perl -w
use strict;
use bin::vcf2;
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use lib dirname (__FILE__) . "/bin";

my %options= ();
getopts("hc:st:b:v:e:o:", \%options);

if ($options{h}) { my $send = options_manager();}

print STDERR "./filter.pl ";
# HANDLE VCF ###############################################################
my $vcf_file;
if ($options{v}) { 
	$vcf_file = $options{v}; 
	print STDERR "-v $vcf_file ";
}
else { 
	print STDERR "MISSING VCF\n";
	exit;
}
my $vcf = CMD::vcf2->new('file'=> $vcf_file );
############################################################################


# HANDLE BAM ###############################################################
my $bam;
my $depth_cutoff;
if (defined $options{b}) { 
	if (defined $options{c}) {
		$bam = $options{b};
		$depth_cutoff = $options{c};
		print STDERR "-b $bam -c $depth_cutoff ";
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
	print STDERR "-e $exac_cutoff "
} #Lägg till kontroll för numeric value
else {
	print STDERR "-e $exac_cutoff(ExAC Default) "
}
############################################################################

# TUMOR CONC ###############################################################
my $tumor_conc;
if (defined $options{t}) {
	$tumor_conc = $options{t};
    print STDERR "-t $tumor_conc "
} #Lägg till kontroll för numeric value
else {
	print STDERR "Need to define tumor concentration\n";
	exit;
}
############################################################################

# SYNONYMOUS MUT? ##########################################################
my $syn = 0;
if (defined $options{s}) {
	print STDERR "-s (non-syn) ";
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
print STDERR "\nLOADING TUMOR SUPPRESSOR GENES $tumor_sup\n";
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
print STDERR "LOADING COSMIC VARIANTS $cosmic\n";
while (my $row = <$fh5>) {
	if ( $row =~ /^([0-9XY]+)\t(\d+)\t\S+\t(\w+)\t(\w+)/) {
		$COSMIC{"$1:$2"} = "$4";
	}
}
close $fh5;
############################################################################



my $CF = 0;
my $TGF = 0;
my $EF = 0;
my $VF = 0;
my $SF = 0;
### MAINSCRIPT ENSUES ###
my $c = 0;
 while ( my $a = $vcf->next_var() ) {
    $c++;
    #print Dumper($a);
    my $ref = $a->{REF};
 	my $alt = $a->{ALT};
    my $pos = $a->{POS};
    my $chrom =  $a->{CHROM};
    #print "$c PEN_______________________ $chrom\t$pos  __________________________\n";
    #if ($c > 10 ) {exit;}
    ## VAF ###########################################################
    ## id of sample is saved as only key for $a->{GT} need to retrieve
    my @keys = keys %{$a -> {GT}};
	my $RO = $a->{GT} -> {$keys[0]} -> {RO};
	my $AO = $a->{GT} -> {$keys[0]} -> {AO};
    my $VAF = $AO/($RO+$AO);
	my $QR = $a->{GT} -> {$keys[0]} -> {QR};
	my $QA = $a->{GT} -> {$keys[0]} -> {QA};
	my $QpA = $QA/$AO;
	if ($a->{QUAL} < 3000 ) {#print "$RO\t$AO\t$VAF\t$QpA\n";
		next;
	}
    if(! defined $VAF) {next;}
    my $check_vaf = vaf($VAF,$tumor_conc);
    if ($check_vaf == 0) {
        $VF++;
        #print "VAF!\n";
        next;
    }
    ##################################################################

    ## REMOVE COSMIC VARIANTS ###########
	if (defined $COSMIC{"$chrom:$pos"}) { 
    	if ($alt eq $COSMIC{"$chrom:$pos"} ) {
            $CF++;
            #print "COSMIC!\n";
            next;
        }
    }
    ######################################

    # FIND CONSENSUS CONSEQUENCE AND SCORE ######
    my $csq = $a->{INFO}->{CSQ};
    my ($max_conq, $score_conq) = CSQ($alt,$csq);
    #############################################

    ## SYNONYMOUS VARIANTS REMOVE IF OPTION -s ###
    if (defined $options{s}) {
        if ($score_conq > 2) {
            $SF++;
            #print "SYN! \n";
            next;
        }
    }
    ###############################################

    ## Check if tumor suppressor gene ###################
    my $hgnc_id = $a->{INFO}->{CSQ}->[0]->{HGNC_ID};
    my $is_tumor = 0;
    if (defined $tumor_SUPP{$hgnc_id}) { $is_tumor = 1;}
    if ($score_conq < 2 && $is_tumor == 1) {
        $TGF++;
        #print "TUM! \n";
        next;
    } 
    #####################################################

    ## ExAC cutoff ##############################################
    my $exac_maf = $a->{INFO}->{CSQ}->[0]->{ExAC_MAF};
    my $exac_check = control_exac($exac_maf,$exac_cutoff,$alt);
    
    if ($exac_check == 0) {
        $EF++;
        #print "EXAC!\n";
        next;
    }
    #############################################################
    my $variant = $chrom.":".$pos;

    if (defined $options{b}) {
        my $depth_filt = coverage($chrom, $pos, $bam, $depth_cutoff);
        if ($depth_filt == 0) { print   "DEPTH! \n"; }
    }

	$pass_filter++;
}

print STDERR "Cosmic $CF\t";
print STDERR "Tumor supp gene $TGF\t";
print STDERR "ExAC $EF\t";
print STDERR "VAF $VF\t";
if (defined $options{s}) {print STDERR "Synonymous $SF\t";}
print "PASSED $pass_filter/$coding_size\n";

sub CSQ {
	my ($alt,$csq,$syn) = @_;
    

    my %conq_tally_canon;
    my %conq_tally_nocanon;
    my %conq_tally;
    my $has_canon = 0;
    foreach my $b (@$csq) {
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
    }

    if ($has_canon == 1) {%conq_tally = %conq_tally_canon; }
    else {%conq_tally = %conq_tally_nocanon;}
    
    ## Find consensus consequence(most common from VEP, and if several equal most severe)
    my ($max_hash, $score_max) = max_hash(%conq_tally);

    return $max_hash, $score_max;

}

sub vaf {
	my ($VAF, $tumorc) = @_;
	my $check = 0;
	if ($VAF/$tumorc > 0.05 && $VAF/$tumorc < 0.4 ) {
			$check++; 
	}
	return $check;

}


sub coverage {
	my ($chrom, $pos, $bam, $depth_cutoff) = @_;
	

	my $depth_filt = 1;
	my $output = "samtools depth ";
	$output = $output.$bam." -r ".$chrom.":".$pos."-".$pos;
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
	my ($exac_maf,$exac_cutoff,$alt) = @_;
	my $exac_filt = 1;
	my $strict = 1;  ##Allele does not have to match exac allele, not implemented in %options

	if ($exac_maf) {
		if ($exac_maf =~ /[ACGT-]:/) {
			
			my @multi_allele = split('&',$exac_maf);
		
			if (@multi_allele > 1) {
						
				foreach my $var (@multi_allele) {
					my @exac_allele = split(':',$var);
					if ( $exac_allele[1] > $exac_cutoff) {
						$exac_filt = 0;
					}
				}						
			}
			else {
				my @exac_allele = split(':',$multi_allele[0]);
				if ( $exac_allele[1] > $exac_cutoff) {
						$exac_filt = 0;
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
	print "-h\tthis help message\n";
    print "-o\toutput filtered vcf\n";
	exit;
	


}


sub max_hash {
    my (%data) = @_;
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

    
    my $max;
    while (my ($key, $value) = each %data) {
        if (not defined $max) {
            $max = $key;
            next;
        }
        if ($data{$max} < $value) {
            $max = $key;
        }
        if ($data{$max} == $value) {
            if ($score{$key} < $score{$max}) {
                $max = $key;
            }
        }
    }
    my $score_max = $score{$max};
    return $max, $score_max;
}
 






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









