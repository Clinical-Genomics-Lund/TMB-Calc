#!/usr/bin/perl -w
use strict;
use bin::vcf2;
use Data::Dumper;
use Getopt::Std;
use File::Basename;
use lib dirname (__FILE__) . "/bin";

my %options= ();
getopts("h:st:d:v:e:o:a:", \%options);

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
if (defined $options{d}) { 
	$depth_cutoff = $options{d};
	print STDERR "-d $depth_cutoff ";
}
############################################################################

# HANDLE EXAC ##############################################################
my $exac_cutoff = 0.0001; #DEFAULT VALUE
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

if (defined $options{o}) {
	$outputvcf = $options{o};
	open(  OUT, '>', $outputvcf) or die "failure failure!";
}
############################################################################

# CALLER ###################################################################
my $caller;
if (defined $options{a}) {
	$caller = $options{a};
	my %valid = ("sent" => "1", "free" => "1");
	if (! defined $valid{$caller}) { 
		print STDERR "WRONG CALLER!!!, valid callers are sent and free\n"; 
		exit;
	}
}
############################################################################

# SGZ FILES ################################################################
my $mut_agg = "test_runs_output/$caller.mutagg";
open (MUT, '>', $mut_agg) or die $!;
print MUT "\#sample\tmutation\tfrequency\tdepth\tpos\tstatus\tstrand\teffect\n";
#
my $purity = "test_runs_output/$caller.purity";
open (PUR, '>', $purity) or die $!;
print PUR $tumor_conc;
close PUR;
############################################################################

# Some variables coding size needs to be incorporated as a real variable. ##
my $pass_filter = 0;
my $coding_size = 1.421418;
############################################################################

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

my $ct = 0;
my $count_indel = 0;
my $CF = 0;
my $TGF = 0;
my $EF = 0;
my $VF = 0;
my $SF = 0;
my $DF = 0;
my $QCF = 0;
### MAINSCRIPT ENSUES ###
my $c = 0;
 while ( my $a = $vcf->next_var() ) {
    $c++;
    
	#print Dumper($a);
    my $ref = $a->{REF};
 	my $alt = $a->{ALT};
    my $pos = $a->{POS};
    my $chrom =  $a->{CHROM};
	my $strand = $a->{INFO}->{CSQ}->[0]->{STRAND};
	if ($strand > 0) { $strand = "+";}
	elsif ($strand < 0) { $strand = "-";}
	my $hgnc_id = $a->{INFO}->{CSQ}->[0]->{HGNC_ID};
	if (length($ref) > 1|| length($alt) > 1) {$count_indel++;}
    #print "$c PEN_______________________ $chrom\t$pos  __________________________\n";
    #if ($c > 10 ) {exit;}
    ## VAF ###########################################################
    ## id of sample is saved as only key for $a->{GT} need to retrieve
	my $VAF;
	my @keys = keys %{$a -> {GT}};
	my $DEPTH; 
	## SENTIEON
	if ($caller eq "sent") {
		$VAF = $a->{GT} -> {$keys[0]} -> {AF};
		my $QC = $a->{GT} -> {$keys[0]} -> {BaseQRankSumPS};
		$DEPTH = $a->{GT} -> {$keys[0]} -> {AFDP};
		if (defined $QC) {
			if($QC < 0) {
				next;
			}
		}
		elsif ($a->{QUAL} < 500 ) {
			$QCF++;
			next;
		}
	}
	## FREEBAYES
	elsif ($caller eq "free") {
		my $RO = $a->{GT} -> {$keys[0]} -> {RO};
		my $AO = $a->{GT} -> {$keys[0]} -> {AO};
    	$VAF = $AO/($RO+$AO);
		$DEPTH = $a->{GT} -> {$keys[0]} -> {DP};
		## QUALITY CUTOFF ########################################
		my $QR = $a->{GT} -> {$keys[0]} -> {QR};
		my $QA = $a->{GT} -> {$keys[0]} -> {QA};
		my $QpA = $QA/$AO;
	
		if ($QpA < 30  ) {
			$QCF++;
			next;
		}
	}


    # FIND CONSENSUS CONSEQUENCE AND SCORE ######
    my $csq = $a->{INFO}->{CSQ};
    my ($max_conq, $score_conq) = CSQ($alt,$csq);
    #############################################


	##########################################################
    if(! defined $VAF) {next;}
	$VAF = $VAF/$tumor_conc;
    my $check_vaf = vaf($VAF,$tumor_conc);
    if ($check_vaf == 0) {
        $VF++;
        #print "VAF!\n";
        next;
    }
    ##################################################################





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
    my $is_tumor = 0;
    if (defined $tumor_SUPP{$hgnc_id}) { $is_tumor = 1;}
    if ($score_conq < 2 && $is_tumor == 1) {
        $TGF++;
        #print "TUM! \n";
        next;
    } 
    #####################################################

    ## REMOVE COSMIC VARIANTS ###########
	if (defined $COSMIC{"$chrom:$pos"}) { 
    	if ($alt eq $COSMIC{"$chrom:$pos"} ) {
            $CF++;
            #print "COSMIC!\n";
            next;
        }
    }
    ######################################

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

	## DEPTH ########################################################
    if (defined $options{d}) {
        #my $depth_filt = coverage($chrom, $pos, $bam, $depth_cutoff);
        #if ($depth_filt == 0) { print   "DEPTH! \n"; }
		if ($DEPTH < $depth_cutoff ) {
			$DF++;
			next;
		}

    }
	#################################################################

	## COUNT C>T SUBSTITIONS ###############
	if ($ref eq "C" && $alt eq "T") {$ct++;}
	########################################
	
	## PASSED #####
	$pass_filter++;
	###############

	## WRITE SGZ OUTPUT #######################################################
	my $sgz = "$keys[0]\t$hgnc_id:DUMMY:$ref.$pos>$alt\_p.DUMMY\t";
	$sgz = $sgz."$VAF\t$DEPTH\tchr$chrom:$pos\tunknown\t$strand\t$max_conq\n";
	print MUT $sgz;
	###########################################################################

	## BED OUTPUT ########################################
	if (defined $options{o}) {
		print OUT "$chrom\t$pos\t$pos\t$VAF\t$max_conq\n";
	}
	######################################################

}

## SGZ FILTERING ##
my $sgz_command = "/data/bnf/proj/twist/SFZ/basicSGZ.py test_runs_output/$caller.mutagg -f $caller.purity";
print STDERR "RUNNING SGZ, calculating somatic and germline variants\n";
my $run_command = `$sgz_command`;
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##






################################
print STDERR "Cosmic $CF\t";
print STDERR "Tumor supp gene $TGF\t";
print STDERR "ExAC $EF\t";
print STDERR "VAF $VF\t";
print STDERR "QC $QCF\t";
if (defined $options{d}) {print STDERR "DEPTH $DF\t";}
if (defined $options{s}) {print STDERR "Synonymous $SF\t";}
print STDERR "PASSED $pass_filter/$coding_size\tC>T: $ct\tindel:$count_indel\n";
################################


#							  #
#          SUBOUTINES         #
#							  #
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
	my ($VAF,$tumor_conc) = @_;
	my $check = 0;
	if ($VAF > 0.05 && $VAF < $tumor_conc ) {
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
 
