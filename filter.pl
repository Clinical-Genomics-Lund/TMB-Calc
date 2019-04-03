#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Std;
my %options= ();
getopts("hc:st:b:v:e:o:", \%options);

if ($options{h}) { my $send = options_manager();}
my $count_canon = 0;  ##Legacy
print STDOUT "./filter.pl ";
# HANDLE VCF ##############################################################
my $filename; #REQUIRED
if ($options{v}) { 
	$filename = $options{v}; 
	print STDOUT "$filename ";
}
else { 
	print STDERR "MISSING VCF\n";
	exit;
}
my @variants;
open (my $fh, '<', $filename)
	or die "Could not open file $filename";
while (my $line = <$fh>) {
	chomp $line;
	push @variants, $line;
}
close $fh;
############################################################################

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
	print STDOUT "-(ExAC Default) $exac_cutoff "
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
		$tumor_SUPP{$2} = $1;
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




### MAINSCRIPT ENSUES ###

foreach my $line (@variants) {
	if (defined $options{o}) {
		if ($line =~ /^#/ ) {
			print $out "$line\n";
		}
	}
	my @vep;
	my %vep_trans;
	my @keep_trans;
	my $ref;
	my $alt;
	my $depth_filt = 1;
	if ($line =~ /^([0-9XY]{1,2})\t(\d+)\t(\S+|\.)\t(\S+)\t(\S+)\t[\.0-9]+\t/) {
		$ref = $4;
		$alt = $5;
		## REMOVE COSMIC VARIANTS ###########
		if (defined $COSMIC{"$1:$2"}) { 
			if ($alt eq $COSMIC{"$1:$2"} ) {
				#next;
			}
		}
		#####################################
		my $variant = $1.":".$2.":".$5;
		if ($line =~ /(PASS)\t(ECNT=[\.0-9]+\;FS=[\.0-9]+\;HCNT=[\.0-9]+\;MAX_ED=[\.0-9]+\;MIN_ED=[\.0-9]+\;NLOD=[\.0-9]+\;NLODF=[\.0-9]+\;SOR=[\.0-9]+\;TLOD=[\.0-9]+\;CSQ=)(\S+)\t(.+$)/) {
			my $VEP = $3;
			my $FILT = $1;
			my $INFO = $2;
			my $SENTION = $4;
			@vep = split(",",$VEP);
			if (@vep > 0) {
				my $count = 0;
				foreach my $item (@vep) {
					$count++;
					my @trans = split('\|',$item);
					$vep_trans{$count} = [@trans];
				}
			}

			########################Send to sub##################
				my ($check_exac, $somatic) = exac(\%vep_trans,$alt,$ref, $exac_cutoff, $syn,\%tumor_SUPP);
				$count_canon = $count_canon + $somatic;
				my $check_vaf = vaf($SENTION,$tumor_conc);
				
				if ($check_exac == 1 && $check_vaf == 1) {

					if (defined $options{b}) {
						$depth_filt = coverage($variant, $bam, $depth_cutoff);
						if ($depth_filt == 1) {
							$pass_filter++;
						}
					}
					else {
						$pass_filter++;
					}
					if (defined $options{o}) {print $out "$line\n";}
				}
				#print "$somatic\n";
			#####################################################
		}
		else { #print $line,"\n";
		}
	}
	
	
}
#print $filt,"\t",$notfilt,"\n";
print "$count_canon\n";

sub exac {
	my ($in, $alt, $ref, $exac_cutoff, $syn,$in2) = @_;
	my %in_hash = %$in;
	my %tumor_supp = %$in2;
	my $check = 0;
	my $coding = 0;
	my $exac_filt = 1;
	my $somatic = 0;
	#print "NEW VARIANT___________________________\n\n\n\n";
	foreach my $key (sort keys %in_hash) {
		my $csq = $in_hash{$key}[0];
		my $exac_maf = $in_hash{$key}[46];
		my $canon = $in_hash{$key}[24];
		
		#print "$canon\t$in_hash{$key}[1]\t$exac_filt\t$csq\t$alt\t$exac_maf OUTERLOOP\n";
		if ($canon) {
			$somatic++;
			if ($in_hash{$key}[1] =~ /missense_variant/) {
				$coding++;
				$exac_filt = control_exac($exac_maf,$alt,$ref,$exac_cutoff,$csq);
				#print "$canon\t$in_hash{$key}[1]\t$exac_filt\t$csq\t$alt\t$exac_maf IFLOOP\n";
			}
			elsif ($in_hash{$key}[1] =~ /stop_gained/) {
				#print "$in_hash{$key}[3] STOP GAINED\n";
				if (defined $tumor_SUPP{$in_hash{$key}[3]}) {
					$exac_filt = 0;
					last;
				}
				$coding++;
				$exac_filt = control_exac($exac_maf,$alt,$ref,$exac_cutoff,$csq);
				#print "$canon\t$in_hash{$key}[1]\t$exac_filt\t$csq\t$alt\t$exac_maf IFLOOP\n";
			}
			elsif ($syn == 1) { 
				if ($in_hash{$key}[1] =~ /synonymous_variant/ ) {
					$coding++;
					$exac_filt = control_exac($exac_maf,$alt,$ref,$exac_cutoff,$csq);
				}
			}
			
			last;
		}
		else {
			
		}	
			
	}
	if ($coding >= 1 && $exac_filt == 1) { 
		$check++; 
	}
	return $check, $somatic;
}

sub vaf {
	my @sent = @_;
	my $tumorc = $sent[1];
	my $check = 0;
	if ($sent[0] =~/(\S+)\t(\S+)/) {
		my @array = split(':', $2);
		if ($array[2]/$tumorc > 0.05 && $array[2]/$tumorc < 0.4 ) {
			#print "$array[2]/$tumorc\n";
			
			if ($array[7] > 0 ) {
				$check++; 
			}
		}
		
		
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


















