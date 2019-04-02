#!/usr/bin/perl -w
use strict;
use CMD::vcf2;
use Data::Dumper;
## Script to change name of family for scout upload. Helpfull if a sample has been rerun without renaming family in first step
## Takes the FULL PATH!! of SNV-VCF from run as first argument and the new family-ID as the second
## From its location the script finds relevant files and changes family id in these, VCF, PED, YAML

my $vcf = CMD::vcf2->new('file'=> $ARGV[0] );

my $c = 0;
while ( my $a = $vcf->next_var() ) {
$c++;
print Dumper($a);
if ($c > 0) {exit;}
}