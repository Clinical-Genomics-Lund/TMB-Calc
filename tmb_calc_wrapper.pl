#!/usr/bin/perl -w
use strict;


my $samples = "test_runs_output/sent_free_samples_tumorconc";
my $coding_size = $ARGV[0];
my $outpath = $ARGV[1];
open(LIST, $samples) or die $!;
my $sent;
my $free;
my $tc;
my $id;
my $c = 0;
while ( <LIST> ) {
    my @input = split /:/;
    $id = $input[0];
    $sent = $input[1];
    $free = $input[2];
    $tc = $input[3];
    chomp $tc;

    my $sent_cmd = "./filter.pl -v $sent -s -t $tc -d 200 -a sent -e 0.0001 -o $outpath/$id.sent.bed";
    my $free_cmd = "./filter.pl -v $free -s -t $tc -d 200 -a free -e 0.0001 -o $outpath/$id.free.bed";
    
    print STDOUT "RUNNING SENTIEON VCF\n";
    my $cmd_sent = `$sent_cmd`;

    print STDOUT "RUNNING FREEBAYES VCF\n";
    my $cmd_free = `$free_cmd`;
    


    my $SGZ_free = "$outpath/free.mu.basic.sgz.txt";
    my $SGZ_sent = "$outpath/sent.mu.basic.sgz.txt";

    my $est_free = sgz_compare($SGZ_free);
    my $est_sent = sgz_compare($SGZ_sent);
    my %consensus;
    my $count = 0;
    open (INT, '>', "$outpath/$id.intersect.bed") or die $!;
    foreach my $key (keys %$est_sent) {
        if (defined $$est_free{$key}) {
           if ($$est_free{$key} eq "somatic" || $$est_sent{$key} eq "somatic") {
               $consensus{$key} = "somatic";
               my @a = split /:/, $key;
               print INT "$a[0]\t$a[1]\t$a[1]\n";
               $count++;
           }
           elsif($$est_free{$key} eq "germline" && $$est_sent{$key} ne "somatic") {
               $consensus{$key} = "germline";
           }
           elsif($$est_free{$key} ne "somatic" && $$est_sent{$key} eq "germline") {
               $consensus{$key} = "germline";
           }
           else {
               $consensus{$key} = "somatic";
               my @a = split /:/, $key;
               print INT "$a[0]\t$a[1]\t$a[1]\n";
               $count++;
           }
        }
    }


    my $TMB = $count / $coding_size;
    print "$id TMB: ",$TMB,"\n";

    $c++;
    #if ($c >= 1) { exit;}
}
close LIST;



sub sgz_compare {
    my ($sgz) = shift;

    open (SGZ, $sgz) or die $!;
    my %est_som_germ;
    while ( <SGZ> ) {
        my @line = split /\t/;
        if ($line[1] =~ /chr/) {
            $line[1] =~ s/chr//;
            chomp $line[4];
            $est_som_germ{$line[1]} = $line[4];
        }
    }
    close SGZ;
    return \%est_som_germ;

}