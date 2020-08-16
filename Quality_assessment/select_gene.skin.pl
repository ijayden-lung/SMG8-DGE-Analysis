#!/usr/bin/perl -w

use Statistics::Descriptive;

my @sample = qw/10DG0934    11DG0060    11DG0268
11DG0165	11DG0840	12DG1794	14DG2098	15DG2154	15DG2530	16DG0559
16DG0676	16DG0790	16DG1353	18DG0180	18DG0295	19DG0151
13DG2283	14DG2019	16DG0144	16DG0518	16DG0932	17DG0349
18DG0348    18DG0464F   18DG0603F   19DG0230
19DG0152F 19DG1391F 19DG2599F
/;

my %sample;
foreach my $sam (@sample){
	$sample{$sam} = '';
}


open TPM,"../tpm/skin_tpm.tsv";
my $header = <TPM>;
chomp $header;
my @header = split /\t/,$header;
my %tpm;
while(<TPM>){
	chomp;
	my @data = split;
	my @tpm;
	for(my$i=6;$i<@data;$i++){
		if(exists $sample{$header[$i]}){
			push @tpm,$data[$i];
		}
	}
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@tpm);
		my $mean = $stat->mean();
		next if $mean < 50;
		my $min = $stat->min();
		next if $min == 0;
		$tpm{$data[5]} = \@tpm;
}

open FILE,"../tpm/skin_tpm.tsv";
$header = <FILE>;
chomp $header;
@header = split /\t/,$header;
my %cv;
while(<FILE>){
	chomp;
	my @data = split;
	my @tpm;
	for(my$i=6;$i<@data;$i++){
		if(exists $sample{$header[$i]}){
			push @tpm,$data[$i];
		}
	}
	if(exists $tpm{$data[5]}){
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@tpm);
		my $mean = $stat->mean();
		my $standard_deviation = $stat->standard_deviation();
		my $cv = $standard_deviation/$mean;
		$cv{$data[5]} = $cv;
	}
}

my $gene_num = keys %tpm;

open OUT,">skin.tpm50.MiVsCV.txt";
print OUT "gene\tM_i\tCV\n";
foreach my $gene1 (keys %tpm){
	my $sum = 0;
	foreach my $gene2 (keys %tpm){
		next if $gene2 eq $gene1;
		my $val1 = $tpm{$gene1}; #Ri
		my $val2 = $tpm{$gene2}; #Rj
		my @Aij;
		for (my $k=0;$k<@$val1;$k++){
			my $Rijk = ($val1->[$k])/($val2->[$k]); #Rijk
			push @Aij, log($Rijk)/log(2);

		}
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@Aij);
		my $Vij = $stat->standard_deviation();
		$sum += $Vij;
	}
	my $Mi = $sum/($gene_num-1);
	print OUT "$gene1\t$Mi\t$cv{$gene1}\n";
}
