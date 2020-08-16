#!/usr/bin/perl -w

use Statistics::Descriptive;


my @sample=qw/18DG0778 09DG00934 15DG2630 17DG0832 
14DG1171 18DG0120 19DG0060 17DG0977 16DG1333 16DG0662
16DG0744 18DG0147 15DG2234 16DG0328 10DG0840 15DG0371
16DG0991 19DG0555 18DG1094 18DG0734 16DG1068
14DG1686 16DG1048 19DG0591 17DG0969 17DG0444
18DG0638 10DG0792 18DG0717 15DG1349
15DG0678 16DG1186 18DG0646 15DG0918 
19DG0152L 14DG1661L 19DG1391L 19DG1424L 19DG2599L
/;
my %sample;
foreach my $sam (@sample){
	$sample{$sam} = '';
}


open TPM,"../tpm/blood_tpm.tsv";
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

open FILE,"../tpm/blood_tpm.tsv";
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

open OUT,">blood.tpm50.MiVsCV.txt";
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
