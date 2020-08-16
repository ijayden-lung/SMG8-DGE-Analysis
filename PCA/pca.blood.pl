#!/usr/bin/perl -w
use Statistics::Descriptive;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::Random;


my @control_sample=qw/18DG0778 09DG00934 15DG2630 17DG0832 
14DG1171 18DG0120 19DG0060 17DG0977 16DG1333 16DG0662
16DG0744 18DG0147 15DG2234 16DG0328 10DG0840 15DG0371
16DG0991 19DG0555 18DG1094 18DG0734 16DG1068
14DG1686 16DG1048 19DG0591 17DG0969 17DG0444
18DG0638 10DG0792 18DG0717 15DG1349
15DG0678 16DG1186 18DG0646 15DG0918 
/;

my @test_sample = qw/19DG0152L 14DG1661L 19DG1391L 19DG1424L 19DG2599L/;
my %control;
foreach my $sample (@control_sample){
	$control{$sample} = '';
}
my %test;
foreach my $sample (@test_sample){
	$test{$sample} = '';
}

open FILE,"../tpm/blood_tpm.tsv";
my $header = <FILE>;
chomp $header;
my @header = split /\t/,$header;
open OUT2,">blood.input.tsv";
print OUT2 "gene_id";
foreach my $sample (@control_sample){
	my $sam = $sample;
	$sam =~ s/L//g;
	print OUT2 "\t $sam";
}
my @t_sample;
foreach my $sample (@test_sample){
	my $sam = $sample;
	$sam =~ s/L//g;
	push @t_sample,$sam;
	print OUT2 "\t $sam";
}
print OUT2 "\n";
while(<FILE>){
	chomp;
	my @data = split;
	my @control;
	my @test;
	my %print_control;
	my %print_test;
	my @all;
	for(my$i=6;$i<@header;$i++){
		$data[$i] += 1;
		if(exists $control{$header[$i]}){
			my $log = log($data[$i])/log(2);
			#push @control,$log;
			$print_control{$header[$i]} = $log;
			push @control,$data[$i];
			#$print_control{$header[$i]} = $data[$i];
			push @all,$data[$i];
		}
		elsif(exists $test{$header[$i]}){
			my $log = log($data[$i])/log(2);
			#push @test,$log;
			$print_test{$header[$i]} = $log;
			push @test,$data[$i];
			#$print_test{$header[$i]} = $data[$i];
			push @all,$data[$i];
		}
	}
	my $max = &max(@all);
	next if $max < 2;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@control);
	my $mean = $stat->mean();
	my $standard_deviation = $stat->standard_deviation();
	my @control_rmOutlier;
	if($standard_deviation == 0){
		@control_rmOutlier = @control;
	}
	else{
		for(my$j=0;$j<@control;$j++){
			my $zscore = ($control[$j]-$mean)/$standard_deviation;
			if(abs($zscore) < 2){
				push @control_rmOutlier,$control[$j];
			}
		}
	}
	my $control_sum = &sum(@control_rmOutlier);
	my $test_sum = &sum(@test);
	my $ave = ($control_sum+$test_sum)/(@control_rmOutlier+@test);
	next if $ave <= 1;
	print OUT2 "$data[4]";
	foreach my $sample (@control_sample){
		my $tpm = $print_control{$sample};
		print OUT2 "\t$tpm";
	}
	foreach my $sample (@test_sample){
		my $tpm = $print_test{$sample};
		print OUT2 "\t$tpm";
	}
	print OUT2 "\n";
}
my @group1 = ("SMG+")x34;
my @group2 = ("SMG8/9-")x5;
my @blank = ("")x34;
print OUT2 join ("\t","Group",@group1,@group2),"\n";
print OUT2 join ("\t","Label",@blank,@t_sample),"\n";
