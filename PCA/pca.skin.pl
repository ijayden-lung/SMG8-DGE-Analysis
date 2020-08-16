#!/usr/bin/perl -w
use Statistics::Descriptive;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::Random;


my @control_sample = qw/10DG0934    11DG0060    11DG0268
11DG0165	11DG0840	12DG1794	14DG2098	15DG2154	15DG2530	16DG0559
16DG0676	16DG0790	16DG1353	18DG0180	18DG0295	19DG0151
13DG2283	14DG2019	16DG0144	16DG0518	16DG0932	17DG0349
18DG0348    18DG0464F   18DG0603F   19DG0230
/;

my @test_sample = qw/19DG0152F 19DG1391F 19DG2599F/;
my %control;
foreach my $sample (@control_sample){
	$control{$sample} = '';
}
my %test;
foreach my $sample (@test_sample){
	$test{$sample} = '';
}


open FILE,"../tpm/skin_tpm.tsv";
my $header = <FILE>;
chomp $header;
my @header = split /\t/,$header;
open OUT2,">skin.input.tsv";
print OUT2 "gene_id";
foreach my $sample (@control_sample){
	my $sam = $sample;
	$sam =~ s/F//g;
	print OUT2 "\t $sam";
}
my @t_sample;
foreach my $sample (@test_sample){
	my $sam = $sample;
	$sam =~ s/F//g;
	push @t_sample, $sam;
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
	print  OUT2 "$data[4]";
	foreach my $sample (@control_sample){
		print OUT2 "\t$print_control{$sample}";
	}
	foreach my $sample (@test_sample){
		print OUT2 "\t$print_test{$sample}";
	}
	print OUT2 "\n";
}
my @group1 = ("SMG+")x26;
my @group2 = ("SMG8/9-")x3;

my @blank = ("")x26;
print OUT2 join ("\t","Group",@group1,@group2),"\n";
print OUT2 join ("\t","Label",@blank,@t_sample),"\n";
