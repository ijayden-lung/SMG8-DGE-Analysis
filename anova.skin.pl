#!/usr/bin/perl -w
use Statistics::ANOVA 0.14;
use Statistics::Descriptive;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);


my ($input,$output) = @ARGV;

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


open FILE,"$input";
my $header = <FILE>;
chomp $header;
my @header = split /\t/,$header;
open OUT,">$output";
print OUT "gene_id\tmean\tlog2foldchange\tp_value";
foreach my $sample (@control_sample){
	my $sam = $sample;
	$sam =~ s/F//g;
	print OUT "\t $sam";
}
foreach my $sample (@test_sample){
	my $sam = $sample;
	$sam =~ s/F//g;
	print OUT "\t $sam";
}
print OUT "\n";
while(<FILE>){
	chomp;
	my @data = split;
	my @control;
	my @test;
	my $sum = 0;
	my %print_control;
	my %print_test;
	my $count = 0;
	my @all;
	my $hava = 0;
	for(my$i=6;$i<@header;$i++){
		$hava++ if $data[$i] ==0;
		$data[$i] += 1;
		if(exists $control{$header[$i]}){
			$count++;
			$sum += $data[$i] ;
			push @control,$data[$i];
			$print_control{$header[$i]} = $data[$i];
			push @all,$data[$i];
		}
		elsif(exists $test{$header[$i]}){
			$count++;
			$sum += $data[$i] ;
			push @test,$data[$i];
			$print_test{$header[$i]} = $data[$i];
			push @all,$data[$i];
		}
	}
	#next if $hava > 25;
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
	my $control_ave = $control_sum/@control_rmOutlier;
	my $test_ave = $test_sum/@test;
	my $ave = ($control_sum+$test_sum)/(@control_rmOutlier+@test);
	my $log2fc = log($test_ave/$control_ave)/log(2);
	next if $ave ==1;

	my $aov = Statistics::ANOVA->new();

	# Load the data:
	$aov->load_data({control => \@control_rmOutlier, test => \@test});

	my %res = $aov->anova(independent => 1, parametric => 1,ordinal => 0);
	print  OUT "$data[4]\t$ave\t$log2fc\t$res{'p_value'}";
	foreach my $sample (@control_sample){
		print OUT "\t$print_control{$sample}";
	}
	foreach my $sample (@test_sample){
		print OUT "\t$print_test{$sample}";
	}
	print OUT "\n";
}
my @group1 = ("SMG+")x26;
my @group2 = ("SMG8/9-")x3;
