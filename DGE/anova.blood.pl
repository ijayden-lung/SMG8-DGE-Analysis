#!/usr/bin/perl -w
use Statistics::ANOVA 0.14;
use Statistics::Descriptive;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#

#Modified May 15
#




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
open OUT,">Anova.blood.tsv";
print OUT "gene_id\tmean\tlog2foldchange\tp_value";
foreach my $sample (@control_sample){
	print OUT "\t$sample";
}
foreach my $sample (@test_sample){
	$sam = $sample;
	$sam =~ s/L//g;
	print OUT "\t$sam";
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
	for(my$i=6;$i<@header;$i++){
		$data[$i] += 1;
		if(exists $control{$header[$i]}){
			$count++;
			$sum += $data[$i] ;
			push @control,$data[$i];
			$print_control{$header[$i]} = $data[$i];
			push @all,$data[$i];
		}
		elsif(exists $test{$header[$i]}){
			$sum += $data[$i] ;
			$count++;
			push @test,$data[$i];
			$print_test{$header[$i]} = $data[$i];
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
	my $control_ave = $control_sum/@control_rmOutlier;
	my $test_ave = $test_sum/@test;
	my $ave = ($control_sum+$test_sum)/(@control_rmOutlier+@test);
	my $log2fc = log($test_ave/$control_ave)/log(2);
	next if $ave == 1;
	my $aov = Statistics::ANOVA->new();
	# Load the data:
	$aov->load_data({control => \@control_rmOutlier, test => \@test});

	my %res = $aov->anova(independent => 1, parametric => 1,ordinal => 0);
	while(my($key,$val) = each %res){
		#print "$key\t$val\n";
	}

	print  OUT "$data[4]\t$ave\t$log2fc\t$res{'p_value'}";
	foreach my $sample (@control_sample){
		print OUT "\t$print_control{$sample}";
	}
	foreach my $sample (@test_sample){
		print OUT "\t$print_test{$sample}";
	}
	print OUT "\n";
}


