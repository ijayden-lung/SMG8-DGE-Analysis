#!/usr/bin/perl -w
#

use Statistics::Descriptive;
use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
use Statistics::Multtest qw(:all);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
my ($TYPE) = @ARGV;


my %brain_gene;
open FILE,"table/tableS3.csv";
<FILE>;
while(<FILE>){
	chomp;
	my ($gene,$type) = (split /,/)[0,-1];
	$brain_gene{$gene} = $type if $type =~ "Stab";
}

my %sub;
open FILE,"table/core_NMD_substrate.txt";
while(<FILE>){
	chomp;
	$sub{$_} = '';
}


my %high_confidence;
open FILE,"table/tableS5.csv";
<FILE>;
while(<FILE>){
	chomp;
	my ($symbol,$type) = (split /,/)[0,-1];
	if($type =~ /siUPF1-Up/){
		$high_confidence{$symbol} = '';
	}
}



my %id2gene;
open FILE,"table/id2symbol.tsv";
while(<FILE>){
	chomp;
	my ($gene_id,$symbol) = split;
	$id2gene{$gene_id} = $symbol;
}


my %pval;
open FILE,"Anova.$TYPE.tsv";
<FILE>;
while(<FILE>){
	chomp;
	my ($gene_id,$mean_val,$log2FoldChange,$pvalue) = (split)[0,1,2,3];
	my ($id) = split /\./,$gene_id;
	$pval{$id} = $pvalue;
}
print "Start calculate qvalue\n";

my $p_ref = \%pval;

my $res;
#eval '$res = qvalue($p_ref)';
$res = BH($p_ref);
print "finish calcualting qvalue\n";

open FILE,"sort -k 4,4n Anova.$TYPE.tsv |";
my $header = <FILE>;
chomp $header;
my @header = split /\t/,$header;
open OUT,">DE.Anova.$TYPE.tsv";
open OUT2,">zscore.$TYPE.tsv";
print OUT "gene_id\tsymbol\tmean_tpm\tlog2FoldChange\tpvalue\tqvalue\tSig\ttableS3\ttableS5\tcoreNMDsubstrates\n";
print OUT2 "gene_id\t",join("\t",@header[4..$#header]),"\n";
my $total = 0;
my $up = 0;
my $down = 0;
my $substrate = 0;
my $sub_up  = 0;
my $tables5 = 0;
my $tables3 = 0;
my $tables3_up = 0;


while(<FILE>){
	chomp;
	my ($gene_id,$mean_val,$log2FoldChange,$pvalue) = (split)[0,1,2,3];
	my ($id) = split /\./,$gene_id;
	my $symbol =  "NA";
	if($id2gene{$id}){
		$symbol = $id2gene{$id};
	}
	else{
		print "$id\n";
		next;
	}
	my @data = split;
	my @tpm = @data[4..$#data];
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@tpm);
	my $mean = $stat->mean();

	my $standard_deviation=$stat->standard_deviation();
	if($standard_deviation == 0){
		print "@tpm\n";
	}
	print OUT2 "$id2gene{$id}";
	for(my$j=0;$j<@tpm;$j++){
		my $zscore = ($tpm[$j]-$mean)/$standard_deviation;
		print OUT2 "\t$zscore";
	}
	print OUT2 "\n";
	
	my $brain_related = "No";
	my $tableS5 = "No";
	my $coreNMD = "No";
	if(exists $brain_gene{$symbol}){
		$tables3++;
		$brain_related = $brain_gene{$symbol};
	}
	if(exists $sub{$symbol}){
		$substrate++;
		$coreNMD = "coreNMDsubstrate";
	}
	if(exists $high_confidence{$symbol}){
		$tables5++;
		$tableS5 = "tableS5";
	}
	my $sig = "No";
	if($pvalue ne "NA" && $res->{$id} < 0.05){
		if($log2FoldChange > 0.5){
			$sig = "Up";
			$up++;
			if(exists $brain_gene{$symbol}){
				$tables3_up++;
			}
			if(exists $sub{$symbol}){
				$sub_up++;
			}
			if(exists $high_confidence{$symbol}){
				$tables5_up++;
			}
		}
		elsif($log2FoldChange < -0.5){
			$sig = "Down";
			$down++;
		}
	}
	$total++;
	print OUT "$gene_id\t$symbol\t$mean_val\t$log2FoldChange\t$pvalue\t$res->{$id}\t$sig\t$brain_related\t$tableS5\t$coreNMD\n";
}

print "total gene: $total, up: $up, down: $down\n";
print "tables3: $tables3_up:$tables3\n";
print "coresub: $sub_up: $substrate\n";
print "tables5: $tables5_up:$tables5\n";
