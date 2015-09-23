#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $cosmicFile;

GetOptions (
	's=s' => \$cosmicFile  		# 
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($cosmicFile))
{
        print STDERR "Error - $0: cosmicFile file not specified\n";
		$p_fail=1;
}

if ($p_fail)
{
	printHelp();
	exit(255);
}

open CF, "<".$cosmicFile or die "INFO - $0: Error opening cosmicFile file. \n";

my %mut = ();
my %fus = ();
my %mutfus = ();
my %mut_dis = ();
my %fus_dis = ();
my %mutfus_dis = ();
my %sampleTOT = ();	
my %sampleTOT_dis = ();	

while(my $line=<CF>)
{
	chomp($line);

	my @fields = split('\t', $line);

	if($line =~ m/^Sample/)
	{
		next;
	}

	$sampleTOT{$fields[0]}++;
	$sampleTOT_dis{$fields[4]}{$fields[0]}++;

	if(defined($fields[5]))
	{
		my $sampleID = $fields[0];
		my $disease = $fields[4];
		my @genes= split("/", $fields[5]);
		my $genelen = scalar @genes;

		# check if it is a fusion
		if ($genelen > 1)
		{
			$mutfus_dis{$genes[0]}{$disease}{$sampleID}++;
			$mutfus_dis{$genes[1]}{$disease}{$sampleID}++;
			$fus_dis{$genes[0]}{$disease}{$sampleID}++;
			$fus_dis{$genes[1]}{$disease}{$sampleID}++;

			$mutfus{$genes[0]}{$sampleID}++;
			$mutfus{$genes[1]}{$sampleID}++;
			$fus{$genes[0]}{$sampleID}++;
			$fus{$genes[1]}{$sampleID}++;
		}
		elsif($genelen == 1)
		{
			$mutfus_dis{$genes[0]}{$disease}{$sampleID}++;
			$mut_dis{$genes[0]}{$disease}{$sampleID}++;

			$mutfus{$genes[0]}{$sampleID}++;
			$mut{$genes[0]}{$sampleID}++;
		}
		else
		{
			next;
		}	
	}
}

printf "GeneName\t#mut\t#fus\t#tot\t#freq_mut\t#freq_fus\t#samples\n";

my $samples = scalar keys %sampleTOT;

foreach my $g (keys %mutfus)
{
	my $tot_max = 0;
	my $mut_max = 0;
	my $fus_max = 0;
	my $tot = 0;
	my $mut = 0;
	my $fus = 0;

	$tot = scalar keys %{$mutfus{$g}};
	$mut = scalar keys %{$mut{$g}} if defined($mut{$g});
	$fus = scalar keys %{$fus{$g}} if defined($fus{$g});

	my $sample_dis = 0;
	my $freq_max = 0;
	my $name_disease = "";
	foreach my $d (keys %{$mutfus_dis{$g}})
	{
		my $tot_tmp = 0;
		$tot_tmp = scalar keys %{$mutfus_dis{$g}{$d}};
		my $mut_tmp = 0;
		$mut_tmp = scalar keys %{$mut_dis{$g}{$d}} if defined($mut_dis{$g}{$d});
		my $fus_tmp = 0;
		$fus_tmp = scalar keys %{$fus_dis{$g}{$d}} if defined($fus_dis{$g}{$d});
		my $sample_dis_tmp = 0;
		$sample_dis_tmp = scalar keys %{$sampleTOT_dis{$d}} if defined($sampleTOT_dis{$d});

		my $freq_tmp = $tot_tmp/$sample_dis_tmp*100;

		if($freq_max < $freq_tmp)
		{
			$name_disease = $d;
			$tot_max = $tot_tmp;			
		}
		if($freq_max < $freq_tmp)		
		{
			$mut_max = $mut_tmp;			
		}
		if($freq_max < $freq_tmp)
		{
			$fus_max = $fus_tmp;			
		}
		if($freq_max < $freq_tmp)
		{
			$sample_dis = $sample_dis_tmp;			
		}
	}

	#printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", $g, $mut, $fus, $tot, $samples, $mut_max, $fus_max, $tot_max, $sample_dis, $name_disease;
	my $stat_mut = 0.0;
	my $stat_fus = 0.0;
	$stat_mut = $mut / $samples * 100;
	$stat_fus = $fus / $samples * 100;
	printf "%s\t", $g;
	printf "%d\t", $mut;
	printf "%d\t", $fus;
	printf "%d\t", $tot;
	printf "%.5f\t", $stat_mut;
	printf "%.5f\t", $stat_fus;
	printf "%d\n", $samples;

}




# foreach my $g (keys %mut)
# {
# 	$genecount_mut{$g} = scalar keys %{$mut{$g}};
# }
# foreach my $g (keys %fus)
# {
# 	$genecount_fus{$g} = scalar keys %{$fus{$g}};
# }


close(CF);

sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
}










