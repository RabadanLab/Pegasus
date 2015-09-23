#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $fusionList;
my $cosmicFile;

GetOptions (
	'r=s' => \$fusionList,  	# 
	's=s' => \$cosmicFile  		# it has to be the cosmic table created with makeCosmicTable.pl
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($fusionList))
{
        print STDERR "Error - $0: fusionList not specified\n";
		$p_fail=1;
}

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



open FL, "<".$fusionList or die "INFO - $0: Error opening fusionList file. \n";
open CF, "<".$cosmicFile or die "INFO - $0: Error opening cosmicFile file. \n";

my %mut = ();
my %fus = ();
	
while(my $line=<CF>)
{
	chomp($line);

	my @fields = split('\t', $line);

	if($line =~ m/^Gene/)
	{
		next;
	}

	my $mut_freq = $fields[4];
	my $fus_freq = $fields[5];
	my $g_name = $fields[0];

	$mut{$g_name} = $mut_freq;
	$fus{$g_name} = $fus_freq;

}

printf "GeneName1\tGeneName2\tGene1_CosmicMutations\tGene2_CosmicMutations\tGene1_CosmicFusions\tGene2_CosmicFusions\n";

while(my $line=<FL>)
{
	chomp($line);

	my @fields = split('\t', $line);

	if($line =~ m/^Gene_Name1/)
	{
		next;
	}

	my @genelist1 = split (',', $fields[0]);
	my @genelist2 = split (',', $fields[1]);

	printf "%s\t%s\t", $fields[0], $fields[1];
	
	#check in cosmic mutation gene1
	my $max = -100;
	for my $g (@genelist1)
	{
		if(defined($mut{$g}))
		{
			if ($mut{$g} > $max)
			{
				$max = $mut{$g};
			}
		}
	}
	if($max < 0)
	{
		$max = 0;
	}
	printf "%f", $max;
	printf "\t";

	#check in cosmic mutation gene2
	$max = 0;
	for my $g (@genelist2)
	{
		if(defined($mut{$g}))
		{
			if ($mut{$g} > $max)
			{
				$max = $mut{$g};
			}
		}
	}
	if($max < 0)
	{
		$max = 0;
	}
	printf "%f", $max;
	printf "\t";

	#check in cosmic mutation gene1
	$max = 0;
	for my $g (@genelist1)
	{
		if(defined($fus{$g}))
		{
			if ($fus{$g} > $max)
			{
				$max = $fus{$g};
			}
		}
	}
	if($max < 0)
	{
		$max = 0;
	}
	printf "%f", $max;
	printf "\t";

	#check in cosmic mutation gene2
	$max = 0;
	for my $g (@genelist2)
	{
		if(defined($fus{$g}))
		{
			if ($fus{$g} > $max)
			{
				$max = $fus{$g};
			}
		}
	}
	if($max < 0)
	{
		$max = 0;
	}
	printf "%f", $max;
	printf "\n";




}

close(CF);
close(FL);





sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
}


