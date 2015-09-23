#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $pfam;

GetOptions (
	's=s' => \$pfam  		# 
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($pfam))
{
        print STDERR "Error - $0: pfam file not specified\n";
		$p_fail=1;
}

if ($p_fail)
{
	printHelp();
	exit(255);
}

my %inter_pfam = ();

open INTERACTION, "<".$pfam or die "INFO - $0: Error opening pfam file. \n";

while(my $line=<INTERACTION>)
{
	chomp($line);

	my @fields = split('\t', $line);

	$inter_pfam{$fields[0]}{$fields[1]}++;
	$inter_pfam{$fields[1]}{$fields[0]}++;
}


while(my $line=<STDIN>)
{
	chomp($line);

	my @fields = split('\t', $line);
	
	if(defined($fields[12]) && defined($line))
	{
		if(defined($inter_pfam{$fields[12]}))
		{
			foreach my $k (keys %{$inter_pfam{$fields[12]}})
			{
				printf $line."\t";
				printf $k."\n";
			}
		}
		else
		{
			printf $line."\t"."\n";
		}
	}
	else
	{
		printf $line."\t"."\n";
	}
}

close(INTERACTION);


sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
}


