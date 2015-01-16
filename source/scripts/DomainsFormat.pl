#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $goodDomains;
my $label;
my $domainfile;

GetOptions (
	'r=s' => \$goodDomains,  	# 
	'l=s' => \$label,  	# 
	's=s' => \$domainfile  		# 
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($goodDomains))
{
        print STDERR "Error - $0: goodDomains not specified\n";
		$p_fail=1;
}

if(!defined($label))
{
        print STDERR "Error - $0: label not specified\n";
		$p_fail=1;
}

if(!defined($domainfile))
{
        print STDERR "Error - $0: domainfile file not specified\n";
		$p_fail=1;
}

if ($p_fail)
{
	printHelp();
	exit(255);
}



open GOOD, "<".$goodDomains or die "INFO - $0: Error opening fusionList file. \n";
open DOMAINS, "<".$domainfile or die "INFO - $0: Error opening cosmicFile file. \n";

my %goodDomains = ();
my %goodRegions = ();
my %goodRepeats = ();
my %goodTD = ();
my %domains = ();
my %regions = ();
my %repeats = ();
my %topolDomain = ();



while(my $line=<GOOD>)
{
	chomp($line);

	my @fields = split('\t', $line);

	if($line =~ m/^Domain/)
	{
		$goodDomains{$fields[1]}++;
	}
	elsif($line =~ m/^Region/)
	{
		$goodRegions{$fields[1]}++;
	}
	elsif($line =~ m/^Repeat/)
	{
		$goodRepeats{$fields[1]}++;
	}
	elsif($line =~ m/^Topological domain/)
	{
		$goodTD{$fields[1]}++;
	}
	else
	{
		next;
	}
}

#printf "orig\t";

foreach my $k (sort keys %goodDomains)
{
	printf "%s_%s\t", $k, $label;
}
foreach my $k (sort keys %goodRegions)
{
	printf "%s_%s\t", $k, $label;
}
foreach my $k (sort keys %goodRepeats)
{
	printf "%s_%s\t", $k, $label;
}
foreach my $k (sort keys %goodTD)
{
	printf "%s_%s\t", $k, $label;
}

printf "\n";


while(my $line=<DOMAINS>)
{
	chomp($line);

	my @fields = split(',', $line);

	if($line =~ m/^Conserved_Domain/)
	{
		next;
	}
	elsif($line =~ m/^Lost_Domain/)
	{
		next;
	}

	%domains = ();
	%regions = ();
	%repeats = ();
	%topolDomain = ();


	foreach my $f (@fields)
	{

		if($f =~ m/Domain=/)
		{
			$f =~ s/Domain=//;
			foreach my $d (sort keys %goodDomains)
			{
				if ($f =~ m/${d}/)
				{
					$domains{$d}++;
				}
			}
		}
		elsif($f =~ m/Region=/)
		{
			$f =~ s/Region=//;
			foreach my $d (sort keys %goodRegions)
			{
				if ($f =~ m/${d}/)
				{
					$regions{$d}++;
				}
			}
		}
		elsif($f =~ m/Repeat=/)
		{
			$f =~ s/Repeat=//;
			foreach my $d (sort keys %goodRepeats)
			{
				if ($f =~ m/${d}/)
				{
					$repeats{$d}++;
				}
			}
		}
		elsif($f =~ m/Topological domain=/)
		{
			$f =~ s/Topological domain=//;
			foreach my $d (sort keys %goodTD)
			{
				if ($f =~ m/${d}/)
				{
					$topolDomain{$d}++;
				}
			}
		}
	}

	foreach my $d (sort keys %goodDomains)
	{
		if(defined($domains{$d}))
		{
			printf "1\t";
		}
		else
		{
			printf "0\t";
		}
	}
	foreach my $d (sort keys %goodRegions)
	{
		if(defined($regions{$d}))
		{
			printf "1\t";
		}
		else
		{
			printf "0\t";
		}
	}
	foreach my $d (sort keys %goodRepeats)
	{
		if(defined($repeats{$d}))
		{
			printf "1\t";
		}
		else
		{
			printf "0\t";
		}
	}
	foreach my $d (sort keys %goodTD)
	{
		if(defined($topolDomain{$d}))
		{
			printf "1\t";
		}
		else
		{
			printf "0\t";
		}
	}

	printf "\n";

}

close(GOOD);
close(DOMAINS);





sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
}


