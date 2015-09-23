#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;

my $DomainFile;
my $OverlappingFile;
my $label;

GetOptions (
	'r=s' => \$DomainFile,  	# List of Domains from Pegasus
	'l=s' => \$label,  			# List of Domains from Pegasus
	's=s' => \$OverlappingFile  # List of relevant domains to be overlapped
);  

# 1 is true, 0 false
my $p_fail;

if(!defined($DomainFile))
{
        print STDERR "Error - $0: DomainFile not specified\n";
		$p_fail=1;
}

if(!defined($label))
{
        print STDERR "Error - $0: label not specified\n";
		$p_fail=1;
}

if(!defined($OverlappingFile))
{
        print STDERR "Error - $0: OverlappingFile file not specified\n";
		$p_fail=1;
}

if ($p_fail)
{
	printHelp();
	exit(255);
}



open DF, "<".$DomainFile or die "INFO - $0: Error opening DomainFile file. \n";
open OF, "<".$OverlappingFile or die "INFO - $0: Error opening OverlappingFile file. \n";

my %targetDomains = ();	# domains to be overlapped to
	
while(my $line=<OF>)
{
	chomp($line);

	my @f = split('\t', $line);

	if($line =~ m/^Gene/)
	{
		next;
	}

	$targetDomains{$f[1]}++;

}


printf "%s_MatchedDomain\t%s_MatchedRegion\t%s_MatchedRepeats\t%s_MatchedTopolDomain\n", $label, $label, $label, $label;
<DF>;
while(my $line=<DF>)
{
	chomp($line);

	my @domain_list = split(',', $line);

	if($line =~ m/^Gene_Name1/)
	{
		next;
	}

	my $num_matchingdomain = 0;
	my $num_matchingregion = 0;
	my $num_matchingrepeat = 0;
	my $num_matchingtopoldomain = 0;
	foreach my $domain (@domain_list)
	{
		if($domain =~ m/^Domain=/)
		{
			$domain =~ s/Domain=//;
			$domain =~ tr/A-Z/a-z/;

			foreach my $d (keys %targetDomains)
			{
				#printf "%s\t%s\n", $domain, $d;
				if($d =~ m/$domain/)
				{
					$num_matchingdomain=1;
				}
			}
		}
		if($domain =~ m/^Repeat=/)
		{
			$domain =~ s/Repeat=//;
			$domain =~ tr/A-Z/a-z/;

			foreach my $d (keys %targetDomains)
			{
				#printf "%s\t%s\n", $domain, $d;
				if($d =~ m/$domain/)
				{
					$num_matchingrepeat=1;
				}
			}
		}
		if($domain =~ m/^Region=/)
		{
			$domain =~ s/Region=//;
			$domain =~ tr/A-Z/a-z/;

			foreach my $d (keys %targetDomains)
			{
				#printf "%s\t%s\n", $domain, $d;
				if($d =~ m/$domain/)
				{
					$num_matchingregion=1;
				}
			}
		}
		if($domain =~ m/^Topological domain=/)
		{
			$domain =~ s/Topological domain=//;
			$domain =~ tr/A-Z/a-z/;

			foreach my $d (keys %targetDomains)
			{
				#printf "%s\t%s\n", $domain, $d;
				if($d =~ m/$domain/)
				{
					$num_matchingtopoldomain=1;
				}
			}
		}
	}
	printf "%d\t", $num_matchingdomain;
	printf "%d\t", $num_matchingregion;
	printf "%d\t", $num_matchingrepeat;
	printf "%d\n", $num_matchingtopoldomain;
}

close(DF);
close(OF);





sub printHelp
{
	print "\t-i [input file]\n";
	print "\t-o [output file]\n";
}


