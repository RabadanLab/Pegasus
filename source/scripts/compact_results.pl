#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my %GeneTranscriptID1_hash;
my %GeneTranscriptID2_hash;
my %FrameType_hash;
my %FusedSequence_hash;
my %ProteinStart1_hash;
my %ProteinStop1_hash;
my %ProteinStart2_hash;
my %ProteinStop2_hash;
my %ProteinSequence_hash;
my %ExonBreak1_hash;
my %ExonBreak2_hash;
my %BpRegion1_hash;
my %BpRegion2_hash;
my %DomainConserved1_hash;
my %DomainLost1_hash;
my %DomainConserved2_hash;
my %DomainLost2_hash;

my $GeneName1="";
my $GeneName2="";;
my $GeneBreakpoin1="";;
my $GeneBreakpoin2="";;

while(my $line = <>)
{
	chomp($line);
	my @fields = split(/\t/, $line);

	$GeneName1 = $fields[0];
	$GeneName2 = $fields[1];
	$GeneBreakpoin1 = $fields[2];
	$GeneBreakpoin2 = $fields[3];
	my $GeneTranscriptID1 = $fields[4];
	my $GeneTranscriptID2 = $fields[5];
	my $FrameType = $fields[6];
	my $FusedSequence = $fields[7];
	my $ProteinStart1 = $fields[8];
	my $ProteinStop1 = $fields[9];
	my $ProteinStart2 = $fields[10];
	my $ProteinStop2 = $fields[11];
	my $ProteinSequence = $fields[12];
	my $ExonBreak1 = $fields[13];
	my $ExonBreak2 = $fields[14];
	my $BpRegion1 = $fields[15];
	my $BpRegion2 = $fields[16];
	my $DomainConserved1 = $fields[17];
	my $DomainLost1 = $fields[18];
	my $DomainConserved2 = $fields[19];
	my $DomainLost2 = $fields[20];
	
	(!defined($GeneTranscriptID1_hash{$GeneTranscriptID1})) ? $GeneTranscriptID1_hash{$GeneTranscriptID1}=1 : $GeneTranscriptID1_hash{$GeneTranscriptID1}++; 
	
	(!defined($GeneTranscriptID2_hash{$GeneTranscriptID2})) ? $GeneTranscriptID2_hash{$GeneTranscriptID2}=1 : $GeneTranscriptID2_hash{$GeneTranscriptID2}++; 

	$FrameType_hash{$FrameType}{"$GeneTranscriptID1-$GeneTranscriptID2"}++;
	$FusedSequence_hash{$FusedSequence}{"$GeneTranscriptID1-$GeneTranscriptID2"}++;
	$ProteinStart1_hash{$ProteinStart1}{"$GeneTranscriptID1"}++;
	$ProteinStop1_hash{$ProteinStop1}{"$GeneTranscriptID1"}++;
	$ProteinStart2_hash{$ProteinStart2}{"$GeneTranscriptID2"}++;
	$ProteinStop2_hash{$ProteinStop2}{"$GeneTranscriptID2"}++;
	$ProteinSequence_hash{$ProteinSequence}{"$GeneTranscriptID1-$GeneTranscriptID2"}++;
	$ExonBreak1_hash{$ExonBreak1}{"$GeneTranscriptID1"}++;
	$ExonBreak2_hash{$ExonBreak2}{"$GeneTranscriptID2"}++;
	$BpRegion1_hash{$BpRegion1}{"$GeneTranscriptID1"}++;
	$BpRegion2_hash{$BpRegion2}{"$GeneTranscriptID2"}++;
	$DomainConserved1_hash{$DomainConserved1}{"$GeneTranscriptID1"}++;
	$DomainLost1_hash{$DomainLost1}{"$GeneTranscriptID1"}++;
	$DomainConserved2_hash{$DomainConserved2}{"$GeneTranscriptID2"}++;
	$DomainLost2_hash{$DomainLost2}{"$GeneTranscriptID2"}++;
}

printf "%s\t%s\t%s\t%s\t", $GeneName1, $GeneName2, $GeneBreakpoin1, $GeneBreakpoin2;

foreach my $GeneTranscriptID1 (keys %GeneTranscriptID1_hash)
{
	printf "%s,", $GeneTranscriptID1;
}

printf "\t";

foreach my $GeneTranscriptID2 (keys %GeneTranscriptID2_hash)
{
	printf "%s,", $GeneTranscriptID2;
}

printf "\t";

foreach my $FrameType (keys %FrameType_hash)
{
	printf "%s", $FrameType;
	#foreach my $GeneTranscriptID (keys %{$FrameType_hash{$FrameType}})
	#{
	#	printf "%s,", $GeneTranscriptID;
	#}
	printf " ";
}

printf "\t";

foreach my $FusedSequence (keys %FusedSequence_hash)
{
	printf "%s:", $FusedSequence;
	foreach my $GeneTranscriptID (keys %{$FusedSequence_hash{$FusedSequence}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ProteinStart1 (keys %ProteinStart1_hash)
{
	printf "%s:", $ProteinStart1;
	foreach my $GeneTranscriptID (keys %{$ProteinStart1_hash{$ProteinStart1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ProteinStop1 (keys %ProteinStop1_hash)
{
	printf "%s:", $ProteinStop1;
	foreach my $GeneTranscriptID (keys %{$ProteinStop1_hash{$ProteinStop1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ProteinStart2 (keys %ProteinStart2_hash)
{
	printf "%s:", $ProteinStart2;
	foreach my $GeneTranscriptID (keys %{$ProteinStart2_hash{$ProteinStart2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ProteinStop2 (keys %ProteinStop2_hash)
{
	printf "%s:", $ProteinStop2;
	foreach my $GeneTranscriptID (keys %{$ProteinStop2_hash{$ProteinStop2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ProteinSequence (keys %ProteinSequence_hash)
{
	printf "%s:", $ProteinSequence;
	foreach my $GeneTranscriptID (keys %{$ProteinSequence_hash{$ProteinSequence}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ExonBreak1 (keys %ExonBreak1_hash)
{
	printf "%s:", $ExonBreak1;
	foreach my $GeneTranscriptID (keys %{$ExonBreak1_hash{$ExonBreak1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $ExonBreak2 (keys %ExonBreak2_hash)
{
	printf "%s:", $ExonBreak2;
	foreach my $GeneTranscriptID (keys %{$ExonBreak2_hash{$ExonBreak2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $BpRegion1 (keys %BpRegion1_hash)
{
	printf "%s:", $BpRegion1;
	foreach my $GeneTranscriptID (keys %{$BpRegion1_hash{$BpRegion1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $BpRegion2 (keys %BpRegion2_hash)
{
	printf "%s:", $BpRegion2;
	foreach my $GeneTranscriptID (keys %{$BpRegion2_hash{$BpRegion2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $DomainConserved1 (sort keys %DomainConserved1_hash)
{
	printf "%s:", $DomainConserved1;
	foreach my $GeneTranscriptID (keys %{$DomainConserved1_hash{$DomainConserved1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $DomainLost1 (sort keys %DomainLost1_hash)
{
	printf "%s:", $DomainLost1;
	foreach my $GeneTranscriptID (keys %{$DomainLost1_hash{$DomainLost1}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $DomainConserved2 (sort keys %DomainConserved2_hash)
{
	printf "%s:", $DomainConserved2;
	foreach my $GeneTranscriptID (keys %{$DomainConserved2_hash{$DomainConserved2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\t";

foreach my $DomainLost2 (sort keys %DomainLost2_hash)
{
	printf "%s:", $DomainLost2;
	foreach my $GeneTranscriptID (keys %{$DomainLost2_hash{$DomainLost2}})
	{
		printf "%s,", $GeneTranscriptID;
	}
	printf " ";
}

printf "\n";



















