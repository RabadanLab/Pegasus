#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $chimerascan_report;
my $fus_annot_report;

GetOptions (
	'i=s' => \$chimerascan_report,
	'r=s' => \$fus_annot_report
);  

if(!defined($chimerascan_report))
{
        print "Error - Input file not specified - Exit......\n";
        exit 1;
}

if(!defined($fus_annot_report))
{
        print "Error - Input file not specified - Exit......\n";
        exit 1;
}


my %fus_annotations;

open IN, "<".$fus_annot_report or die $!;

while(my $line = <IN>)
{
	#printf "$line";
	chomp($line);
	my @fields = split(/\t/, $line);

	if($line =~ m/^$/ || $line =~ m/^\t*$/)
	{
		next;
	}

	my $GeneName1 = $fields[0];
	my $GeneName2 = $fields[1];
	my $GeneBreakpoin1 = $fields[2];
	my $GeneBreakpoin2 = $fields[3];
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
	my $DomainConserved1 = $fields[15];
	my $DomainLost1 = $fields[16];
	my $DomainConserved2 = $fields[17];
	my $DomainLost2 = $fields[18];

	if(!defined($fus_annotations{$GeneName1}{$GeneName2}{$GeneBreakpoin1}{$GeneBreakpoin2}))
	{
		#$fus_annotations{$GeneName1}{$GeneName2}{$GeneBreakpoin1}{$GeneBreakpoin2} = ( $GeneTranscriptID1, $GeneTranscriptID2, $FusedSequence, $ProteinStart1, $ProteinStop1, $ProteinStart2, $ProteinStop2, $ProteinSequence, $ExonBreak1, $ExonBreak2, $DomainConserved1, $DomainLost1, $DomainConserved2, $DomainLost2);
		$fus_annotations{$GeneName1}{$GeneName2}{$GeneBreakpoin1}{$GeneBreakpoin2} = $line;

	}
	#else
	#{
	#	printf "Double match for %s %s %s %s\n", $GeneName1, $GeneName2, $GeneBreakpoin1, $GeneBreakpoin2;
	#}
}

close(IN);
#sleep 50;
#foreach my $GeneName1 (keys %fus_annotations)
#{
#	foreach my $GeneName2 (keys %{$fus_annotations{$GeneName1}})
#	{
#		foreach my $GeneBreakpoin1 (keys %{$fus_annotations{$GeneName1}{$GeneName2}})
#		{
#			foreach my $GeneBreakpoin2 (keys %{$fus_annotations{$GeneName1}{$GeneName2}{$GeneBreakpoin1}})
#			{
#				printf "%s\t%s\t%s\t%s\t%s\n", $GeneName1, $GeneName2, $GeneBreakpoin1, $GeneBreakpoin2, $fus_annotations{$GeneName1}{$GeneName2}{$GeneBreakpoin1}{$GeneBreakpoin2};
#			}
#		}
#	}
#}

open IN, "<".$chimerascan_report or die $!; 

my $header = "FusionID" . "\t";
$header   .= "Sample_Name" . "\t";
$header   .= "Program" . "\t";
$header   .= "Tot/span_reads" . "\t";
$header   .= "Split_reads" . "\t";
$header   .= "Chr1" . "\t";
$header   .= "Chr2" . "\t";
$header   .= "Gene_Start1" . "\t";
$header   .= "Gene_End1" . "\t";
$header   .= "Gene_Start2" . "\t";
$header   .= "Gene_End2" . "\t";
$header   .= "Strand1" . "\t";
$header   .= "Strand2" . "\t";
$header   .= "Gene_Name1" . "\t";
$header   .= "Gene_Name2" . "\t";
$header   .= "Gene_Breakpoint1" . "\t";
$header   .= "Gene_Breakpoint2" . "\t";
$header   .= "Gene_ID1" . "\t";
$header   .= "Gene_ID2" . "\t";
$header   .= "Sample_Type" . "\t";
$header   .= "Sample_occurrency_list" . "\t";
$header   .= "Sample_Type_occurency_list" . "\t";
$header   .= "Kinase_info" . "\t";
$header   .= "Gene_Name1" . "\t";
$header   .= "Gene_Name2" . "\t";
$header   .= "Gene_Breakpoint1" . "\t";
$header   .= "Gene_Breakpoint2" . "\t";
$header   .= "Transcript_ID1" . "\t";
$header   .= "Transcript_ID2" . "\t";
$header   .= "Reading_Frame_Info" . "\t";
$header   .= "Fused_Nucleotide_Seq" . "\t";
$header   .= "Protein_Start1" . "\t";
$header   .= "Protein_End1" . "\t";
$header   .= "Protein_Start2" . "\t";
$header   .= "Protein_End2" . "\t";
$header   .= "Protein_Sequence" . "\t";
$header   .= "Exon_Gene1" . "\t";
$header   .= "Exon_Gene2" . "\t";
$header   .= "Breakpoint_Region1" . "\t";
$header   .= "Breakpoint_Region2" . "\t";
$header   .= "Conserved_Domain1" . "\t";
$header   .= "Lost_Domain1" . "\t";
$header   .= "Conserved_Domain2" . "\t";
$header   .= "Lost_Domain2" . "\n";

printf "%s", $header;


while(my $line = <IN>)
{
	chomp($line);
	
	my @fields = split(/\t/, $line);

	my $g1_strand = $fields[11];
	my $g2_strand = $fields[12];
	my $g1_name = $fields[13];
	my $g2_name = $fields[14];
	my $g1_start = $fields[7];
	my $g1_stop = $fields[8];
	my $g2_start = $fields[9];
	my $g2_stop = $fields[10];
	my $program = $fields[2];

	my $g1_bp;
	my $g2_bp;
	if($program eq "chimerascan")
	{
		if($g1_strand eq "+")
		{
			$g1_bp=$g1_stop;
			$g1_stop++;
			$g1_bp++;
		}
		else
		{
			$g1_bp=$g1_start;
			$g1_start++;
			$g1_bp++;
		}
		if($g2_strand eq "+")
		{
			$g2_bp=$g2_start;
			$g2_start++;
			$g2_bp++;
		}
		else
		{
			$g2_bp=$g2_stop;
			$g2_stop++;
			$g2_bp++;
		}	
	}
	elsif($program eq "bellerophontes")
	{
		if($g1_strand eq "+")
		{
			$g1_bp=$g1_stop;
			$g1_stop++;
			$g1_bp++;
		}
		else
		{
			$g1_bp=$g1_start;
			$g1_start++;
			$g1_bp++;
		}
		if($g2_strand eq "+")
		{
			$g2_bp=$g2_start;
			$g2_start++;
			$g2_bp++;
		}
		else
		{
			$g2_bp=$g2_stop;
			$g2_stop++;
			$g2_bp++;
		}	
	}
	elsif($program eq "general")
	{
		if($g1_strand eq "+")
		{
			$g1_bp=$g1_stop;
			$g1_stop++;
			$g1_bp++;
		}
		else
		{
			$g1_bp=$g1_start;
			$g1_start++;
			$g1_bp++;
		}
		if($g2_strand eq "+")
		{
			$g2_bp=$g2_start;
			$g2_start++;
			$g2_bp++;
		}
		else
		{
			$g2_bp=$g2_stop;
			$g2_stop++;
			$g2_bp++;
		}	
	}
	else
	{
		$g1_bp = $fields[15];
		$g2_bp = $fields[16];
	}
	my $null_str = "NULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL";
	
	#printf "%s\t%s\t%s\t%s\n", $g1_name, $g2_name, $g1_bp, $g2_bp;
	
	if(defined($fus_annotations{$g1_name}{$g2_name}{$g1_bp}{$g2_bp}))
	{
		printf "%s\t%s\n", $line, $fus_annotations{$g1_name}{$g2_name}{$g1_bp}{$g2_bp};
	}
	else
	{
		printf "%s\t%s\n", $line, $null_str;
	}

}

close(IN);














