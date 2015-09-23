#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

OPTIONS:
   -h      Show this message
   -r      Reference fasta file - eg. /path/to/hg19.fa
   -g      GTF Annotation file
   -i      Input gene list file
   -o      Output file
   -t      Temporary directory
   -s      Script directory
   -j      Jar directory
EOF
}

while getopts “hr:g:i:o:j:s:t:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         r)
	     #Reference fasta file - eg. /path/to/hg19.fa
             reference_fasta_file=$OPTARG
             ;;
         g)
	     #GTF Annotation file
             gtf_annotation=$OPTARG
             ;;
         i)
	     #Input gene list file
             input_geneList=$OPTARG
             ;;
         s)
	     #Script directory for finding executables
             scripts_dir=$OPTARG
             ;;
         j)
	     #Jar directory for finding executables
             jar_directory=$OPTARG
             ;;
         o)
	    #output_file
             output_file=$OPTARG
             ;;
         t)
	     #temp file
             tmp_folder=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [ -z $scripts_dir ] || [ -z $reference_fasta_file ] || [ -z $gtf_annotation ] || [ -z $input_geneList ] || [ -z $jar_directory ] || [ -z $output_file ] || [ -z $tmp_folder ]
then
     usage
     exit 1	  
fi

if [ ! -e $tmp_folder ]
then
	mkdir -p $tmp_folder
fi

if [ -e $tmp_folder/tmp_fusout.txt ]
then
	rm $tmp_folder/tmp_fusout.txt
else
	touch $tmp_folder/tmp_fusout.txt
fi

n_genes=0;
tot_genes=$(cat $input_geneList | wc -l);
tmp_gtf=${tmp_folder}/tmp.gtf;
cat $input_geneList | while read line
do
	ngenes=$((ngenes+1));
	perc=$((ngenes*100/tot_genes))
	gene1=`echo $line | cut -d" " -f1 | cut -d":" -f1 | sed -e 's/,/\\\|/g'`;
	gene2=`echo $line | cut -d" " -f2 | cut -d":" -f1 | sed -e 's/,/\\\|/g'`;
	echo -ne "Elaborating: \t$gene1 - $gene2\t\tTotal Elaborated: $perc % \n";
	
	cat $gtf_annotation | grep "gene_name \"${gene1}\"\|gene_name \"${gene2}\"" > $tmp_gtf;
	#cat $gtf_annotation | grep "gene_name \"${gene1}\"\|gene_name \"${gene2}\"" > $(echo ${gene1} | sed -e 's/\\\|/_/g')_$(echo ${gene2} | sed -e 's/\\\|/_/g').gtf;
	
	echo $line | tr ' ' '\t' > $tmp_folder/tmp_input.txt
	#echo "java -Xmx1048m -jar $jar_directory/MultiFusionSequenceFromGTF.jar -g $tmp_gtf -r $reference_fasta_file -i $tmp_folder/tmp_input.txt -e -v -d 2> /dev/null | grep -v Breakpoin | perl $scripts_dir/compact_results.pl >> $tmp_folder/tmp_fusout.txt"
	java -Xmx1048m -jar $jar_directory/MultiFusionSequenceFromGTF.jar -g $tmp_gtf -r $reference_fasta_file -i $tmp_folder/tmp_input.txt -e -v -d 2> /dev/null | grep -v Breakpoin | perl $scripts_dir/compact_results.pl >> $tmp_folder/tmp_fusout.txt
done
echo
echo

#head -1 $tmp_folder/tmp_fusout.txt > $tmp_folder/tmp.txt

#header="GeneName1\tGeneName2\tGeneBreakpoin1\tGeneBreakpoin2\tGeneTranscriptID1\tGeneTranscriptID2\tFrameType\tFusedSequence\tProteinStart1\tProteinStop1\tProteinStart2\tProteinStop2\tProteinSequence\tExonBreak1\tExonBreak2\tBpRegion1\tBpRegion2\tDomainConserved1\tDomainLost1\tDomainConserved2\tDomainLost2"
header="";
echo $header > $tmp_folder/tmp.txt
cat $tmp_folder/tmp_fusout.txt | grep -v Breakpoint | sort -u >> $tmp_folder/tmp.txt

mv $tmp_folder/tmp.txt $output_file







