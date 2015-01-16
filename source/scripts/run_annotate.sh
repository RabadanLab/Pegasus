#!/bin/bash

if [ $#	== 4 ];	then 

	# scripts directory
	d=$( dirname $( readlink -m $0 ) )
	
	split_in=$( readlink -m $1 )
	split_out=$( readlink -m $2 )
	basedir=$( readlink -m $3 )
	sge_params=$4

	n_files=$( ls ${split_in}/*.txt | wc -l);
	n_files=$((n_files-1));
	echo $n_files	
	mkdir ${d}/logs
	for i in `seq 0 1 $n_files`
	do
		#echo "qsub -N FuAnn_${i} -e ${d}/logs -o ${d}/logs ${d}/annotate.qsub ${split_out}/res_${i}.txt ${split_in}/splitted_file_${i}.txt ${split_in}/tmp_${i};"
		h=`qsub -N FuAnn_${i} -e ${d}/logs -o ${d}/logs ${sge_params} ${d}/annotate.qsub ${split_out}/res_${i}.txt ${split_in}/splitted_file_${i}.txt ${split_in}/tmp_${i} $basedir | cut -d" " -f3`;
		htot=$htot","$h
	done
	htot=`echo $htot | sed -e 's/^,//'`;

	qsub -N FuAnn_fin -hold_jid $htot -sync y -e ${d}/logs -o ${d}/logs ${sge_params} -b y echo 
else

	echo "error - need to supply an argument"

fi

