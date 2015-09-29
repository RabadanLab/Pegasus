#!/bin/bash

if [ $#	== 5 ];	then 

	split_in=$( readlink -m $1 )
	split_out=$( readlink -m $2 )
	basedir=$( readlink -m $3 )
	log_dir=$4
	sge_params=$5

	# scripts directory
	d=${basedir}/source/scripts
	

	n_files=$( ls ${split_in}/*.txt | wc -l);
	n_files=$((n_files-1));
	echo $n_files	
	mkdir -p ${log_dir}/run_ann
	for i in `seq 0 1 $n_files`
	do
		#echo "qsub -N FuAnn_${i} -e ${log_dir}/run_ann -o ${log_dir}/run_ann ${d}/annotate.qsub ${split_out}/res_${i}.txt ${split_in}/splitted_file_${i}.txt ${split_in}/tmp_${i};"
		h=`qsub -N FuAnn_${i} -e ${log_dir}/run_ann -o ${log_dir}/run_ann ${sge_params} ${d}/annotate.qsub ${split_out}/res_${i}.txt ${split_in}/splitted_file_${i}.txt ${split_in}/tmp_${i} $basedir | grep "has been submitted" | cut -d" " -f3  `;
		htot=$htot","$h
		echo "qsub -N FuAnn_${i} -e ${log_dir}/run_ann -o ${log_dir}/run_ann ${sge_params} ${d}/annotate.qsub ${split_out}/res_${i}.txt ${split_in}/splitted_file_${i}.txt ${split_in}/tmp_${i} $basedir"
		echo $htot
	done
	htot=`echo $htot | sed -e 's/^,//' -e 's/variable //g'`;
	echo $htot

	qsub -N FuAnn_fin -hold_jid $htot -sync y -e ${log_dir}/run_ann -o ${log_dir}/run_ann ${sge_params} -b y echo 
else

	echo "error - need to supply an argument"

fi

