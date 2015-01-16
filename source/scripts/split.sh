fileinput=$1
splitted_dir=$2
n_split=$3

n=0;
c=0;
cat ${fileinput} | while read line
do 
	if [ $c -lt $n_split ] 
	then 
		echo $line | tr ' ' '\t' >> ${splitted_dir}/splitted_file_${n}.txt; 
		mkdir -p ${splitted_dir}/tmp_${n};
		c=$[$c+1]; 
	else 
		echo $line | tr ' ' '\t' >> ${splitted_dir}/splitted_file_${n}.txt; 
		mkdir -p ${splitted_dir}/tmp_${n};
		c=0;
		n=$[$n+1]; 
	fi
done

