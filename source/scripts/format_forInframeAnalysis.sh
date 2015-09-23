# fusion_file_list=resources/overlap_new_ALCL_tab_formatted.txt
fusion_file_list=$1
#cat fusion_file_list |grep $2 | grep $3| awk '
cat $fusion_file_list | awk '
{
	program=$3;
	g1_strand=$12;
	g2_strand=$13;
	g1_name=$14;
	g2_name=$15;
	g1_start=$8;
	g1_stop=$9;
	g2_start=$10;
	g2_stop=$11;
	g1_bp=$16;
	g2_bp=$17;

	if(program=="defuse")
	{
		printf g1_name":"g1_strand":"g1_start":"g1_stop":"g1_bp"\t"g2_name":"g2_strand":"g2_start":"g2_stop":"g2_bp"\n"
	}
	else
	{
		if(g1_strand=="+")
		{
			g1_stop++;
			g1_bp++;
		}
		else
		{
			g1_start++;
			g1_bp++;
		}
		if(g2_strand=="+")
		{
			g2_start++;
			g2_bp++;
		}
		else
		{
			g2_stop++;
			g2_bp++;
		}
		printf g1_name":"g1_strand":"g1_start":"g1_stop":"g1_bp"\t"g2_name":"g2_strand":"g2_start":"g2_stop":"g2_bp"\n"
	}

}'
