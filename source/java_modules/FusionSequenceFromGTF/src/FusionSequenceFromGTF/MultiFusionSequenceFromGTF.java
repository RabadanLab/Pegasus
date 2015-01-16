package FusionSequenceFromGTF;
import example.UniprotWebServiceClient;
import example.UniprotWebServiceClient.ParameterNameValue;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;

import gnu.getopt.Getopt;

public class MultiFusionSequenceFromGTF {
	
	public MultiFusionSequenceFromGTF(
			String GTFfilepath, 
			String input_file, 
			String reference_FASTA_file, 
			boolean verbose, 
			boolean print_exons,
			boolean print_domains)
	{
		try 
		{

			ArrayList<String> domains = new ArrayList<String>();
			ArrayList<String> domains_all = new ArrayList<String>();
			ArrayList<String> domains_excluded = new ArrayList<String>();
			
			GTFInterface gtfinterface = new GTFInterface(GTFfilepath);
			String five_prime_gene_info = "";
			String three_prime_gene_info = "";
			
			System.out.print(PrintHeadings());
			
			FileReader f_in = new FileReader(new File(input_file));
			BufferedReader br_in = new BufferedReader(f_in);
			FusionSequenceFromGTF fsfgtf;
			String line = "";
			while((line = br_in.readLine()) != null)
			{
				StringBuffer sb = new StringBuffer();
				String[] fields = line.split("\t");
				five_prime_gene_info = fields[0];
				three_prime_gene_info = fields[1];
				fsfgtf = new FusionSequenceFromGTF(gtfinterface, five_prime_gene_info, three_prime_gene_info, reference_FASTA_file, verbose, print_exons);
				 
				int max_five_prime_amino_end = -1;
				int max_three_prime_amino_end = -1;
				boolean all_frameshifted = true;
				boolean all_five_bp_introns=true;
				boolean all_three_bp_introns=true;
				ArrayList<FusionSequenceInfo> fusionsequenceinfolist_complete = fsfgtf.getFusionSequencesInfo();
				ArrayList<FusionSequenceInfo> fusionsequenceinfolist_codingcoding = new ArrayList<FusionSequenceInfo>();
				ArrayList<FusionSequenceInfo> fusionsequenceinfolist_UTR = new ArrayList<FusionSequenceInfo>();
				ArrayList<FusionSequenceInfo> fusionsequenceinfolist_Intergenic = new ArrayList<FusionSequenceInfo>();
				
				for (FusionSequenceInfo f : fusionsequenceinfolist_complete)
				{
					if(f.GetType5P().equals("Coding") && f.GetType3P().equals("Coding"))
					{
						fusionsequenceinfolist_codingcoding.add(f);
					}
					else if(
							(f.GetType5P().equals("5UTR") && f.GetType3P().equals("Coding")) ||
							(f.GetType5P().equals("3UTR") && f.GetType3P().equals("Coding")) ||
							(f.GetType5P().equals("Coding") && f.GetType3P().equals("5UTR")) ||
							(f.GetType5P().equals("Coding") && f.GetType3P().equals("3UTR")) ||
							(f.GetType5P().equals("5UTR") && f.GetType3P().equals("5UTR")) ||
							(f.GetType5P().equals("5UTR") && f.GetType3P().equals("3UTR")) ||
							(f.GetType5P().equals("3UTR") && f.GetType3P().equals("5UTR")) ||
							(f.GetType5P().equals("3UTR") && f.GetType3P().equals("3UTR"))
					) 
					{
						// either one or both the genes is UTR and the other is a coding or UTR, but not a intergenic
						fusionsequenceinfolist_UTR.add(f);
					}
					else if(
							f.GetType5P().equals("Intergenic") ||
							f.GetType3P().equals("Intergenic")
					)
					{
						// one of the gene is an Intergenic
						fusionsequenceinfolist_Intergenic.add(f);
					}
				}
				for (FusionSequenceInfo f : fusionsequenceinfolist_codingcoding)
				{
					if(f.GetFrame().equals("InFrame"))
					{
						all_frameshifted = false;
					}
					if(f.GetExomeNumber1()!=0)
					{
						all_five_bp_introns = false;
					}
					if(f.GetExomeNumber2()!=0)
					{
						all_three_bp_introns = false;
					}
//					if(f.GetFivePrimeAminoAcidStop() > max_five_prime_amino_end)
//					{
//						max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
//					}
//					if(f.GetThreePrimeAminoAcidStop() > max_three_prime_amino_end)
//					{
//						max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
//					}
				}

				//FusionSequenceInfo f_sel = null;
				ArrayList<FusionSequenceInfo> f_sel = new ArrayList<FusionSequenceInfo>();

				//I print the longest transcript
				for (FusionSequenceInfo f : fusionsequenceinfolist_codingcoding)
				{
					if(all_frameshifted)
					{
						if(all_five_bp_introns && all_three_bp_introns)
						{
							if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
							{
								max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
								max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
								//f_sel = f;
								f_sel.add(f);
							}
						}
						else if(!all_five_bp_introns && all_three_bp_introns)
						{
							if(f.GetExomeNumber1()!=0)
							{
								if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
								{
									max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
									max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
									//f_sel = f;
									f_sel.add(f);
								}
							}
						}
						else if(all_five_bp_introns && !all_three_bp_introns)
						{
							if(f.GetExomeNumber2()!=0)
							{
								if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
								{
									max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
									max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
									//f_sel = f;
									f_sel.add(f);
								}
							}
						}
						else if(!all_five_bp_introns && !all_three_bp_introns)
						{
							if(f.GetExomeNumber1()!=0 && f.GetExomeNumber2()!=0)
							{
								if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
								{
									max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
									max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
									//f_sel = f;
									f_sel.add(f);
								}
							}
						}
					}
					else //there are both frameshift and inframe, I chose Inframe
					{
						if(f.GetFrame().equals("InFrame"))
						{
							if(all_five_bp_introns && all_three_bp_introns)
							{
								if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
								{
									max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
									max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
									//f_sel = f;
									f_sel.add(f);
								}
							}
							else if(!all_five_bp_introns && all_three_bp_introns)
							{
								if(f.GetExomeNumber1()!=0)
								{
									if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
									{
										max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
										max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
										//f_sel = f;
										f_sel.add(f);
									}
								}
							}
							else if(all_five_bp_introns && !all_three_bp_introns)
							{
								if(f.GetExomeNumber2()!=0)
								{
									if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end && f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
									{
										max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
										max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
										//f_sel = f;
										f_sel.add(f);
									}
								}
							}
							else if(!all_five_bp_introns && !all_three_bp_introns)
							{
								if(f.GetExomeNumber1()!=0 && f.GetExomeNumber2()!=0)
								{
									if(f.GetFivePrimeAminoAcidStop() >= max_five_prime_amino_end &&	f.GetThreePrimeAminoAcidStop() >= max_three_prime_amino_end)
									{
										max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
										max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
										//f_sel = f;
										f_sel.add(f);
									}
								}
							}
						}
					}
				}
				
				fusionsequenceinfolist_codingcoding = new ArrayList<FusionSequenceInfo>();
//				if(f_sel!=null)
//					fusionsequenceinfolist_codingcoding.add(f_sel);
				if(f_sel!=null)
					fusionsequenceinfolist_codingcoding.addAll(f_sel);
				
				//I print the longest transcript
				for (FusionSequenceInfo f : fusionsequenceinfolist_codingcoding)
				{
					sb.append(five_prime_gene_info.split(":")[0]);
					sb.append("\t");
					sb.append(three_prime_gene_info.split(":")[0]);
					sb.append("\t");
					sb.append(five_prime_gene_info.split(":")[4]);
					sb.append("\t");
					sb.append(three_prime_gene_info.split(":")[4]);
					sb.append("\t");
					if(print_exons)
					{
						//System.out.print(f.toString());
						sb.append(f.toString());
					}
					else
					{
						//System.out.print(f.toStringNoExons());	
						sb.append(f.toStringNoExons());
					}
					sb.append("\t");
					sb.append("CDS");
					sb.append("\t");
					sb.append("CDS");
					
					
					// START - Domain
					if(print_domains)
					{
						//System.out.println("\t");
						//domain print information
						sb.append("\t");
						domains = GetOverlappingProteinDomain(f.GetFivePrimeTranscritp_Id(), f.GetFivePrimeAminoAcidStart(), f.GetFivePrimeAminoAcidStop());
						domains_all = GetAllProteinDomain(f.GetFivePrimeTranscritp_Id());
						domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
						if(domains.size() > 0)
						{
							for (String d : domains)
							{
								//System.out.print(d + ",");
								sb.append(d + ",");
							}
						}
						else
						{
							//System.out.print("NO_DOMAIN");
							sb.append("NO_DOMAIN");
						}
						sb.append("\t");
						
						if(domains_excluded.size() > 0)
						{
							for (String d : domains_excluded)
							{
								//System.out.print(d + ",");
								sb.append(d + ",");
							}
						}
						else
						{
							//System.out.print("NO_DOMAIN");
							sb.append("NO_DOMAIN");
						}
	
						//System.out.print("\t");
						sb.append("\t");
						domains = GetOverlappingProteinDomain(f.GetThreePrimeTranscritp_Id(), f.GetThreePrimeAminoAcidStart(), f.GetThreePrimeAminoAcidStop());
						domains_all = GetAllProteinDomain(f.GetThreePrimeTranscritp_Id());
						domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
	
						if(domains.size() > 0)
						{
							for (String d : domains)
							{
								//System.out.print(d + ",");
								sb.append(d + ",");
							}
						}
						else
						{
							//System.out.print("NO_DOMAIN");
							sb.append("NO_DOMAIN");
						}					
						sb.append("\t");
						
						if(domains_excluded.size() > 0)
						{
							for (String d : domains_excluded)
							{
								//System.out.print(d + ",");
								sb.append(d + ",");
							}
						}
						else
						{
							//System.out.print("NO_DOMAIN");
							sb.append("NO_DOMAIN");
						}
					}
					//END - Domain
					
					//System.out.print("\n");
					sb.append("\n");
					System.out.print(sb.toString());
				}
				
				
				
				
				/*
				 * Print UTR information if exists
				 * and only if a coding coding does not exist
				 * */
				if(fusionsequenceinfolist_codingcoding.size()==0)
				{
					if(fusionsequenceinfolist_UTR.size()!=0)
					{
						FusionSequenceInfo f = fusionsequenceinfolist_UTR.get(0);
						sb = new StringBuffer();
						
						if(f.GetType5P().equals("5UTR") && f.GetType3P().equals("Coding"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("5UTR");
							sb.append("\t");
							sb.append("CDS");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							
							domains = GetOverlappingProteinDomain(f.GetThreePrimeTranscritp_Id(), f.GetThreePrimeAminoAcidStart(), f.GetThreePrimeAminoAcidStop());
							domains_all = GetAllProteinDomain(f.GetThreePrimeTranscritp_Id());
							domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
							if(domains.size() > 0)
							{
								for (String d : domains)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}
							sb.append("\t");
							
							if(domains_excluded.size() > 0)
							{
								for (String d : domains_excluded)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}

//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
						}
						
						if(f.GetType5P().equals("Coding") && f.GetType3P().equals("5UTR"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("CDS");
							sb.append("\t");
							sb.append("5UTR");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}

						if(f.GetType5P().equals("Coding") && f.GetType3P().equals("3UTR"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("CDS");
							sb.append("\t");
							sb.append("3UTR");

							sb.append("\t");
//							sb.append("DOMAIN_UNK");
//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
							
							domains = GetOverlappingProteinDomain(f.GetFivePrimeTranscritp_Id(), f.GetFivePrimeAminoAcidStart(), f.GetFivePrimeAminoAcidStop());
							domains_all = GetAllProteinDomain(f.GetFivePrimeTranscritp_Id());
							domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
							if(domains.size() > 0)
							{
								for (String d : domains)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}
							sb.append("\t");
							
							if(domains_excluded.size() > 0)
							{
								for (String d : domains_excluded)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}

							
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}

						if(f.GetType5P().equals("3UTR") && f.GetType3P().equals("Coding"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("3UTR");
							sb.append("\t");
							sb.append("CDS");
							
							sb.append("\t");
							domains = GetOverlappingProteinDomain(f.GetFivePrimeTranscritp_Id(), f.GetFivePrimeAminoAcidStart(), f.GetFivePrimeAminoAcidStop());
							domains_all = GetAllProteinDomain(f.GetFivePrimeTranscritp_Id());
							domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
							if(domains.size() > 0)
							{
								for (String d : domains)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}
							sb.append("\t");
							
							if(domains_excluded.size() > 0)
							{
								for (String d : domains_excluded)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}


//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}

						if(f.GetType5P().equals("3UTR") && f.GetType3P().equals("5UTR"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("3UTR");
							sb.append("\t");
							sb.append("5UTR");
							sb.append("\t");

							//entire five prime domain
							domains = GetOverlappingProteinDomain(f.GetFivePrimeTranscritp_Id(), f.GetFivePrimeAminoAcidStart(), f.GetFivePrimeAminoAcidStop());
							domains_all = GetAllProteinDomain(f.GetFivePrimeTranscritp_Id());
							domains_excluded = GetNotOverlappingProteinDomains(domains_all, domains);
							if(domains.size() > 0)
							{
								for (String d : domains)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}
							sb.append("\t");
							
							if(domains_excluded.size() > 0)
							{
								for (String d : domains_excluded)
								{
									//System.out.print(d + ",");
									sb.append(d + ",");
								}
							}
							else
							{
								//System.out.print("NO_DOMAIN");
								sb.append("NO_DOMAIN");
							}

//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
//							sb.append("\t");
//							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						
						System.out.print(sb.toString() + "\n");
					}
					else if(fusionsequenceinfolist_Intergenic.size()!=0)
					{
						// the fusionsequenceinfolist_UTR at this point is empty
						// if fusionsequenceinfolist_Intergenic is not empty check and print the info
						FusionSequenceInfo f = fusionsequenceinfolist_Intergenic.get(0);
						sb = new StringBuffer();
						
						if(f.GetType5P().equals("Coding") && f.GetType3P().equals("Intergenic"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("CDS");
							sb.append("\t");
							sb.append("Intergenic");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						else if(f.GetType5P().equals("5UTR") && f.GetType3P().equals("Intergenic"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("5UTR");
							sb.append("\t");
							sb.append("Intergenic");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						else if(f.GetType5P().equals("3UTR") && f.GetType3P().equals("Intergenic"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("3UTR");
							sb.append("\t");
							sb.append("Intergenic");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						if(f.GetType5P().equals("Intergenic") && f.GetType3P().equals("Coding"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("Intergenic");
							sb.append("\t");
							sb.append("CDS");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						else if(f.GetType5P().equals("Intergenic") && f.GetType3P().equals("5UTR"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("Intergenic");
							sb.append("\t");
							sb.append("5UTR");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
						else if(f.GetType5P().equals("Intergenic") && f.GetType3P().equals("3UTR"))
						{
							sb.append(five_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[0]);
							sb.append("\t");
							sb.append(five_prime_gene_info.split(":")[4]);
							sb.append("\t");
							sb.append(three_prime_gene_info.split(":")[4]);
							sb.append("\t");
							if(print_exons)
							{
								//System.out.print(f.toString());
								sb.append(f.toString());
							}
							else
							{
								//System.out.print(f.toStringNoExons());	
								sb.append(f.toStringNoExons());
							}
							sb.append("\t");
							sb.append("Intergenic");
							sb.append("\t");
							sb.append("3UTR");

							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
							sb.append("\t");
							sb.append("DOMAIN_UNK");
						}
					}
					else
					{
						/* Both the coding-coding list and the UTR are empty. */
					}
					System.out.print(sb.toString() + "\n");
				}
				
			}
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private ArrayList<String> GetProteinId(String transcript_id)
	{
		UniprotWebServiceClient uwsc = new UniprotWebServiceClient();		
		TreeMap<String, String> p_ids = new TreeMap<String, String>();
		ArrayList<String> tmp = new ArrayList<String>();

		try {
			String s;
			s = uwsc.run_public("mapping", new ParameterNameValue[] {
	            	new ParameterNameValue("from", "ENSEMBL_TRS_ID"),
	            	new ParameterNameValue("to", "ACC"),
	                new ParameterNameValue("format", "tab"),
	                new ParameterNameValue("query", transcript_id)});
			//System.out.println(s);

			String[] lines = s.split("\n");
			for (int i=1; i<lines.length; i++)
			{
				String[] fields = lines[i].split("\t");
				p_ids.put(fields[1], null);
			}

		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for (String s : p_ids.keySet())
		{
			tmp.add(s);
		}
		
		return tmp;
	}
	
	private ArrayList<String> GetNotOverlappingProteinDomains(ArrayList<String> all_domains, ArrayList<String> overlapping_domains)
	{
		ArrayList<String> not_overlapping = new ArrayList<String>();
		for (String all_domain : all_domains)
		{
			boolean found = false;
			for(String overlapping_domain : overlapping_domains)
			{
				if (overlapping_domain.compareTo(all_domain)==0)
				{
					found = true;
					break;
				}
			}
			if (!found)
			{
				not_overlapping.add(all_domain);
			}
		}
		
		return not_overlapping;
	}

	private ArrayList<String> GetAllProteinDomain(String transcript_id)
	{

		UniprotWebServiceClient uwsc = new UniprotWebServiceClient();
		ArrayList<String> domains = new ArrayList<String>();
		ArrayList<String> prot_ids = GetProteinId(transcript_id);

		try {
			for (String prot_id : prot_ids)
			{
				String s;
				
				s = uwsc.run_public("batch", new ParameterNameValue[] {
				  	  new ParameterNameValue("query", prot_id),
				    	  new ParameterNameValue("columns", "features"),
				    	  new ParameterNameValue("columns", "families"),
				    	  new ParameterNameValue("columns", "comments"),
				    	  new ParameterNameValue("format", "gff")});
				//System.out.println(s);
				String[] lines = s.split("\n");
				for (int i=0; i<lines.length; i++)
				{
					if(lines[i].matches("^#.*"))
					{
						continue;
					}
					String[] fields = lines[i].split("\t");
					
					if(fields[2].compareTo("Domain")==0)
					{
						domains.add("Domain=" + fields[8].split("Note=")[1]);
					}					
					else if(fields[2].compareTo("Topological domain")==0)
					{
						domains.add("Topological domain=" + fields[8].split("Note=")[1]);
					}					
					else if(fields[2].compareTo("Repeat")==0)
					{
						domains.add("Repeat=" + fields[8].split("Note=")[1]);
					}					
					else if(fields[2].compareTo("Region")==0)
					{
						domains.add("Region=" + fields[8].split("Note=")[1]);
					}

				}
				
			}

		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return domains;
	}

	private ArrayList<String> GetOverlappingProteinDomain(String transcript_id, int p_start, int p_end)
	{
		int prot_start;
		int prot_end;
		if (p_start > p_end)
		{
			prot_start = p_end;
			prot_end = p_start;
		}
		else 
		{
			prot_start = p_start;
			prot_end = p_end;
		}
		UniprotWebServiceClient uwsc = new UniprotWebServiceClient();
		ArrayList<String> domains = new ArrayList<String>();
		ArrayList<String> prot_ids = GetProteinId(transcript_id);

		try {
			for (String prot_id : prot_ids)
			{
				String s;
				
				s = uwsc.run_public("batch", new ParameterNameValue[] {
				  	  new ParameterNameValue("query", prot_id),
				    	  new ParameterNameValue("columns", "features"),
				    	  new ParameterNameValue("columns", "families"),
				    	  new ParameterNameValue("columns", "comments"),
				    	  new ParameterNameValue("format", "gff")});
				//System.out.println(s);
				String[] lines = s.split("\n");
				for (int i=0; i<lines.length; i++)
				{
					if(lines[i].matches("^#.*"))
					{
						continue;
					}
					String[] fields = lines[i].split("\t");
					int gff_start = Integer.parseInt(fields[3]);
					int gff_stop = Integer.parseInt(fields[4]);
					
					if(fields[2].compareTo("Domain")==0)
					{
						if(gff_start <= prot_end && gff_stop >= prot_start)
						{
							domains.add("Domain=" + fields[8].split("Note=")[1]);
						}
					}									
					else if(fields[2].compareTo("Topological domain")==0)
					{
						if(gff_start <= prot_end && gff_stop >= prot_start)
						{
							domains.add("Topological domain=" + fields[8].split("Note=")[1]);
						}
					}					
					else if(fields[2].compareTo("Repeat")==0)
					{
						if(gff_start <= prot_end && gff_stop >= prot_start)
						{
							domains.add("Repeat=" + fields[8].split("Note=")[1]);
						}
					}					
					else if(fields[2].compareTo("Region")==0)
					{
						if(gff_start <= prot_end && gff_stop >= prot_start)
						{
							domains.add("Region=" + fields[8].split("Note=")[1]);
						}
					}
//					if(	fields[2].compareTo("Domain")==0 || 
//						fields[2].compareTo("Repeat")==0 || 
//						fields[2].compareTo("Region")==0)
//					{
//						if(gff_start <= prot_end && gff_stop >= prot_start)
//						{
//							domains.add(fields[8]);
//						}
//					}
				}
				
			}

		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return domains;
	}

	private String PrintHeadings()
	{
		StringBuffer sb = new  StringBuffer();
		
		sb.append("GeneName1");
		sb.append("\t");
		sb.append("GeneName2");
		sb.append("\t");
		sb.append("GeneBreakpoin1");
		sb.append("\t");
		sb.append("GeneBreakpoin2");
		sb.append("\t");
		sb.append("GeneTranscriptID1");
		sb.append("\t");
		sb.append("GeneTranscriptID2");
		sb.append("\t");
		sb.append("FrameType");
		sb.append("\t");
		sb.append("FusedSequence");
		sb.append("\t");
		sb.append("ProteinStart1");
		sb.append("\t");
		sb.append("ProteinStop1");
		sb.append("\t");
		sb.append("ProteinStart2");
		sb.append("\t");
		sb.append("ProteinStop2");
		sb.append("\t");
		sb.append("ProteinSequence");
		sb.append("\t");
		sb.append("ExonBreak1");
		sb.append("\t");
		sb.append("ExonBreak2");
		sb.append("\t");
		sb.append("BpRegion1");
		sb.append("\t");
		sb.append("BpRegion2");
		sb.append("\t");
		sb.append("DomainConserved1");
		sb.append("\t");
		sb.append("DomainLost1");
		sb.append("\t");
		sb.append("DomainConserved2");
		sb.append("\t");
		sb.append("DomainLost2");
		sb.append("\n");
		
		return sb.toString();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String GTFfilepath = "";
		String input_file = "";
		String reference_FASTA_file = "";
		boolean verbose = false;
		boolean print_exons = false;
		boolean print_domain = false;

        Getopt g = new Getopt("MultiFusionSequenceFromGTF", args, "evdg:r:i:");
        //
        int c;

        while ((c = g.getopt()) != -1)
        {
            switch(c)
            {
    			case 'e':
    				print_exons = true;
    				break;

    			case 'v':
    				verbose = true;
    				break;

        		case 'd':
        			print_domain = true;                	 
        			break;

        		case 'g':
        			GTFfilepath = g.getOptarg();                	 
        			break;

        		case 'i':
        			input_file = g.getOptarg();                	 
        			break;

        		case 'r':
            		reference_FASTA_file = g.getOptarg();                	 
            		break;

                 case '?':
                	 System.exit(1);
                	 break;

                 default:
                	 System.exit(1);
                	 System.out.println("default");
            }
        }
		
        if (	
        		GTFfilepath.compareTo("")==0 || 
        		input_file.compareTo("")==0 || 
        		reference_FASTA_file.compareTo("")==0
        	)
        {
        	System.err.println("Error!\nusage: java -jar MultiFusionSequenceFromGTF" +
        			" -g [GTF file list] " +
        			" -i [Fusion Input File ] " +
        			" -r [reference FASTA file] " +
        			" -e [Print Exons (Optional)]" +
        			" -d [Print Domains (Optional)]" +
        			" -v [Verbose (Optional)]"
        	);
        	System.exit(1);
        }

        MultiFusionSequenceFromGTF mfsfgtf = new MultiFusionSequenceFromGTF(GTFfilepath, input_file, reference_FASTA_file, verbose, print_exons, print_domain);

		
	}
	
	

}
