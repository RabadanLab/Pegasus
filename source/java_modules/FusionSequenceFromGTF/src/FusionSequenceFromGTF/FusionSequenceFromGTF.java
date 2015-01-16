package FusionSequenceFromGTF;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.util.SequenceUtil;

import javagene.seq.AminoAcid;
import javagene.seq.Location;

import gnu.getopt.Getopt;

public class FusionSequenceFromGTF {
	private ArrayList<FusionSequenceInfo> fusion_sequences_info;

	public FusionSequenceFromGTF(
		GTFInterface gtfi,
		String fiveprime_gene_info,
		String threeprime_gene_info,
		String referenceFastafile,
		boolean verbose,
		boolean print_exons
	)
	{
		fusion_sequences_info = new ArrayList<FusionSequenceInfo>();
		File genomeReferenceFastaFile = new File(referenceFastafile);
		String genomeReferenceIndex = genomeReferenceFastaFile.getAbsolutePath() + ".fai";
		File genomeReferenceIndexFile = new File(genomeReferenceIndex);

		FastaSequenceIndex fastaseqindex = new FastaSequenceIndex(genomeReferenceIndexFile);
		
		IndexedFastaSequenceFile indexedfasta = new IndexedFastaSequenceFile(genomeReferenceFastaFile, fastaseqindex);
		
		if (indexedfasta.isIndexed())
		{
			System.err.println("File " + " has been indexed.");				
		}	
		else
		{
			System.err.println("FusionSequenceFromGTF: Error. Unable to index the " + genomeReferenceFastaFile + " file.");
		}

		GeneInfo fiveprimegene = new GeneInfo(fiveprime_gene_info);
		
		GeneInfo threeprimegene = new GeneInfo(threeprime_gene_info);
		
		String[] gene_names1 = fiveprimegene.getGeneName().split(",");
		
		// these two vector are meant to mantain the list of five prime and three prime genes
		// that are in the UTR region
		Vector<String> fp_gene_fp_UTR_transcripts = new Vector<String>();
		Vector<String> fp_gene_tp_UTR_transcripts = new Vector<String>();
		Vector<String> tp_gene_fp_UTR_transcripts = new Vector<String>();
		Vector<String> tp_gene_tp_UTR_transcripts = new Vector<String>();
		
		Vector<String> five_prime_transcripts;
		if(gene_names1.length != 1)
		{
			five_prime_transcripts = gtfi.getTranscriptIdByMultipleGeneName(fiveprimegene.getGeneName());
		}
		else
		{
			five_prime_transcripts = gtfi.getTranscriptIdByGeneName(fiveprimegene.getGeneName());
		}
		Vector<String> tmp = new Vector<String>();
		for (String five_prime_t : five_prime_transcripts)
		{
			if(!isCodingFusionSequence(five_prime_t, fiveprimegene, gtfi))
			{
				if(is5UtrFusionSequence(five_prime_t, fiveprimegene, gtfi))
				{
					fp_gene_fp_UTR_transcripts.add(five_prime_t);
//					System.out.println("DEBUG: 5' UTR for fiveprimegene");
//					if(fiveprimegene.getGeneStrand().equals("+"))
//					{
//						String s = getFpGeneFpFUTRposSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);		
//						System.out.println("DEBUG: 5' UTR pos sequence transcript " + five_prime_t + ": " + s);
//					}
//					else
//					{
//						String s = getFpGeneFpFUTRnegSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);						
//						System.out.println("DEBUG: 5' UTR neg sequence: " + s);
//					}
				}
				else if(is3UtrFusionSequence(five_prime_t, fiveprimegene, gtfi))
				{
					fp_gene_tp_UTR_transcripts.add(five_prime_t);
//					System.out.println("DEBUG: 3' UTR for 5' transcript " + five_prime_t + " of " + fiveprimegene.getGeneName());
				}
				else if(verbose)
					System.out.println("Breakpoint does not fall into a coding region for 5 prime Transcript " + five_prime_t + " of " + fiveprimegene.getGeneName());
			}
			else
			{
				tmp.add(five_prime_t);
			}
		}
		five_prime_transcripts = new Vector<String>();
		for (String s : tmp)
		{
			five_prime_transcripts.add(s);
		}
		
		
		
		
		
		String[] gene_names2 = threeprimegene.getGeneName().split(",");
		
		Vector<String> three_prime_transcripts;
		if(gene_names2.length != 1)
		{
			three_prime_transcripts = gtfi.getTranscriptIdByMultipleGeneName(threeprimegene.getGeneName());
		}
		else
		{
			three_prime_transcripts = gtfi.getTranscriptIdByGeneName(threeprimegene.getGeneName());
		}
		tmp = new Vector<String>();
		for (String three_prime_t : three_prime_transcripts)
		{
			if(!isCodingFusionSequence(three_prime_t, threeprimegene, gtfi))
			{
				if(is5UtrFusionSequence(three_prime_t, threeprimegene, gtfi))
				{
					tp_gene_fp_UTR_transcripts.add(three_prime_t);
//					System.out.println("DEBUG: 5' UTR for threeprimegene");
//					if(threeprimegene.getGeneStrand().equals("+"))
//					{
//						String s = getTpGeneFpFUTRposSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);		
//						System.out.println("DEBUG: 5' UTR pos sequence transcript " + three_prime_t + ": " + s);
//					}
//					else
//					{
//						String s = getTpGeneFpFUTRnegSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);						
//						System.out.println("DEBUG: 5' UTR neg sequence: " + s);
//					}
				}
				else if(is3UtrFusionSequence(three_prime_t, threeprimegene, gtfi))
				{
					tp_gene_tp_UTR_transcripts.add(three_prime_t);
//					System.out.println("DEBUG: 3' UTR for 3' transcript " + three_prime_t + " of " + threeprimegene.getGeneName());
				}
				else if(verbose)
					System.out.println("Breakpoint does not fall into a coding region for 3 prime Transcript " + three_prime_t + " of " + threeprimegene.getGeneName());
			}
			else
			{
				tmp.add(three_prime_t);
			}
		}
		three_prime_transcripts = new Vector<String>();
		for (String s : tmp)
		{
			three_prime_transcripts.add(s);
		}
		

		/*
		 * this is the Coding Coding breakpoint case
		 * */
		for (String five_prime_t : five_prime_transcripts)
		{
			for (String three_prime_t : three_prime_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("Coding");
				fusion_sequence_info.SetType3P("Coding");
				
				
				//String five_seq = getFivePrimeSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
				fusion_sequence_info.SetFive_NucleotideSeq(getFivePrimeSequenceFromTranscriptId(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta));
				//String three_seq = getThreePrimeSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);
				fusion_sequence_info.SetThree_NucleotideSeq(getThreePrimeSequenceFromTranscriptId(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta));

				//String complete_threeprimeseq = getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta);
				fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq(getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta));
				
				StringBuffer sb = new StringBuffer();

				//System.out.println(five_prime_t + " " + three_prime_t + ":");
				//sb.append(five_prime_t + "\t" + three_prime_t);
				//sb.append("\t");
				if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0 && fusion_sequence_info.GetThree_NucleotideSeq().compareTo("")!=0 )
				{
					String chimeric_seq = fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq();
					if((chimeric_seq.length() % 3)==0)
					{
						//sb.append("InFrame\t");
						fusion_sequence_info.SetFrame("InFrame");
					}
					else
					{
						//sb.append("FrameShift\t");
						fusion_sequence_info.SetFrame("FrameShift");
					}
					
					fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
					fusion_sequence_info.SetExomeNumber2(gtfi.getExonNumberByTranscriptIdAndBreakPoint(three_prime_t, threeprimegene.getGeneBreakpoint()));
						
					String threeprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetThree_NucleotideSeq());					
					String fiveprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq());

					String threeprime_aminoacid_compl = AminoAcid.translate(fusion_sequence_info.GetCompleteThreeprime_NucleotideSeq());

					fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
					fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid.length());
					if(fusion_sequence_info.GetExomeNumber2()!=0)
					{
						// there's a 3' exonic breakpoint 
						fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid.length());
						fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());
					}
					else
					{
						// There is an intronic breakpoint in the 3' gene. Because of the intron presence the protein sequence is longer.
						// I have to compute the protein sequence only for the CDS region if I want to have the correct domain
						String threeprime_nucleotidesequence_cds = getThreePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta);
						String threeprime_aminoacid_cds = AminoAcid.translate(threeprime_nucleotidesequence_cds);
						fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid_cds.length());
						fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());

					}

					
					String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
					if(chimeric_aminoacid_sequence.contains("X"))
					{
						fusion_sequence_info.SetFrame("FrameShift");
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
					}
					else
					{
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
					}

					fusion_sequences_info.add(fusion_sequence_info);
					
				}
			}			
		}
		
		
		/*
		 * this is the 5G-5UTR / Coding breakpoint case
		 * */
		for (String five_prime_t : fp_gene_fp_UTR_transcripts)
		{
			for (String three_prime_t : three_prime_transcripts)
			{
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("5UTR");
				fusion_sequence_info.SetType3P("Coding");
				
				//String s = getFpGeneFpFUTRposSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
				String s;
				if (fiveprimegene.getGeneStrand().equals("+"))
				{
					s = getFpGeneFpFUTRposSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
				}
				else
				{
					s = getFpGeneFpFUTRnegSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
				}
				fusion_sequence_info.SetFive_NucleotideSeq(s);

				fusion_sequence_info.SetThree_NucleotideSeq(getThreePrimeSequenceFromTranscriptId(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta));

				fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq(getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta));
				
				if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0 && fusion_sequence_info.GetThree_NucleotideSeq().compareTo("")!=0 )
				{
					String chimeric_seq = fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq();
					
					String ATG_chimeric = FindATG_sequence(chimeric_seq);
					
					if((ATG_chimeric.length() % 3)==0)
					{
						fusion_sequence_info.SetFrame("InFrame");
					}
					else
					{
						fusion_sequence_info.SetFrame("FrameShift");
					}
					
					fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
					fusion_sequence_info.SetExomeNumber2(gtfi.getExonNumberByTranscriptIdAndBreakPoint(three_prime_t, threeprimegene.getGeneBreakpoint()));
						
					//String threeprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetThree_NucleotideSeq());					
					String threeprime_aminoacid = AminoAcid.translate(ATG_chimeric);					
					String fiveprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq());

					String threeprime_aminoacid_compl = AminoAcid.translate(fusion_sequence_info.GetCompleteThreeprime_NucleotideSeq());

					fusion_sequence_info.SetFivePrimeAminoAcidStart(-1);
					fusion_sequence_info.SetFivePrimeAminoAcidStop(-1);
					if(fusion_sequence_info.GetExomeNumber2()!=0)
					{
						// there's a 3' exonic breakpoint 
						fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid.length());
						fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());
					}
					else
					{
						// There is an intronic breakpoint in the 3' gene. Because of the intron presence the protein sequence is longer.
						// I have to compute the protein sequence only for the CDS region if I want to have the correct domain
						String threeprime_nucleotidesequence_cds = getThreePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta);
						String threeprime_nucleotidesequence_cds_fromATG=FindATG_sequence(threeprime_nucleotidesequence_cds);
						String threeprime_aminoacid_cds = AminoAcid.translate(threeprime_nucleotidesequence_cds_fromATG);
						fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid_cds.length());
						fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());

					}

					
					//String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
					String chimeric_aminoacid_sequence = AminoAcid.translate(ATG_chimeric);
					if(chimeric_aminoacid_sequence.contains("X"))
					{
						fusion_sequence_info.SetFrame("FrameShift");
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
					}
					else
					{
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
					}
					fusion_sequence_info.SetFive_NucleotideSeq("");
					fusion_sequence_info.SetThree_NucleotideSeq(ATG_chimeric);

					fusion_sequences_info.add(fusion_sequence_info);
					
				}
			}			
		}
		
		
		/*
		 * this is the 5G-5UTR / 3G-5UTR breakpoint case
		 * */
		for (String five_prime_t : fp_gene_fp_UTR_transcripts)
		{
			for (String three_prime_t : tp_gene_fp_UTR_transcripts)
			{				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("5UTR");
				fusion_sequence_info.SetType3P("5UTR");
				
//				String s = "";
//				if(fiveprimegene.getGeneStrand().equals("+"))
//				{
//					s = getFpGeneFpFUTRposSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
//				}
//				else
//				{
//					s = getFpGeneFpFUTRnegSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta);
//				}
//				fusion_sequence_info.SetFive_NucleotideSeq(s);
//
//				
//				if(threeprimegene.getGeneStrand().equals("+"))
//				{
//					s = getTpGeneFpFUTRposSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);
//				}
//				else
//				{
//					s = getTpGeneFpFUTRnegSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);
//				}
//				fusion_sequence_info.SetThree_NucleotideSeq(s + getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta));
//
//				fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq(getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta));
//				
//
//				if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0 && fusion_sequence_info.GetThree_NucleotideSeq().compareTo("")!=0 )
//				{
//					fusion_sequence_info.SetFrame("InFrame");
//					
//					fusion_sequence_info.SetExomeNumber1(-1);
//					fusion_sequence_info.SetExomeNumber2(-1);
//						
//					String threeprime_aminoacid_compl = AminoAcid.translate(fusion_sequence_info.GetCompleteThreeprime_NucleotideSeq());
//
//					fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
//					fusion_sequence_info.SetFivePrimeAminoAcidStop(threeprime_aminoacid_compl.length());
//
//					String chimeric_aminoacid_sequence = threeprime_aminoacid_compl;
//					if(chimeric_aminoacid_sequence.contains("X"))
//					{
//						fusion_sequence_info.SetFrame("FrameShift");
//						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
//						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
//						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
//						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
//					}
//					else
//					{
//						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
//					}
//
//					fusion_sequences_info.add(fusion_sequence_info);
//					
//				}
			}			
		}
		
		
		/*
		 * this is the Coding / 3G-3UTR breakpoint case
		 * */
		for (String five_prime_t : five_prime_transcripts)
		{
			for (String three_prime_t : tp_gene_tp_UTR_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("Coding");
				fusion_sequence_info.SetType3P("3UTR");
				
				
				fusion_sequence_info.SetFive_NucleotideSeq(getFivePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta));
				
				String s = "";
				if(threeprimegene.getGeneStrand().equals("+"))
				{
					s = "";
				}
				else
				{
					s = "";
				}
				fusion_sequence_info.SetThree_NucleotideSeq(s);
				
				if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0)
				{
					String chimeric_seq = fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq();
					if((chimeric_seq.length() % 3)==0)
					{
						//sb.append("InFrame\t");
						fusion_sequence_info.SetFrame("InFrame");
					}
					else
					{
						//sb.append("FrameShift\t");
						fusion_sequence_info.SetFrame("FrameShift");
					}
					
					fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
					fusion_sequence_info.SetExomeNumber2(-1);
						
					String fiveprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq());

					fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
					fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid.length());
					
					String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
					if(chimeric_aminoacid_sequence.contains("X"))
					{
						fusion_sequence_info.SetFrame("FrameShift");
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
					}
					else
					{
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
					}

					fusion_sequences_info.add(fusion_sequence_info);
					
				}
			}			
		}
		
		/*
		 * this is the Coding 3G-5UTR breakpoint case
		 * */
		for (String five_prime_t : five_prime_transcripts)
		{
			for (String three_prime_t : tp_gene_fp_UTR_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("Coding");
				fusion_sequence_info.SetType3P("5UTR");
				
				
				fusion_sequence_info.SetFive_NucleotideSeq(getFivePrimeSequenceFromTranscriptId(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta));
				
				String s = "";
				if(threeprimegene.getGeneStrand().equals("+"))
				{
					s = getTpGeneFpFUTRposSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);
				}
				else
				{
					s = getTpGeneFpFUTRnegSequenceFromTranscriptId(three_prime_t, threeprimegene, gtfi, indexedfasta);
				}

				fusion_sequence_info.SetThree_NucleotideSeq(s + getThreePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta));

				fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq(getThreePrimeSequenceFromTranscriptId_ONLYCDS(three_prime_t, threeprimegene, gtfi, indexedfasta));
				
				if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0 && fusion_sequence_info.GetThree_NucleotideSeq().compareTo("")!=0 )
				{
					String chimeric_seq = fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq();
					if((chimeric_seq.length() % 3)==0)
					{
						fusion_sequence_info.SetFrame("InFrame");
					}
					else
					{
						fusion_sequence_info.SetFrame("FrameShift");
					}
					
					fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
					fusion_sequence_info.SetExomeNumber2(-1);
						
					String fiveprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq());

					fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
					fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid.length());

					
					String threeprime_aminoacid_cds = AminoAcid.translate(fusion_sequence_info.GetCompleteThreeprime_NucleotideSeq());
					fusion_sequence_info.SetThreePrimeAminoAcidStart(1);
					fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_cds.length());
					
					String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
					if(chimeric_aminoacid_sequence.contains("X"))
					{
						fusion_sequence_info.SetFrame("FrameShift");
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
					}
					else
					{
						fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
					}

					fusion_sequences_info.add(fusion_sequence_info);
					
				}
			}			
		}
		
		
		/*
		 * this is the 5G-3UTR Coding breakpoint case
		 * Only the first gene is mantained
		 * */
		for (String five_prime_t : fp_gene_tp_UTR_transcripts)
		{
			for (String three_prime_t : three_prime_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("3UTR");
				fusion_sequence_info.SetType3P("Coding");
				
				String fiveprime_nucleotidesequence_cds = getFivePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta);
				String fiveprime_aminoacid_cds = AminoAcid.translate(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
				fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid_cds.length());
				fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
				fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
				fusion_sequence_info.SetFive_NucleotideSeq(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetThree_NucleotideSeq("");
				fusion_sequence_info.SetChimericAminoAcidSeq(fiveprime_aminoacid_cds);
				
//				fusion_sequence_info.SetFive_NucleotideSeq(getFivePrimeSequenceFromTranscriptId(five_prime_t, fiveprimegene, gtfi, indexedfasta));
//				fusion_sequence_info.SetThree_NucleotideSeq("");
//				
//				fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
//				fusion_sequence_info.SetExomeNumber2(-1);
//				
//				fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
//				fusion_sequence_info.SetFivePrimeAminoAcidStop(fusion_sequence_info.GetFive_NucleotideSeq().length());
//				fusion_sequence_info.SetThreePrimeAminoAcidStart(0);
//				fusion_sequence_info.SetThreePrimeAminoAcidStop(0);
//				
//				fusion_sequence_info.SetChimericAminoAcidSeq(AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq()));
//				
//				if((fusion_sequence_info.GetFive_NucleotideSeq().length() % 3)==0)
//				{
					fusion_sequence_info.SetFrame("InFrame");
//				}
//				else
//				{
//					fusion_sequence_info.SetFrame("FrameShift");
//				}
				fusion_sequences_info.add(fusion_sequence_info);
			}
		}
		
		/*
		 * this is the 5G-3UTR 3G-5UTR breakpoint case
		 * Only the first gene is mantained
		 * */
		for (String five_prime_t : fp_gene_tp_UTR_transcripts)
		{
			for (String three_prime_t : tp_gene_fp_UTR_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("3UTR");
				fusion_sequence_info.SetType3P("5UTR");
				
				String fiveprime_nucleotidesequence_cds = getFivePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta);
				String fiveprime_aminoacid_cds = AminoAcid.translate(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
				fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid_cds.length());
				fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
				fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
				fusion_sequence_info.SetFive_NucleotideSeq(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetThree_NucleotideSeq("");
				fusion_sequence_info.SetChimericAminoAcidSeq(fiveprime_aminoacid_cds);

				
				fusion_sequences_info.add(fusion_sequence_info);
			}
		}
		
		/*
		 * this is the 5G-3UTR 3G-5UTR breakpoint case
		 * Only the first gene is mantained
		 * */
		for (String five_prime_t : fp_gene_tp_UTR_transcripts)
		{
			for (String three_prime_t : tp_gene_tp_UTR_transcripts)
			{
				//HashSet<String> aminoacid_sequences = new HashSet<String>();
				
				FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
				
				fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
				fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
				fusion_sequence_info.SetType5P("3UTR");
				fusion_sequence_info.SetType3P("3UTR");
				
				String fiveprime_nucleotidesequence_cds = getFivePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta);
				String fiveprime_aminoacid_cds = AminoAcid.translate(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
				fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid_cds.length());
				fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
				fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
				fusion_sequence_info.SetFive_NucleotideSeq(fiveprime_nucleotidesequence_cds);
				fusion_sequence_info.SetThree_NucleotideSeq("");
				fusion_sequence_info.SetChimericAminoAcidSeq(fiveprime_aminoacid_cds);

				
				fusion_sequences_info.add(fusion_sequence_info);
			}
		}
		
		
		/*
		 * this is the case that the five prime gene breakpoint
		 * is in an intergenic region
		 * */
		if(
				five_prime_transcripts.size()==0 &&
				fp_gene_fp_UTR_transcripts.size()==0 &&
				fp_gene_tp_UTR_transcripts.size()==0
		)
		{
			if(three_prime_transcripts.size() != 0)
			{
				/*
				 * this is 5G-CDS / 3G-Intergenic breakpoint case
				 * Only the 5G gene is mantained
				 * */
				for (String three_prime_t : three_prime_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetThreePrimeTranscritp_Id(three_prime_t);
					fusion_sequence_info.SetType5P("Intergenic");
					fusion_sequence_info.SetType3P("Coding");
					
					fusion_sequence_info.SetFive_NucleotideSeq("");

					fusion_sequence_info.SetThree_NucleotideSeq(getThreePrimeSequenceFromTranscriptId(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta));

					fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq(getEntireSequenceFromTranscript(three_prime_t, threeprimegene, gtfi, indexedfasta));
					
					if(fusion_sequence_info.GetThree_NucleotideSeq().compareTo("")!=0 )
					{
						String chimeric_seq = fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq();
						if((chimeric_seq.length() % 3)==0)
						{
							fusion_sequence_info.SetFrame("InFrame");
						}
						else
						{
							fusion_sequence_info.SetFrame("FrameShift");
						}
						
						fusion_sequence_info.SetExomeNumber1(-1);
						fusion_sequence_info.SetExomeNumber2(gtfi.getExonNumberByTranscriptIdAndBreakPoint(three_prime_t, threeprimegene.getGeneBreakpoint()));
							
						String threeprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetThree_NucleotideSeq());					
						String fiveprime_aminoacid = "";

						String threeprime_aminoacid_compl = AminoAcid.translate(fusion_sequence_info.GetCompleteThreeprime_NucleotideSeq());

						fusion_sequence_info.SetFivePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetFivePrimeAminoAcidStop(-1);
						if(fusion_sequence_info.GetExomeNumber2()!=0)
						{
							// there's a 3' exonic breakpoint 
							fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid.length());
							fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());
						}
						else
						{
							// There is an intronic breakpoint in the 3' gene. Because of the intron presence the protein sequence is longer.
							// I have to compute the protein sequence only for the CDS region if I want to have the correct domain
							String threeprime_nucleotidesequence_cds = getThreePrimeSequenceFromTranscriptId_ONLYCDS(fusion_sequence_info.GetThreePrimeTranscritp_Id(), threeprimegene, gtfi, indexedfasta);
							String threeprime_aminoacid_cds = AminoAcid.translate(threeprime_nucleotidesequence_cds);
							fusion_sequence_info.SetThreePrimeAminoAcidStart(threeprime_aminoacid_compl.length()-threeprime_aminoacid_cds.length());
							fusion_sequence_info.SetThreePrimeAminoAcidStop(threeprime_aminoacid_compl.length());

						}

						
						String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
						if(chimeric_aminoacid_sequence.contains("X"))
						{
							fusion_sequence_info.SetFrame("FrameShift");
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
							fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
							fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
						}
						else
						{
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
						}

						fusion_sequences_info.add(fusion_sequence_info);
						
					}
					
				}
			}
			else if(fp_gene_fp_UTR_transcripts.size()!=0)
			{
				/*
				 * this is 5G-Intergenic / 3G-5UTR breakpoint case
				 * Only the 3G gene is mantained from the 5UTR of the 3G, the intergenic region does make difference
				 * 
				 * This case has higher priority compared to the 5G-5UTR
				 * */
				for (String five_prime_t : fp_gene_fp_UTR_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
					fusion_sequence_info.SetThreePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetFrame("UnSpec");
					fusion_sequence_info.SetType5P("Intergenic");
					fusion_sequence_info.SetType3P("5UTR");
					
					fusion_sequences_info.add(fusion_sequence_info);
				}

			}
			else if(tp_gene_tp_UTR_transcripts.size()!=0)
			{
				/*
				 * this is 5G-Intergenic / 3G-3UTR breakpoint case
				 * nothing is actually produced as coding sequence
				 * Only the 3G gene is mantained from the 5UTR, the intergenic region does make difference
				 * because a stop codon is introduced after the 3UTR
				 * 
				 * This case has higher priority compared to the 5G-5UTR
				 * */
				for (String five_prime_t : fp_gene_tp_UTR_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
					fusion_sequence_info.SetThreePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetFrame("UnSpec");
					fusion_sequence_info.SetType5P("Intergenic");
					fusion_sequence_info.SetType3P("3UTR");
					
					fusion_sequences_info.add(fusion_sequence_info);
				}

			}
		}

		

		/*
		 * this is the case that the three prime gene breakpoint
		 * is in an intergenic region
		 * */
		if(
				three_prime_transcripts.size()==0 &&
				tp_gene_fp_UTR_transcripts.size()==0 &&
				tp_gene_tp_UTR_transcripts.size()==0
		)
		{
			if(five_prime_transcripts.size() != 0)
			{
				/*
				 * this is 5G-CDS / 3G-Intergenic breakpoint case
				 * Only the 5G gene is mantained
				 * */
				for (String five_prime_t : five_prime_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
					fusion_sequence_info.SetThreePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetType5P("Coding");
					fusion_sequence_info.SetType3P("Intergenic");
					
					
					fusion_sequence_info.SetFive_NucleotideSeq(getFivePrimeSequenceFromTranscriptId(fusion_sequence_info.GetFivePrimeTranscritp_Id(), fiveprimegene, gtfi, indexedfasta));

					fusion_sequence_info.SetThree_NucleotideSeq("");

					fusion_sequence_info.SetCompleteThreeprime_NucleotideSeq("");
					
					if(fusion_sequence_info.GetFive_NucleotideSeq().compareTo("")!=0)
					{
						fusion_sequence_info.SetFrame("UnSpec");
						
						fusion_sequence_info.SetExomeNumber1(gtfi.getExonNumberByTranscriptIdAndBreakPoint(five_prime_t, fiveprimegene.getGeneBreakpoint()));
						fusion_sequence_info.SetExomeNumber2(-1);
							
						String fiveprime_aminoacid = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq());

						fusion_sequence_info.SetFivePrimeAminoAcidStart(1);
						fusion_sequence_info.SetFivePrimeAminoAcidStop(fiveprime_aminoacid.length());

						fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
						fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
						
						String chimeric_aminoacid_sequence = AminoAcid.translate(fusion_sequence_info.GetFive_NucleotideSeq() + fusion_sequence_info.GetThree_NucleotideSeq());
						if(chimeric_aminoacid_sequence.contains("X"))
						{
							fusion_sequence_info.SetFrame("FrameShift");
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X", "\\*"));
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence.replaceAll("X[A-Z]*", "*"));	
							fusion_sequence_info.SetThreePrimeAminoAcidStart(-1);
							fusion_sequence_info.SetThreePrimeAminoAcidStop(-1);
						}
						else
						{
							fusion_sequence_info.SetChimericAminoAcidSeq(chimeric_aminoacid_sequence);
						}

						fusion_sequences_info.add(fusion_sequence_info);
						
					}
				}
			}
			else if(fp_gene_tp_UTR_transcripts.size()!=0)
			{
				/*
				 * this is 5G-3UTR / 3G-Intergenic breakpoint case
				 * Only the 5G gene is mantained up to the 3UTR, the intergenic region does make difference
				 * because a stop codon is introduced after the 3UTR
				 * 
				 * This case has higher priority compared to the 5G-5UTR
				 * */
				for (String five_prime_t : fp_gene_tp_UTR_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
					fusion_sequence_info.SetThreePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetFrame("UnSpec");
					fusion_sequence_info.SetType5P("3UTR");
					fusion_sequence_info.SetType3P("Intergenic");
					
					fusion_sequences_info.add(fusion_sequence_info);
				}

			}
			else if(fp_gene_fp_UTR_transcripts.size()!=0)
			{
				/*
				 * this is 5G-5UTR / 3G-Intergenic breakpoint case
				 * nothing is actually produced as coding sequence
				 * */
				for (String five_prime_t : fp_gene_fp_UTR_transcripts)
				{
					FusionSequenceInfo fusion_sequence_info = new FusionSequenceInfo();
					
					fusion_sequence_info.SetFivePrimeTranscritp_Id(five_prime_t);
					fusion_sequence_info.SetThreePrimeTranscritp_Id("UnSpec");
					fusion_sequence_info.SetFrame("UnSpec");
					fusion_sequence_info.SetType5P("5UTR");
					fusion_sequence_info.SetType3P("Intergenic");
					
					fusion_sequences_info.add(fusion_sequence_info);
				}

			}
		}
		
		
		
		
		
		
		
		
	}
	
	private String getTpGeneTpFUTRposSequenceFromTranscriptId(
			String threePrimeT, GeneInfo threeprimegene, GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta) {
		// TODO Auto-generated method stub
		return null;
	}

	public ArrayList<FusionSequenceInfo> getFusionSequencesInfo()
	{
		return fusion_sequences_info;
	}
	
	private String FindATG_sequence(String seq)
	{
		String subseq = "";
		
		int index_ATG = seq.toUpperCase().indexOf("ATG");
		
		subseq = seq.substring(index_ATG);
		
		return subseq;
	}
	
	//this method return true if the fusion breakpoint of the transcritp id t_id
	//is within start and stop codon of the CDS sequence
	private boolean isCodingFusionSequence(
			String t_id,
			GeneInfo g,
			GTFInterface gtfi
	)
	{
		
		Vector<Location> locationList_startCodon = gtfi.getStartCodonLocationListByTranscriptId(t_id);
		Vector<Location> locationList_stopCodon = gtfi.getStopCodonLocationListByTranscriptId(t_id);
		
		for (Location l_start : locationList_startCodon)
		{
			for(Location l_stop : locationList_stopCodon)
			{
				if(g.getGeneStrand().compareTo("+")==0)
				{
					if(g.getGeneBreakpoint() >= l_start.start() && g.getGeneBreakpoint() <= l_stop.end())
					{
						return true;
					}
				}
				if(g.getGeneStrand().compareTo("-")==0)
				{
					if(g.getGeneBreakpoint() <= l_start.end() && g.getGeneBreakpoint() >= l_stop.start())
					{
						return true;
					}
				}
			}
		}
		
		return false;
	}
	
	//this method return true if the fusion breakpoint of the transcritp id t_id
	//is within the 5' UTR of the gene
	private boolean is5UtrFusionSequence(
			String t_id,
			GeneInfo g,
			GTFInterface gtfi
	)
	{
		if(g.getGeneStrand().compareTo("+")==0)
		{
			Vector<Location> locationList_startCodon = gtfi.getStartCodonLocationListByTranscriptId(t_id);
			
			Vector<Location> locationList_allTranscripts = gtfi.getLocationListByTranscriptId(t_id);
			
			for (Location l_startcodon : locationList_startCodon)
			{
				//carefull: in the second condition has to ben > not >=
				if(locationList_allTranscripts.get(0).start() <= g.getGeneBreakpoint() &&  l_startcodon.start() > g.getGeneBreakpoint() )
				{
					// Breakpoint in the UTR
					return true;
				}
			}
		}
		
		if(g.getGeneStrand().compareTo("-")==0)
		{
			Vector<Location> locationList_startCodon = gtfi.getStartCodonLocationListByTranscriptId(t_id);
			
			Vector<Location> locationList_allTranscripts = gtfi.getLocationListByTranscriptId(t_id);
			
			for (Location l_startcodon : locationList_startCodon)
			{
				//carefull: in the second condition has to ben < not <=
				if(locationList_allTranscripts.get(0).end() >= g.getGeneBreakpoint() &&  l_startcodon.end() < g.getGeneBreakpoint() )
				{
					// Breakpoint in the UTR
					return true;
				}
			}
		}

		return false;
	}
	
	//this method return true if the fusion breakpoint of the transcritp id t_id
	//is within the 3' UTR of the gene
	private boolean is3UtrFusionSequence(
			String t_id,
			GeneInfo g,
			GTFInterface gtfi
	)
	{
		if(g.getGeneStrand().compareTo("+")==0)
		{
			Vector<Location> locationList_stopCodon = gtfi.getStopCodonLocationListByTranscriptId(t_id);
			
			Vector<Location> locationList_allTranscripts = gtfi.getLocationListByTranscriptId(t_id);
			Location first_loc = GetLastLocationFromList(locationList_allTranscripts);
			for (Location l_stopCodon : locationList_stopCodon)
			{
				//carefull: in the second condition has to ben < not <=
				if(first_loc.end() >= g.getGeneBreakpoint() &&  l_stopCodon.end() < g.getGeneBreakpoint() )
				{
					// Breakpoint in the UTR
					return true;
				}
			}
		}
		
		if(g.getGeneStrand().compareTo("-")==0)
		{
			Vector<Location> locationList_stopCodon = gtfi.getStopCodonLocationListByTranscriptId(t_id);
			
			Vector<Location> locationList_allTranscripts = gtfi.getLocationListByTranscriptId(t_id);
			Location first_loc = GetFirstLocationFromList(locationList_allTranscripts);
			for (Location l_stopcodon : locationList_stopCodon)
			{
				//carefull: in the second condition has to ben < not <=
				if(first_loc.start() <= g.getGeneBreakpoint() &&  l_stopcodon.start() > g.getGeneBreakpoint() )
				{
					// Breakpoint in the UTR
					return true;
				}
			}
		}

		return false;
	}
	
	private Location GetFirstLocationFromList (Vector<Location> list)
	{
		Location l_res = list.firstElement();
		for(Location l:list)
		{
			if(l.isBefore(l_res))
			{
				l_res=l;
			}
		}
		return l_res;
	}
	
	private Location GetLastLocationFromList (Vector<Location> list)
	{
		Location l_res = list.firstElement();
		for(Location l:list)
		{
			if(l.isAfter(l_res))
			{
				l_res=l;
			}
		}
		return l_res;
	}
	
	private String getFpGeneFpFUTRposSequenceFromTranscriptId(
			String t_id,
			GeneInfo fiveprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String fiveprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Location lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), fiveprimegene.getGeneBreakpoint());

		for (Location l : locationList)
		{
			if(l.start()<=fiveprimegene.getGeneBreakpoint())
			{
				if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
				{
					lastlocation = new Location(l.start(),fiveprimegene.getGeneBreakpoint());
					overlap_breakpoint=true;
					sb.append(getSequenceStringFromLocation(
							indexedfasta, 
							fiveprime_ref, 
							fiveprimegene.getGeneStrand(), 
							l.start()+1, 
							fiveprimegene.getGeneBreakpoint()));
				}
				else
				{
					lastlocation = new Location(l);
					sb.append(getSequenceStringFromLocation(
							indexedfasta, 
							fiveprime_ref, 
							fiveprimegene.getGeneStrand(), 
							l.start()+1, 
							l.end()));			
				}
			}	
		}
		if(lastlocation.end()<fiveprimegene.getGeneBreakpoint())
		{
			String intronsequence = getSequenceStringFromLocation(
					indexedfasta, 
					fiveprime_ref, 
					fiveprimegene.getGeneStrand(), 
					lastlocation.end()+1, 
					fiveprimegene.getGeneBreakpoint());
			sb.append(intronsequence);
			overlap_breakpoint=true;
		}


		
		return (overlap_breakpoint) ? sb.toString() : "";
	}

	private String getFpGeneFpFUTRnegSequenceFromTranscriptId(
			String t_id,
			GeneInfo fiveprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String fiveprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Location lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), fiveprimegene.getGeneBreakpoint());

		for (Location l : locationList)
		{
			
			if(l.end()>=fiveprimegene.getGeneBreakpoint())
			{
				if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
				{
					lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), l.end());
					overlap_breakpoint=true;
					sb.append(getSequenceStringFromLocation(
							indexedfasta, 
							fiveprime_ref, 
							fiveprimegene.getGeneStrand(), 
							fiveprimegene.getGeneBreakpoint(), 
							l.end()));
				}
				else
				{
					lastlocation = new Location(l);
					sb.append(getSequenceStringFromLocation(
							indexedfasta, 
							fiveprime_ref, 
							fiveprimegene.getGeneStrand(), 
							l.start()+1, 
							l.end()));			
				}
			}	
		}
		if(lastlocation.start()>fiveprimegene.getGeneBreakpoint())
		{
			String intronsequence = getSequenceStringFromLocation(
					indexedfasta, 
					fiveprime_ref, 
					fiveprimegene.getGeneStrand(), 
					fiveprimegene.getGeneBreakpoint(), 
					lastlocation.start());
			sb.append(intronsequence);
			overlap_breakpoint=true;
		}			

		return (overlap_breakpoint) ? sb.toString() : "";
	}
	
	private String getTpGeneFpFUTRposSequenceFromTranscriptId(
			String t_id,
			GeneInfo threeprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String threeprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t_id);
		Vector<Location> locationList_startcodon = gtfi.getStartCodonLocationListByTranscriptId(t_id);
		
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Vector<Location> locationSequence = new Vector<Location>();
		
		for (Location l_startcodon : locationList_startcodon)
		{
			for (Location l : locationList)
			{
				
				if(l.end()>=threeprimegene.getGeneBreakpoint() && l.start() <= l_startcodon.start())
				{
					// the breakpoint is not in the same exon of the start_codon
					// the location is not in the same exon of the start_codon
					// the location is the the same exon of the breakpoint
					if(		l.end()>=threeprimegene.getGeneBreakpoint() && 
							l.start()<=threeprimegene.getGeneBreakpoint() && 
							l.end() <= l_startcodon.start())
					{
						Location location = new Location(threeprimegene.getGeneBreakpoint(), l.end());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								threeprimegene.getGeneBreakpoint(), 
								l.end()));
					}
					// the breakpoint is not in the same exon of the start_codon
					// the location is in a different exon of the start_codon and breakpoint
					// (an exon in the middle)
					else if(l.start()>threeprimegene.getGeneBreakpoint() && 
							l.end()<l_startcodon.end())
					{
						Location location = new Location(l);
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));					
					}
					// breakpoint and startcodon are in different locations
					// the location overlap the startcodon but not the breakpoint
					else if(l.start()>=threeprimegene.getGeneBreakpoint() && 
							l.start()<l_startcodon.start() &&
							l.end()>=l_startcodon.start())
					{
						Location location = new Location(l.start(), l_startcodon.start());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l_startcodon.start()));			
					}
					// the startcodon and the breakpoint are on the same location
					// the location overlaps both
					else if(l.end()>=l_startcodon.start() && 
							l.start()<=threeprimegene.getGeneBreakpoint())
					{
						Location location = new Location(threeprimegene.getGeneBreakpoint(), l_startcodon.end());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								threeprimegene.getGeneBreakpoint(), 
								l_startcodon.start()));			
					}

					
					
					
					
					
//					if(		l.start()<=threeprimegene.getGeneBreakpoint() && 
//							l.end()>=threeprimegene.getGeneBreakpoint() && 
//							l.end() <= l_startcodon.start())
//					{
//						Location location = new Location(threeprimegene.getGeneBreakpoint(), l.end());
//						locationSequence.add(location);
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								threeprimegene.getGeneBreakpoint(), 
//								l.end()));
//					}
//					else if(l.start()<=threeprimegene.getGeneBreakpoint() && 
//							l.end()>=threeprimegene.getGeneBreakpoint() && 
//							l.end() >= l_startcodon.start())
//					{
//						Location location = new Location(threeprimegene.getGeneBreakpoint(), l_startcodon.start());
//						locationSequence.add(location);
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								threeprimegene.getGeneBreakpoint(), 
//								l_startcodon.start()));					
//					}
//					else if(l.start()>=threeprimegene.getGeneBreakpoint() && 
//							l.end()<= l_startcodon.start())
//					{
//						Location location = new Location(threeprimegene.getGeneBreakpoint(), l.end());
//						locationSequence.add(location);
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));				
//					}
					
					
//					else
//					{
//						Location location = new Location(l);
//						locationSequence.add(location);
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
				}	
			}
			Location l_firsexon = locationSequence.firstElement();
			if(l_firsexon.start()>threeprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						threeprime_ref, 
						threeprimegene.getGeneStrand(), 
						threeprimegene.getGeneBreakpoint(), 
						l_firsexon.start());
				String tmp = sb.toString();
				sb = new StringBuffer();
				sb.append(intronsequence);
				sb.append(tmp);
				overlap_breakpoint=true;
			}				
		}
		

		return (overlap_breakpoint) ? sb.toString() : "";
	}
	
	private String getTpGeneFpFUTRnegSequenceFromTranscriptId(
			String t_id,
			GeneInfo threeprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String threeprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t_id);
		Vector<Location> locationList_startcodon = gtfi.getStartCodonLocationListByTranscriptId(t_id);
		
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Vector<Location> locationSequence = new Vector<Location>();
		
		for (Location l_startcodon : locationList_startcodon)
		{
			for (Location l : locationList)
			{
				
				if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end() > l_startcodon.end())
				{
					// the breakpoint is not in the same exon of the start_codon
					// the location is not in the same exon of the start_codon
					// the location is the the same exon of the breakpoint
					if(		l.start()<=threeprimegene.getGeneBreakpoint() && 
							l.end()>=threeprimegene.getGeneBreakpoint() && 
							l.start() >= l_startcodon.end())
					{
						Location location = new Location(l.start()+1, threeprimegene.getGeneBreakpoint());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								threeprimegene.getGeneBreakpoint()));
					}
					// the breakpoint is not in the same exon of the start_codon
					// the location is in a different exon of the start_codon and breakpoint
					// (an exon in the middle)
					else if(l.start()>l_startcodon.end() && 
							l.end()<threeprimegene.getGeneBreakpoint())
					{
						Location location = new Location(l);
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));					
					}
					// breakpoint and startcodon are in different locations
					// the location overlap the startcodon but not the breakpoint
					else if(l.start()<=l_startcodon.start() && 
							l.end()>=l_startcodon.start() &&
							l.end()<=threeprimegene.getGeneBreakpoint())
					{
						Location location = new Location(l_startcodon.end(), l.end());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l_startcodon.end()+1, 
								l.end()));			
					}
					// the startcodon and the breakpoint are on the same location
					// the location overlaps both
					else if(l.start()<=l_startcodon.start() && 
							l.end()>=threeprimegene.getGeneBreakpoint())
					{
						Location location = new Location(l_startcodon.end(), threeprimegene.getGeneBreakpoint());
						locationSequence.add(location);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l_startcodon.end()+1, 
								threeprimegene.getGeneBreakpoint()));			
					}
//					else
//					{
//						Location location = new Location(l);
//						locationSequence.add(location);
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
				}	
			}
			Location l_firsexon = locationSequence.firstElement();
			if(l_firsexon.end()<threeprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						threeprime_ref, 
						threeprimegene.getGeneStrand(), 
						l_firsexon.end()+1, 
						threeprimegene.getGeneBreakpoint());
				String tmp = sb.toString();
				sb = new StringBuffer();
				sb.append(intronsequence);
				sb.append(tmp);
				overlap_breakpoint=true;
			}				
		}
		

		return (overlap_breakpoint) ? sb.toString() : "";
	}
	

	
	

	
	private String getEntireSequenceFromTranscript(
			String t_id,
			GeneInfo geneinfo,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta					
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getCDSLocationListByTranscriptId(t_id);

		for (Location l : locationList)
		{
			sb.append(getSequenceStringFromLocation(
					indexedfasta, 
					ref, 
					geneinfo.getGeneStrand(), 
					l.start()+1, 
					l.end()));

		}
		
		return sb.toString();
	}
	
	// the five prime onlycds is different than three prime. it does not consider the introns.
	// so it will take the usual gene sequence. usually in case of intron in CSD we have to consider the sequence up to the breakpoint
	// if the breakpoint is in the 3UTR of the 5pG you want simply take the entire gene sequence without the sequence up to the breakpoint.
	private String getFivePrimeSequenceFromTranscriptId_ONLYCDS(
			String t_id,
			GeneInfo fiveprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String fiveprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getCDSLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Vector<Location> locationSequence = new Vector<Location>();
		Location lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), fiveprimegene.getGeneBreakpoint());
		//5' positive strand
		if(fiveprimegene.getGeneStrand().compareTo("+")==0)
		{

			for (Location l : locationList)
			{
				if(l.end()<=fiveprimegene.getGeneBreakpoint())
				{
					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
					{
						lastlocation = new Location(l.start(),fiveprimegene.getGeneBreakpoint());
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								fiveprimegene.getGeneBreakpoint()));
					}
					else
					{
						lastlocation = new Location(l);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}	
			}
//			if(lastlocation.end()<fiveprimegene.getGeneBreakpoint())
//			{
//				String intronsequence = getSequenceStringFromLocation(
//						indexedfasta, 
//						fiveprime_ref, 
//						fiveprimegene.getGeneStrand(), 
//						lastlocation.end()+1, 
//						fiveprimegene.getGeneBreakpoint());
//				sb.append(intronsequence);
//				overlap_breakpoint=true;
//			}
		}
		else //5' negative strand
		{
			for (Location l : locationList)
			{
				
				if(l.end()>=fiveprimegene.getGeneBreakpoint())
				{
					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
					{
						lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), l.end());
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								fiveprimegene.getGeneBreakpoint(), 
								l.end()));
					}
					else
					{
						lastlocation = new Location(l);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}	
			}
//			if(lastlocation.start()>fiveprimegene.getGeneBreakpoint())
//			{
//				String intronsequence = getSequenceStringFromLocation(
//						indexedfasta, 
//						fiveprime_ref, 
//						fiveprimegene.getGeneStrand(), 
//						fiveprimegene.getGeneBreakpoint(), 
//						lastlocation.start());
//				sb.append(intronsequence);
//				overlap_breakpoint=true;
//			}			
		}
		return sb.toString();
	}
	
	private String getFivePrimeSequenceFromTranscriptId(
			String t_id,
			GeneInfo fiveprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String fiveprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getCDSLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Vector<Location> locationSequence = new Vector<Location>();
		
//		//5' positive strand
//		if(fiveprimegene.getGeneStrand().compareTo("+")==0)
//		{
//			for (Location l : locationList)
//			{
//				if(l.end()<=fiveprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
//					{
//						Location l_seq = new Location(l.start(), fiveprimegene.getGeneBreakpoint());
//						locationSequence.add(l_seq);						
//						overlap_breakpoint=true;
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								fiveprime_ref, 
////								fiveprimegene.getGeneStrand(), 
////								l.start()+1, 
////								fiveprimegene.getGeneBreakpoint()));
//					}
//					else
//					{
//						Location l_seq = new Location(l);
//						locationSequence.add(l_seq);
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								fiveprime_ref, 
////								fiveprimegene.getGeneStrand(), 
////								l.start()+1, 
////								l.end()));			
//					}
//				}	
//			}
//			if(locationSequence.lastElement().end() < fiveprimegene.getGeneBreakpoint())
//			{
//				// it means that the 5' breakpoint is in an Intron
//				Location l_last_old = locationSequence.lastElement();
//				Location l_last_new = new Location(l_last_old.start(), fiveprimegene.getGeneBreakpoint());
//				locationSequence.removeElementAt(locationSequence.size()-1);
//				locationSequence.add(l_last_new);
//				overlap_breakpoint=true;
//			}
//			for (Location l_seq : locationSequence)
//			{
//				sb.append(getSequenceStringFromLocation(
//						indexedfasta, 
//						fiveprime_ref, 
//						fiveprimegene.getGeneStrand(), 
//						l_seq.start()+1, 
//						l_seq.end()));		
//			}
//		}
//		else //5' negative strand
//		{
//			for (Location l : locationList)
//			{
//				
//				if(l.end()>=fiveprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
//					{						
//						Location l_seq = new Location(fiveprimegene.getGeneBreakpoint(), l.end());
//						locationSequence.add(l_seq);						
//						overlap_breakpoint=true;
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								fiveprime_ref, 
////								fiveprimegene.getGeneStrand(), 
////								fiveprimegene.getGeneBreakpoint(), 
////								l.end()));
//					}
//					else
//					{
//						Location l_seq = new Location(l);
//						locationSequence.add(l_seq);
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								fiveprime_ref, 
////								fiveprimegene.getGeneStrand(), 
////								l.start()+1, 
////								l.end()));			
//					}
//				}	
//			}
//			if(locationSequence.lastElement().start() > fiveprimegene.getGeneBreakpoint())
//			{
//				// it means that the 5' breakpoint is in an Intron
//				Location l_last_old = locationSequence.lastElement();
//				Location l_last_new = new Location(fiveprimegene.getGeneBreakpoint(), l_last_old.end());
//				locationSequence.removeElementAt(locationSequence.size()-1);
//				locationSequence.add(l_last_new);
//				overlap_breakpoint=true;
//			}
//			for (Location l_seq : locationSequence)
//			{
//				sb.append(getSequenceStringFromLocation(
//						indexedfasta, 
//						fiveprime_ref, 
//						fiveprimegene.getGeneStrand(), 
//						l_seq.start()+1, 
//						l_seq.end()));
//			}
//		}

		Location lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), fiveprimegene.getGeneBreakpoint());
		//5' positive strand
		if(fiveprimegene.getGeneStrand().compareTo("+")==0)
		{

			for (Location l : locationList)
			{
				if(l.end()<=fiveprimegene.getGeneBreakpoint())
				{
					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
					{
						lastlocation = new Location(l.start(),fiveprimegene.getGeneBreakpoint());
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								fiveprimegene.getGeneBreakpoint()));
					}
					else
					{
						lastlocation = new Location(l);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}	
			}
			if(lastlocation.end()<fiveprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						fiveprime_ref, 
						fiveprimegene.getGeneStrand(), 
						lastlocation.end()+1, 
						fiveprimegene.getGeneBreakpoint());
				sb.append(intronsequence);
				overlap_breakpoint=true;
			}
		}
		else //5' negative strand
		{
			for (Location l : locationList)
			{
				
				if(l.end()>=fiveprimegene.getGeneBreakpoint())
				{
					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
					{
						lastlocation = new Location(fiveprimegene.getGeneBreakpoint(), l.end());
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								fiveprimegene.getGeneBreakpoint(), 
								l.end()));
					}
					else
					{
						lastlocation = new Location(l);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								fiveprime_ref, 
								fiveprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}	
			}
			if(lastlocation.start()>fiveprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						fiveprime_ref, 
						fiveprimegene.getGeneStrand(), 
						fiveprimegene.getGeneBreakpoint(), 
						lastlocation.start());
				sb.append(intronsequence);
				overlap_breakpoint=true;
			}			
		}
		
		
		
//		//5' positive strand
//		if(fiveprimegene.getGeneStrand().compareTo("+")==0)
//		{
//
//			for (Location l : locationList)
//			{
//				if(l.end()<=fiveprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
//					{
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								fiveprime_ref, 
//								fiveprimegene.getGeneStrand(), 
//								l.start()+1, 
//								fiveprimegene.getGeneBreakpoint()));
//					}
//					else
//					{
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								fiveprime_ref, 
//								fiveprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
//				}	
//			}
//		}
//		else //5' negative strand
//		{
//			for (Location l : locationList)
//			{
//				
//				if(l.end()>=fiveprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
//					{
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								fiveprime_ref, 
//								fiveprimegene.getGeneStrand(), 
//								fiveprimegene.getGeneBreakpoint(), 
//								l.end()));
//					}
//					else
//					{
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								fiveprime_ref, 
//								fiveprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
//				}	
//			}
//		}
		
		return (overlap_breakpoint) ? sb.toString() : "";
	}
	
	private String getThreePrimeSequenceFromTranscriptId(
			String t_id,
			GeneInfo threeprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String threeprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getCDSLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		boolean overlap_breakpoint = false; 
		
		Vector<Location> locationSequence = new Vector<Location>();
		
//		//3' positive strand
//		if(threeprimegene.getGeneStrand().compareTo("+")==0)
//		{
//			for (Location l : locationList)
//			{
//				if(l.end()>=threeprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
//					{
//						Location l_seq = new Location(threeprimegene.getGeneBreakpoint()-1, l.end());
//						locationSequence.add(l_seq);						
//						overlap_breakpoint=true;
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								threeprime_ref, 
////								threeprimegene.getGeneStrand(), 
////								threeprimegene.getGeneBreakpoint(), 
////								l.end()));
//					}
//					else
//					{
//						Location l_seq = new Location(l);
//						locationSequence.add(l_seq);
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								threeprime_ref, 
////								threeprimegene.getGeneStrand(), 
////								l.start()+1, 
////								l.end()));			
//					}
//				}
//			}
//			if(locationSequence.firstElement().start() > threeprimegene.getGeneBreakpoint())
//			{
//				// it means that the 5' breakpoint is in an Intron
//				Location l_first_old = locationSequence.firstElement();
//				Location l_first_new = new Location(l_first_old.start(), threeprimegene.getGeneBreakpoint());
//				locationSequence.removeElementAt(0);
//				locationSequence.add(0, l_first_new);
//				overlap_breakpoint=true;
//			}
//			for (Location l_seq : locationSequence)
//			{
//				sb.append(getSequenceStringFromLocation(
//						indexedfasta, 
//						threeprime_ref, 
//						threeprimegene.getGeneStrand(), 
//						l_seq.start()+1, 
//						l_seq.end()));		
//			}
//		}
//		else //3' negative strand
//		{
//			for (Location l : locationList)
//			{
//				if(l.end()<=threeprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
//					{
//						Location l_seq = new Location(l.start()-1, threeprimegene.getGeneBreakpoint());
//						locationSequence.add(l_seq);						
//						overlap_breakpoint=true;
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								threeprime_ref, 
////								threeprimegene.getGeneStrand(), 
////								l.start()+1, 
////								threeprimegene.getGeneBreakpoint()));			
//					}
//					else
//					{
//						Location l_seq = new Location(l);
//						locationSequence.add(l_seq);
////						sb.append(getSequenceStringFromLocation(
////								indexedfasta, 
////								threeprime_ref, 
////								threeprimegene.getGeneStrand(), 
////								l.start()+1, 
////								l.end()));			
//					}
//				}
//			}
//			if(locationSequence.firstElement().end() < threeprimegene.getGeneBreakpoint())
//			{
//				// it means that the 5' breakpoint is in an Intron
//				Location l_first_old = locationSequence.firstElement();
//				Location l_first_new = new Location(l_first_old.start(), threeprimegene.getGeneBreakpoint());
//				locationSequence.removeElementAt(0);
//				locationSequence.add(0, l_first_new);
//				overlap_breakpoint=true;
//			}
//			for (Location l_seq : locationSequence)
//			{
//				sb.append(getSequenceStringFromLocation(
//						indexedfasta, 
//						threeprime_ref, 
//						threeprimegene.getGeneStrand(), 
//						l_seq.start()+1, 
//						l_seq.end()));		
//			}			
//		}
		
		//3' positive strand
		if(threeprimegene.getGeneStrand().compareTo("+")==0)
		{
			for (Location l : locationList)
			{
				if(l.end()>=threeprimegene.getGeneBreakpoint())
				{
					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
					{
						Location l_seq = new Location(threeprimegene.getGeneBreakpoint(), l.end());
						locationSequence.add(l_seq);
						overlap_breakpoint=true;
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								threeprimegene.getGeneBreakpoint(), 
								l.end()));
					}
					else
					{
						Location l_seq = new Location(l);
						locationSequence.add(l_seq);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}
			}
			Location l_firsexon = locationSequence.firstElement();
			if(l_firsexon.start()>threeprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						threeprime_ref, 
						threeprimegene.getGeneStrand(), 
						threeprimegene.getGeneBreakpoint(), 
						l_firsexon.end());
				String tmp = sb.toString();
				sb = new StringBuffer();
				sb.append(intronsequence);
				sb.append(tmp);
				overlap_breakpoint=true;
			}			
			
		}
		else //3' negative strand
		{
			for (Location l : locationList)
			{
				if(l.end()<=threeprimegene.getGeneBreakpoint())
				{
					overlap_breakpoint=true;
					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
					{
						Location l_seq = new Location(l.start(), threeprimegene.getGeneBreakpoint());
						locationSequence.add(l_seq);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								threeprimegene.getGeneBreakpoint()));			
					}
					else
					{
						Location l_seq = new Location(l);
						locationSequence.add(l_seq);
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}
			}
			Location l_firsexon = locationSequence.firstElement();
			if(l_firsexon.end()<threeprimegene.getGeneBreakpoint())
			{
				String intronsequence = getSequenceStringFromLocation(
						indexedfasta, 
						threeprime_ref, 
						threeprimegene.getGeneStrand(), 
						l_firsexon.end()+1, 
						threeprimegene.getGeneBreakpoint());
				String tmp = sb.toString();
				sb = new StringBuffer();
				sb.append(intronsequence);
				sb.append(tmp);
				overlap_breakpoint=true;
			}			

		}
		
		
		
//		//3' positive strand
//		if(threeprimegene.getGeneStrand().compareTo("+")==0)
//		{
//			for (Location l : locationList)
//			{
//				if(l.end()>=threeprimegene.getGeneBreakpoint())
//				{
//					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
//					{
//						overlap_breakpoint=true;
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								threeprimegene.getGeneBreakpoint(), 
//								l.end()));
//					}
//					else
//					{
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
//				}
//			}
//		}
//		else //3' negative strand
//		{
//			for (Location l : locationList)
//			{
//				if(l.end()<=threeprimegene.getGeneBreakpoint())
//				{
//					overlap_breakpoint=true;
//					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
//					{
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								threeprimegene.getGeneBreakpoint()));			
//					}
//					else
//					{
//						sb.append(getSequenceStringFromLocation(
//								indexedfasta, 
//								threeprime_ref, 
//								threeprimegene.getGeneStrand(), 
//								l.start()+1, 
//								l.end()));			
//					}
//				}
//			}
//		}
		
		return (overlap_breakpoint) ? sb.toString() : "";

	}
	
	private String getThreePrimeSequenceFromTranscriptId_ONLYCDS(
			String t_id,
			GeneInfo threeprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta		
	)
	{
		StringBuffer sb = new StringBuffer();
		
		String threeprime_ref = gtfi.getRefByTranscriptID(t_id);
		
		Vector<Location> locationList = gtfi.getCDSLocationListByTranscriptId(t_id);
		// if it does not overlap bp means that is no good transcript
		
		Vector<Location> locationSequence = new Vector<Location>();
		
		//3' positive strand
		if(threeprimegene.getGeneStrand().compareTo("+")==0)
		{
			for (Location l : locationList)
			{
				if(l.end()>=threeprimegene.getGeneBreakpoint())
				{
					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
					{
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								threeprimegene.getGeneBreakpoint(), 
								l.end()));
					}
					else
					{
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}
			}
		}
		else //3' negative strand
		{
			for (Location l : locationList)
			{
				if(l.end()<=threeprimegene.getGeneBreakpoint())
				{
					if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
					{
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								threeprimegene.getGeneBreakpoint()));			
					}
					else
					{
						sb.append(getSequenceStringFromLocation(
								indexedfasta, 
								threeprime_ref, 
								threeprimegene.getGeneStrand(), 
								l.start()+1, 
								l.end()));			
					}
				}
			}
		}
		
		return sb.toString();

	}
	
	
	private void PrintFivePrimeSequence(
			GeneInfo fiveprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta
	)
	{
		Vector<String> transcripts = gtfi.getTranscriptIdByGeneName(fiveprimegene.getGeneName());
		
		for (String t : transcripts)
		{
			System.out.println(t + ":");
			
			String fiveprime_ref = gtfi.getRefByTranscriptID(t);
			
			Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t);

			//5' positive strand
			if(fiveprimegene.getGeneStrand().compareTo("+")==0)
			{

				for (Location l : locationList)
				{
					
					if(l.end()<=fiveprimegene.getGeneBreakpoint())
					{
						if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									fiveprime_ref, 
									fiveprimegene.getGeneStrand(), 
									l.start()+1, 
									fiveprimegene.getGeneBreakpoint()));
						}
						else
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									fiveprime_ref, 
									fiveprimegene.getGeneStrand(), 
									l.start()+1, 
									l.end()));			
						}
					}	
				}
			}
			else //5' negative strand
			{
				for (Location l : locationList)
				{
					
					if(l.end()>=fiveprimegene.getGeneBreakpoint())
					{
						if(l.start()<=fiveprimegene.getGeneBreakpoint() && l.end()>=fiveprimegene.getGeneBreakpoint())
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									fiveprime_ref, 
									fiveprimegene.getGeneStrand(), 
									fiveprimegene.getGeneBreakpoint()+1, 
									l.end()));
						}
						else
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									fiveprime_ref, 
									fiveprimegene.getGeneStrand(), 
									l.start()+1, 
									l.end()));			
						}
					}	
				}
			}
		}

	}
	
	
	private void PrintThreePrimeSequence(
			GeneInfo threeprimegene,
			GTFInterface gtfi,
			IndexedFastaSequenceFile indexedfasta
	)
	{
		Vector<String> transcripts = gtfi.getTranscriptIdByGeneName(threeprimegene.getGeneName());
		
		for (String t : transcripts)
		{
			System.out.println(t + ":");
			
			String threeprime_ref = gtfi.getRefByTranscriptID(t);
			
			Vector<Location> locationList = gtfi.getLocationListByTranscriptId(t);
			
			//3' positive strand
			if(threeprimegene.getGeneStrand().compareTo("+")==0)
			{
				for (Location l : locationList)
				{
					if(l.end()>=threeprimegene.getGeneBreakpoint())
					{
						if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									threeprime_ref, 
									threeprimegene.getGeneStrand(), 
									threeprimegene.getGeneBreakpoint()+1, 
									l.end()));
						}
						else
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									threeprime_ref, 
									threeprimegene.getGeneStrand(), 
									l.start()+1, 
									l.end()));			
						}
					}
				}
			}
			else //3' negative strand
			{
				for (Location l : locationList)
				{
					if(l.end()<=threeprimegene.getGeneBreakpoint())
					{
						if(l.start()<=threeprimegene.getGeneBreakpoint() && l.end()>=threeprimegene.getGeneBreakpoint())
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									threeprime_ref, 
									threeprimegene.getGeneStrand(), 
									l.start()+1, 
									threeprimegene.getGeneBreakpoint()));			
						}
						else
						{
							System.out.println(getSequenceStringFromLocation(
									indexedfasta, 
									threeprime_ref, 
									threeprimegene.getGeneStrand(), 
									l.start()+1, 
									l.end()));			
						}
					}
				}
			}
		}

	}

	
	
	private String getSequenceStringFromLocation(
			IndexedFastaSequenceFile indexedfasta,
			String ref, 
			String strand, 
			long start, 
			long end
	)
	{
		StringBuffer sb = new StringBuffer();
		
		ReferenceSequence rs1 = indexedfasta.getSubsequenceAt(ref, start, end);
		
		if (strand.compareTo("-")==0)
		{
			sb.append(printReferenceSequenceInFastaFormat(rs1, true));
		}
		else
		{
			sb.append(printReferenceSequenceInFastaFormat(rs1, false));
		}
		
		return sb.toString();
	}
	
	private String printReferenceSequenceInFastaFormat(ReferenceSequence rs, boolean reverse)
	{
		StringBuffer sb = new StringBuffer();
		
		byte[] bases = rs.getBases();
		
		if (reverse)
		{
			SequenceUtil.reverseComplement(bases);
		}

		for (int i=0; i<bases.length;i++)
		{
			//if (i > 0 && i % 60 == 0) sb.append("");
			sb.append((char) bases[i]);
		}
		//sb.append("\n");

		return sb.toString();
		//System.out.print(sb);
	}

	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String GTFfilepath = "";
		String five_prime_gene_info = "";
		String three_prime_gene_info = "";
		String reference_FASTA_file = "";
		boolean print_exons = false;
		boolean verbose = false;

        Getopt g = new Getopt("FusionSequenceFromGTF", args, "evg:f:t:r:");
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

        		case 'g':
        			GTFfilepath = g.getOptarg();                	 
        			break;

            	case 'f':
            		five_prime_gene_info = g.getOptarg();                	 
            		break;

            	case 'r':
            		reference_FASTA_file = g.getOptarg();                	 
            		break;

                 case 't':
                	 three_prime_gene_info = g.getOptarg();
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
        		five_prime_gene_info.compareTo("")==0 || 
        		reference_FASTA_file.compareTo("")==0 || 
        		three_prime_gene_info.compareTo("")==0
        	)
        {
        	System.err.println("Error!\nusage: java -jar FusionSequenceFromGTF" +
        			" -g [GTF file list] " +
        			" -f [five prime gene info] " +
        			" -r [reference FASTA file] " +
        			"-v [Verbose (Optional)]" +
        			"-e [Print Exons (Optional)]" +
        			"-t [three prime gene info]"
        			
        	);
        	System.exit(1);
        }

		try 
		{
			
			GTFInterface gtfinterface = new GTFInterface(GTFfilepath);
			
			FusionSequenceFromGTF fsfgtf = new FusionSequenceFromGTF(gtfinterface, five_prime_gene_info, three_prime_gene_info, reference_FASTA_file, verbose, print_exons);
			FusionSequenceInfo f_max = new FusionSequenceInfo(); 
			int max_five_prime_amino_end = 0;
			int max_three_prime_amino_end = 0;
			boolean all_frameshifted = true;
			ArrayList<FusionSequenceInfo> fusionsequenceinfolist = fsfgtf.getFusionSequencesInfo();
			for (FusionSequenceInfo f : fusionsequenceinfolist)
			{
				if(f.GetFrame().equals("InFrame"))
				{
					all_frameshifted = false;
				}
				if(f.GetFivePrimeAminoAcidStop() > max_five_prime_amino_end)
				{
					max_five_prime_amino_end = f.GetFivePrimeAminoAcidStop();
				}
				if(f.GetThreePrimeAminoAcidStop() > max_three_prime_amino_end)
				{
					max_three_prime_amino_end = f.GetThreePrimeAminoAcidStop();
				}
			}
			
			String results = "";
			//I print the longest transcript
			for (FusionSequenceInfo f : fusionsequenceinfolist)
			{
				if(all_frameshifted)
				{
					if(f.GetFivePrimeAminoAcidStop() == max_five_prime_amino_end &&
					   f.GetThreePrimeAminoAcidStop() == max_three_prime_amino_end )
					{
						if(print_exons)
						{
							//System.out.println(f.toString());
							results = f.toString();
						}
						else
						{
							//System.out.println(f.toStringNoExons());
							results = f.toStringNoExons();
						}
					}
				}
				else
				{
					if(f.GetFivePrimeAminoAcidStop() == max_five_prime_amino_end &&
							   f.GetThreePrimeAminoAcidStop() == max_three_prime_amino_end && f.GetFrame().equals("InFrame"))
					{
						if(print_exons)
						{
							//System.out.println(f.toString());
							results = f.toString();
						}
						else
						{
							//System.out.println(f.toStringNoExons());		
							results = f.toStringNoExons();
						}
					}
				}
			}
			
			System.out.println(results);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
        
	}

}
