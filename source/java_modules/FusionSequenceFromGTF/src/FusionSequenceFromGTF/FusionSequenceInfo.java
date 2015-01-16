package FusionSequenceFromGTF;

public class FusionSequenceInfo {
	private String five_prime_transcript_id = "";
	private String three_prime_transcript_id = "";
	private String five_seq = "";
	private String three_seq = "";
	private String complete_threeprimeseq = "";
	private String frame_type = "";
	private int five_prime_aminoacid_start = 0;
	private int five_prime_aminoacid_stop = 0;
	private int three_prime_aminoacid_start = 0;
	private int three_prime_aminoacid_stop = 0;
	private String chimeric_aminoacid_sequence;
	private int exon_number1 = 0;
	private int exon_number2 = 0;
	private String type5p = "";
	private String type3p = "";
	
	public void SetType5P(String type5p)
	{
		this.type5p = type5p;
	}
	
	public void SetType3P(String type3p)
	{
		this.type3p = type3p;
	}
	
	public String GetType5P()
	{
		return this.type5p;
	}
	
	public String GetType3P()
	{
		return this.type3p;
	}
	
	public void SetFivePrimeTranscritp_Id(String five_prime_transcript_id)
	{
		this.five_prime_transcript_id = five_prime_transcript_id;
	}
	
	public void SetThreePrimeTranscritp_Id(String three_prime_transcript_id)
	{
		this.three_prime_transcript_id = three_prime_transcript_id;
	}
	
	public void SetFive_NucleotideSeq(String five_seq)
	{
		this.five_seq = five_seq;
	}
	
	public void SetThree_NucleotideSeq(String three_seq)
	{
		this.three_seq = three_seq;
	}
	
	public void SetCompleteThreeprime_NucleotideSeq(String complete_threeprimeseq)
	{
		this.complete_threeprimeseq = complete_threeprimeseq;
	}
	
	public void SetFrame(String frame_type)
	{
		this.frame_type = frame_type;
	}
	
	public void SetFivePrimeAminoAcidStart(int five_prime_aminoacid_start)
	{
		this.five_prime_aminoacid_start = five_prime_aminoacid_start;
	}
	
	public void SetFivePrimeAminoAcidStop(int five_prime_aminoacid_stop)
	{
		this.five_prime_aminoacid_stop = five_prime_aminoacid_stop;
	}
	
	public void SetThreePrimeAminoAcidStart(int three_prime_aminoacid_start)
	{
		this.three_prime_aminoacid_start = three_prime_aminoacid_start;
	}
	
	public void SetThreePrimeAminoAcidStop(int three_prime_aminoacid_stop)
	{
		this.three_prime_aminoacid_stop = three_prime_aminoacid_stop;
	}
	
	public void SetChimericAminoAcidSeq(String chimeric_aminoacid_sequence)
	{
		this.chimeric_aminoacid_sequence = chimeric_aminoacid_sequence;
	}
	
	public void SetExomeNumber1(int exon_number1)
	{
		this.exon_number1 = exon_number1;
	}
	
	public void SetExomeNumber2(int exon_number2)
	{
		this.exon_number2 = exon_number2;
	}
	
	public String GetFivePrimeTranscritp_Id()
	{
		return this.five_prime_transcript_id;
	}
	
	public String GetThreePrimeTranscritp_Id()
	{
		return this.three_prime_transcript_id;
	}
	
	public String GetFive_NucleotideSeq()
	{
		return this.five_seq;
	}
	
	public String GetThree_NucleotideSeq()
	{
		return this.three_seq;
	}
	
	public String GetCompleteThreeprime_NucleotideSeq()
	{
		return this.complete_threeprimeseq;
	}
	
	public String GetFrame()
	{
		return this.frame_type;
	}
	
	public int GetFivePrimeAminoAcidStart()
	{
		return this.five_prime_aminoacid_start;
	}
	
	public int GetFivePrimeAminoAcidStop()
	{
		return this.five_prime_aminoacid_stop;
	}
	
	public int GetThreePrimeAminoAcidStart()
	{
		return this.three_prime_aminoacid_start;
	}
	
	public int GetThreePrimeAminoAcidStop()
	{
		return this.three_prime_aminoacid_stop;
	}
	
	public String GetChimericAminoAcidSeq()
	{
		return this.chimeric_aminoacid_sequence;
	}
	
	public int GetExomeNumber1()
	{
		return this.exon_number1;
	}
	
	public int GetExomeNumber2()
	{
		return this.exon_number2;
	}
	
	public String toStringNoExons()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(this.five_prime_transcript_id);
		sb.append("\t");
		sb.append(this.three_prime_transcript_id);
		sb.append("\t");
		sb.append(this.frame_type);
		sb.append("\t");
		sb.append(this.five_seq + "|" + this.three_seq);
		sb.append("\t");
		sb.append(this.five_prime_aminoacid_start);
		sb.append("\t");
		sb.append(this.five_prime_aminoacid_stop);
		sb.append("\t");
		sb.append(this.three_prime_aminoacid_start);
		sb.append("\t");
		sb.append(this.three_prime_aminoacid_stop);
		sb.append("\t");
		sb.append(this.chimeric_aminoacid_sequence);
		return sb.toString();
	}

	
	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(this.five_prime_transcript_id);
		sb.append("\t");
		sb.append(this.three_prime_transcript_id);
		sb.append("\t");
		sb.append(this.frame_type);
		sb.append("\t");
		sb.append(this.five_seq + "|" + this.three_seq);
		sb.append("\t");
		sb.append(this.five_prime_aminoacid_start);
		sb.append("\t");
		sb.append(this.five_prime_aminoacid_stop);
		sb.append("\t");
		sb.append(this.three_prime_aminoacid_start);
		sb.append("\t");
		sb.append(this.three_prime_aminoacid_stop);
		sb.append("\t");
		sb.append(this.chimeric_aminoacid_sequence);
		sb.append("\t");
		if(this.exon_number1!=0)
		{
			sb.append(this.exon_number1);
		}
		else
		{
			sb.append("Intron");
		}
		sb.append("\t");
		if(this.exon_number2!=0)
		{
			sb.append(this.exon_number2);
		}
		else
		{
			sb.append("Intron");
		}			
		return sb.toString();
	}
	
	
}
