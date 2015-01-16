package FusionSequenceFromGTF;

public class GeneInfo {
	private String gene_name;
	private String gene_strand;
	private int gene_start;
	private int gene_end;
	private int gene_breakpoint;
	
	public GeneInfo(
			String gene_name,
			String gene_strand,
			int gene_start,
			int gene_end,
			int gene_breakpoint
	)
	{
		this.gene_name = gene_name;
		this.gene_strand = gene_strand;
		this.gene_start = gene_start;
		this.gene_end = gene_end;
		this.gene_breakpoint = gene_breakpoint;
	}
	
	public GeneInfo(String geneInfoString)
	{
		extractGeneInfoFromString(geneInfoString);
	}
	
	private void extractGeneInfoFromString(String geneInfoString)
	{
		String[] geneInfoFields = geneInfoString.split(":");
		this.gene_name = geneInfoFields[0];
		this.gene_strand = geneInfoFields[1];
		this.gene_start = Integer.parseInt(geneInfoFields[2]);
		this.gene_end = Integer.parseInt(geneInfoFields[3]);
		this.gene_breakpoint = Integer.parseInt(geneInfoFields[4]);
	}

	public String getGeneName()
	{
		return this.gene_name;
	}
	
	public String getGeneStrand()
	{
		return this.gene_strand;
	}
	
	public int getGeneStart()
	{
		return this.gene_start;
	}
	
	public int getGeneEnd()
	{
		return this.gene_end;
	}
	
	public int getGeneBreakpoint()
	{
		return this.gene_breakpoint;
	}
	
	public String toString()
	{
		return 
		this.gene_name + "\t" +
		this.gene_strand + "\t" +
		this.gene_start + "\t" +
		this.gene_end + "\t" +
		this.gene_breakpoint;
	}
}
