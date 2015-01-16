package FusionSequenceFromGTF;


import java.util.*;

import javagene.seq.*;
import javagene.io.*;

public class GTFInterface extends Thread{

	private FeatureList gtffeatures;
		
	public GTFInterface(String gtfpath) throws Exception
	{
		gtffeatures = Gff.read(gtfpath);
	}
	
	public FeatureList getGTF()
	{
		return gtffeatures;
	}
	
	public Vector<String> getGeneIdByLocation (String ref, Location l) 
	{
		Vector<String> geneIds = new Vector<String>();
		Collection<String> col = new Vector<String>();
		try
		{
			col = gtffeatures.selectOverlapping(ref, l, true).attributeValues("gene_id");
		}
		catch(ConcurrentModificationException c)
		{
			System.out.println("beccata1");
			System.exit(1);
		}
		catch (Exception e) {
			// TODO: handle exception
			System.out.println("nada");
			System.exit(1);

		}
		
		try
		{
			geneIds = new Vector<String>(col);
		}
		catch(ConcurrentModificationException c)
		{
			System.out.println("beccata2");
			System.exit(1);
		}
		
		return geneIds;
	}

	// this method returns a Vector<String> of geneNames because it could be possible (maybe) that 
	//a gene_id refers to multiple geneName (the contrary is more probable)
	public Vector<String> getGeneNameByGeneId (String geneid) 
	{
		Vector<String> geneNames = new Vector<String>();
		
		Collection<String> col = new Vector<String>();
		
		col = gtffeatures.selectByAttribute("gene_id", geneid).attributeValues("gene_name");
		
		geneNames = new Vector<String>(col);
				
		return geneNames;
	}

	// this method returns a Vector<String> of geneNames because it could be possible (maybe) that 
	// a gene_id refers to multiple geneName (the contrary is more probable)
	// For each gene it also return the coordinates. For ex.: VEGFA:startcoord:endcoord
	//metodo NON COMPLETO!!!!!!
	public Vector<String> getGeneNameAndCoordinatesByGeneId (String geneid) 
	{
		Vector<String> geneNames = new Vector<String>();
		
		Vector<String> geneInfos = new Vector<String>();
		
		Collection<String> col = new Vector<String>();
		
		//col = gtffeatures.selectByAttribute("gene_id", geneid).attributeValues("gene_name");
		
		FeatureList fl = gtffeatures.selectByAttribute("gene_id", geneid);
			
		for (FeatureI f : fl)
		{
			col = fl.attributeValues("gene_name");
			
			geneNames = new Vector<String>(col);
			
			for (String genename:geneNames)
			{
				Location l = f.location();
				
				String geneinfo = genename + ":" + l.start() + ":" + l.end();
				
				geneInfos.add(geneinfo);
			}
		}
		
//		for (String genename:geneNames)
//		{
//			FeatureList loclist = fl.selectByAttribute("gene_name", genename);
//			for (FeatureI f : loclist)
//			{
//				Location l = f.location();
//				String geneinfo = genename + ":" + l.start() + ":" + l.end();
//				geneInfos.add(geneinfo);
//			}
//		}
				
		return geneInfos;
	}

	//this method returns only a single reference because if a gene_id maps on multiple reference something is strange
	public Vector<String> getReferenceNameByGeneId (String geneid) 
	{
		Vector<String> geneReferences = new Vector<String>();
		
		
		FeatureList fl = gtffeatures.selectByAttribute("gene_id", geneid);
		
		for (FeatureI f : fl)
		{
			//add seqname uniquely
			if (geneReferences.isEmpty())
			{				
				geneReferences.add(f.seqname());
			}
			else
			{
				Iterator<String> g_i = geneReferences.iterator();
				boolean found = false;
				while (g_i.hasNext() && !found)
				{
					String g = g_i.next();
					if(g.compareTo(f.seqname())==0)
					{
						found = true;
					}
				}
				if (!found)
				{
					geneReferences.add(f.seqname());
				}
			}
		}
		
		return geneReferences;
	}
	
	//multiple gene name, comma separated as in the case of chimerascan results
	public Vector<String> getTranscriptIdByMultipleGeneName (String multi_g_name) 
	{
		TreeMap<String, String> transcripts_map = new TreeMap<String, String>();
		
		Vector<String> transcript_ids = new Vector<String>();
		
		String [] gene_names = multi_g_name.split(",");
		for (String gene_name : gene_names)
		{	
			FeatureList fl = gtffeatures.selectByAttribute("gene_name", gene_name);
			
			for (FeatureI f : fl)
			{
				//String transcript = f.getAttribute("transcript_id");
				if(f.getAttribute("gene_name").compareTo(gene_name)==0)
				{
					String transcript = f.getAttribute("transcript_id");
					transcripts_map.put(transcript,null);				
				}
			}
			
			for (String t : transcripts_map.keySet())
			{
				transcript_ids.add(t);
			}
		}
		return transcript_ids;
	}
	
	//multiple gene name, comma separated as in the case of chimerascan results
	public int getExonNumberByTranscriptIdAndBreakPoint (String t_id, int breakpoint) 
	{

		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", t_id);
		
		for (FeatureI f : fl)
		{
			if(f.location().isNegative())
			{
				if(f.location().opposite().end() >= breakpoint && f.location().opposite().start() <= breakpoint)
				{
					return Integer.parseInt(f.getAttribute("exon_number"));
				}
			}
			else
			{
				if(f.location().end() >= breakpoint && f.location().start() <= breakpoint)
				{
					return Integer.parseInt(f.getAttribute("exon_number"));
				}
			}
		}
		return 0;
	}
	
	
	
	public Vector<String> getTranscriptIdByGeneName (String gene_name) 
	{
		TreeMap<String, String> transcripts_map = new TreeMap<String, String>();
		
		Vector<String> transcript_ids = new Vector<String>();
		
		FeatureList fl = gtffeatures.selectByAttribute("gene_name", gene_name);
		
		for (FeatureI f : fl)
		{
			//String transcript = f.getAttribute("transcript_id");
			if(f.getAttribute("gene_name").compareTo(gene_name)==0)
			{
				String transcript = f.getAttribute("transcript_id");
				transcripts_map.put(transcript,null);				
			}
		}
		
		for (String t : transcripts_map.keySet())
		{
			transcript_ids.add(t);
		}
		
		return transcript_ids;
	}
	
	public Location getLocationByGeneId (String geneid) 
	{
		int start = 0;
		int end = 0;
		
		FeatureList fl = gtffeatures.selectByAttribute("gene_id", geneid);
		
		Iterator<FeatureI> fl_i = fl.iterator();
		
		if (fl_i.hasNext())
		{
			Location temp_l = fl_i.next().location();
			if (temp_l.isNegative())
			{
				temp_l = temp_l.opposite();
			}
			start = temp_l.start() ;
			end = temp_l.end();
		}
		while (fl_i.hasNext())
		{
			Location temp_l = fl_i.next().location();

			if (temp_l.isNegative())
			{
				temp_l = temp_l.opposite();
			}

			if(start > temp_l.start())
			{
				 start = temp_l.start();
			}
			
			if(end < temp_l.end())
			{
				end = temp_l.end();
			}
		}
		
		return new Location(start, end);
	}
	
	public char CheckGeneIDLocationStrand(String geneid)
	{
		FeatureList fl = gtffeatures.selectByAttribute("gene_id", geneid);
		
		boolean positive = false;
		boolean negative = false;
		
		Iterator<FeatureI> fl_i = fl.iterator();
		
//		if (fl_i.hasNext())
//		{
//			Location temp_l = fl_i.next().location();
//			if (temp_l.isNegative())
//			{
//				negative = true;
//			}
//			else
//			{
//				positive = true;
//			}
//
//		}
		while (fl_i.hasNext())
		{
			Location temp_l = fl_i.next().location();

			if (temp_l.isNegative())
			{
				negative = true;
			}
			else
			{
				positive = true;
			}
		}
		
		if( positive && !negative)
		{
			return '+';
		}
		else if( !positive && negative)
		{
			return '-';
		}
		else if( !positive && !negative)
		{
			return '.';
		}
		else if( positive && negative)
		{
			return '.';
		}
		return '.';
	}
	
	public String getRefByTranscriptID(String transcriptID)
	{
		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", transcriptID);
		
		Iterator<FeatureI> fl_i = fl.iterator();
		
		if (fl_i.hasNext())
		{
			return fl_i.next().seqname();
		}
		else 
		{
			return null;
		}
		
	}
	
	public Vector<Location> getLocationListByTranscriptId (String geneid) 
	{
		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", geneid);

		Vector<Location> LocationList = new Vector<Location>();
		
		fl=fl.sortByStart();

		Iterator<FeatureI> fl_i = fl.iterator();
		
		while (fl_i.hasNext())
		{
			Location temp_l = fl_i.next().location();
			LocationList.add(temp_l);
		}
		
		if(LocationList.get(0).isNegative())
		{
			Vector<Location> LocationList_reversed = new Vector<Location>();
//			for (int i=LocationList.size()-1; i>=0; i--)
//			{
//				LocationList_reversed.add(LocationList.get(i));
//			}
			for(Location l : LocationList)
			{
				LocationList_reversed.add(l.opposite());
			}
			
			return LocationList_reversed;
		}
		else
		{
			return LocationList;
		}
		
	}

	public Vector<Location> getCDSLocationListByTranscriptId (String geneid) 
	{
		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", geneid);

		Vector<Location> LocationList = new Vector<Location>();
		
		fl=fl.sortByStart();

		Iterator<FeatureI> fl_i = fl.iterator();
		
		while (fl_i.hasNext())
		{
			FeatureI f = fl_i.next();
			if(f.type().compareTo("CDS")==0)
			{
				Location temp_l = f.location();
				LocationList.add(temp_l);
			}
		}
		if(LocationList.size()!=0)
		{
			if(LocationList.get(0).isNegative())
			{
				Vector<Location> LocationList_reversed = new Vector<Location>();
//				for (int i=LocationList.size()-1; i>=0; i--)
//				{
//					LocationList_reversed.add(LocationList.get(i));
//				}
				for(Location l : LocationList)
				{
					LocationList_reversed.add(l.opposite());
				}
				
				LocationList = new Vector<Location>();
				for(Location l : LocationList_reversed)
				{
					LocationList.add(l);
				}
			}
//			else
//			{
//				return LocationList;
//			}
		}
		return LocationList;
	}

	public Vector<Location> getStartCodonLocationListByTranscriptId (String geneid) 
	{
		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", geneid);

		Vector<Location> LocationList = new Vector<Location>();
		
		fl=fl.sortByStart();

		Iterator<FeatureI> fl_i = fl.iterator();
		
		while (fl_i.hasNext())
		{
			FeatureI f = fl_i.next();
			if(f.type().compareTo("start_codon")==0)
			{
				Location temp_l = f.location();
				LocationList.add(temp_l);
			}
		}
		if(LocationList.size()!=0)
		{
			if(LocationList.get(0).isNegative())
			{
				Vector<Location> LocationList_reversed = new Vector<Location>();
//				for (int i=LocationList.size()-1; i>=0; i--)
//				{
//					LocationList_reversed.add(LocationList.get(i));
//				}
				for(Location l : LocationList)
				{
					LocationList_reversed.add(l.opposite());
				}
				
				LocationList = new Vector<Location>();
				for(Location l : LocationList_reversed)
				{
					LocationList.add(l);
				}
			}
//			else
//			{
//				return LocationList;
//			}
		}
		return LocationList;
	}

	public Vector<Location> getStopCodonLocationListByTranscriptId (String geneid) 
	{
		FeatureList fl = gtffeatures.selectByAttribute("transcript_id", geneid);

		Vector<Location> LocationList = new Vector<Location>();
		
		fl=fl.sortByStart();

		Iterator<FeatureI> fl_i = fl.iterator();
		
		while (fl_i.hasNext())
		{
			FeatureI f = fl_i.next();
			if(f.type().compareTo("stop_codon")==0)
			{
				Location temp_l = f.location();
				LocationList.add(temp_l);
			}
		}
		if(LocationList.size()!=0)
		{
			if(LocationList.get(0).isNegative())
			{
				Vector<Location> LocationList_reversed = new Vector<Location>();
//				for (int i=LocationList.size()-1; i>=0; i--)
//				{
//					LocationList_reversed.add(LocationList.get(i));
//				}
				for(Location l : LocationList)
				{
					LocationList_reversed.add(l.opposite());
				}
				
				LocationList = new Vector<Location>();
				for(Location l : LocationList_reversed)
				{
					LocationList.add(l);
				}
			}
//			else
//			{
//				return LocationList;
//			}
		}
		return LocationList;
	}


//	//this method returns only a single reference because if a gene_id maps on multiple reference something is strange
//	public String getReferenceNameByGeneId (String geneid) 
//	{
//		String geneReference = "";
//		Collection<String> col = new Vector<String>();
//			
//		col = gtffeatures.selectByAttribute("gene_id", geneid).attributeValues("seqname");
//		
//		if (col.size()>1)
//		{
//			System.out.println("gene_id " + geneid + " corresponds to multiple gene_name");	
//		}
//		else if(col.size()==1)
//		{
//			geneReference = (new Vector<String>(col)).elementAt(0);			
//		}
//			
//		return geneReference;
//	}


}
