package javagene.seq;

import java.util.*;

/**
* An amino acid object with codon/residue utility methods. 
* Implements the classic 20, plus STOP.
*
* @author Hanno Hinsch 
*/
public class AminoAcid
{
	private int mIndex;			//internal index
	private String mS;			//the print representation
	private String mLookups;	//all lookup representations for this symbol (ie ":GCG:GCC:" )
	
	private static HashMap<String,AminoAcid> mHash;
	
	private static AminoAcid[] mIndexTable;		//NB must agree with indices in intializers (ie, be in same order)
	
	//data from Mount book
	//NB indices are those of AminoAcid.mIndex
	private static int[][] mBlosum62= {			
	/* C */ { 9},
	/* S */ {-1,   4},
	/* T */ {-1,   1,  5},
	/* P */ {-1,  -1, -1,  7},
	/* A */ { 0,   1,  0, -1,  4},
	/* G */ {-3,   0, -2, -2,  0,  6},
	/* N */ {-3,   1,  0, -2, -2,  0,   6},
	/* D */ {-3,   0, -1, -1, -2, -1,   1,  6},
	/* E */ {-4,   0, -1, -1, -1,  2,   0,  2,  5},
	/* Q */ {-3,   0, -1, -1, -1,  2,   0,  0,  2,  5},
	/* H */ {-3,  -1, -2, -2, -2, -2,   1, -1,  0,  0,   8},
	/* R */ {-3,  -1, -1, -2, -1, -2,   0, -2,  0,  1,   0,  5},
	/* K */ {-3,   0, -1, -1, -1, -2,   0, -1,  1,  1,  -1,  2,  5},
	/* M */ {-1,  -1, -1, -2, -1, -3,  -2, -3, -2,  0,  -2, -1, -1,   5},
	/* I */ {-1,  -2, -1, -3, -1, -4,  -3, -3, -3, -3,  -3, -3, -3,   1,  4},
	/* L */ {-1,  -2, -1, -3, -1, -4,  -3, -4, -3, -2,  -3, -2, -2,   2,  2,  4},
	/* V */ {-1,  -2,  0, -2,  0, -3,  -3, -3, -2, -2,  -3, -3, -2,   1,  3,  1,  4},
	/* F */ {-2,  -2, -2, -4, -2, -3,  -3, -3, -3, -3,  -1, -3, -3,   0,  0,  0, -1, 6},
	/* Y */ {-2,  -2, -2, -3, -2, -3,  -2, -3, -2, -1,   2, -2, -2,  -1, -1, -1, -1, 3, 7},
	/* W */ {-2,  -3, -2, -4, -3, -2,  -4, -4, -3, -2,  -2, -3, -3,  -1, -3, -2, -3, 1, 2, 11 }};
	
	public static final AminoAcid ALA;
	public static final AminoAcid ARG;
	public static final AminoAcid ASN;
	public static final AminoAcid ASP;
	public static final AminoAcid CYS;
	public static final AminoAcid GLU;
	public static final AminoAcid GLN;
	public static final AminoAcid GLY;
	public static final AminoAcid HIS;
	public static final AminoAcid ILE;
	public static final AminoAcid LEU;
	public static final AminoAcid LYS;
	public static final AminoAcid MET;
	public static final AminoAcid PHE;
	public static final AminoAcid PRO;
	public static final AminoAcid SER;
	public static final AminoAcid THR;
	public static final AminoAcid TRP;
	public static final AminoAcid TYR;
	public static final AminoAcid VAL;
	public static final AminoAcid STOP;
	public static final AminoAcid UNDEFINED;

	static
	{
		mHash= new HashMap<String,AminoAcid>();
		mIndexTable= new AminoAcid[22];
		
		ALA= new AminoAcid( 4, "A", ":GCA:GCG:GCT:GCC:" );
		mIndexTable[ ALA.mIndex ]= ALA;
		
		ARG= new AminoAcid( 11, "R", ":CGT:CGC:CGA:CGG:AGA:AGG" );
		mIndexTable[ ARG.mIndex ]= ARG;
		
		ASN= new AminoAcid( 6, "N", ":AAT:AAC:" );
		mIndexTable[ ASN.mIndex ]= ASN;
		
		ASP= new AminoAcid( 7, "D", ":GAT:GAC:" );
		mIndexTable[ ASP.mIndex ]= ASP;
		
		CYS= new AminoAcid( 0, "C", ":TGT:TGC:" );
		mIndexTable[ CYS.mIndex ]= CYS;
		
		GLU= new AminoAcid( 8, "E", ":GAA:GAG:" );
		mIndexTable[ GLU.mIndex ]= GLU;
		
		GLN= new AminoAcid( 9, "Q", ":CAA:CAG:" );
		mIndexTable[ GLN.mIndex ]= GLN;
		
		GLY= new AminoAcid( 5, "G", ":GGT:GGC:GGA:GGG:" );
		mIndexTable[ GLY.mIndex ]= GLY;
		
		HIS= new AminoAcid( 10, "H", ":CAT:CAC" );
		mIndexTable[ HIS.mIndex ]= HIS;
		
		ILE= new AminoAcid( 14, "I", ":ATT:ATC:ATA:" );
		mIndexTable[ ILE.mIndex ]= ILE;
		
		LEU= new AminoAcid( 15, "L", ":TTA:TTG:CTT:CTC:CTA:CTG:" );
		mIndexTable[ LEU.mIndex ]= LEU;
		
		LYS= new AminoAcid( 12, "K", ":AAA:AAG:" );
		mIndexTable[ LYS.mIndex ]= LYS;
		
		MET= new AminoAcid( 13, "M", ":ATG:" );
		mIndexTable[ MET.mIndex ]= MET;
		
		PHE= new AminoAcid( 17, "F", ":TTT:TTC:" );
		mIndexTable[ PHE.mIndex ]= PHE;
		
		PRO= new AminoAcid( 3, "P", ":CCA:CCT:CCG:CCC:" );
		mIndexTable[ PRO.mIndex ]= PRO;
		
		SER= new AminoAcid( 1, "S", ":AGT:AGC:TCT:TCC:TCA:TCG:" );
		mIndexTable[ SER.mIndex ]= SER;
		
		THR= new AminoAcid( 2, "T", ":ACA:ACT:ACG:ACC" );
		mIndexTable[ THR.mIndex ]= THR;
		
		TRP= new AminoAcid( 19, "W", ":TGG:" );
		mIndexTable[ TRP.mIndex ]= TRP;
		
		TYR= new AminoAcid( 18, "Y", ":TAT:TAC:" );
		mIndexTable[ TYR.mIndex ]= TYR;
		
		VAL= new AminoAcid( 16, "V", ":GTA:GTC:GTG:GTT:" );
		mIndexTable[ VAL.mIndex ]= VAL;
		
		STOP= new AminoAcid( 20, "X", ":TAA:TAG:TGA:" );		//NB "X" for stop codon (non-standard?)
		mIndexTable[ STOP.mIndex ]= STOP;
		
		UNDEFINED= new AminoAcid( 21, "X", "::" );				//NB "X" for unknown (non-standard?)
		mIndexTable[ UNDEFINED.mIndex ]= UNDEFINED;
	
	}
		
	/**
	*
	*/
	private AminoAcid( int index, String s, String lookups )
	{
		mIndex= index;
		mS= s;
		mLookups= lookups;
		register( lookups, this );
	}
		
	
	/**
	* Find the AminoAcid object corresponding to a three-letter Dna/Rna codon.
	*
	* @param str A string of exactly three letters to look up. Letters can be
	* upper or lower case. Both U and T are accepted.
	* @return An AminoAcid object. Unknown codons return AminoAcid.UNDEFINED. 
	* Stop codons (TGA, TAG, TAA ) return AminoAcid.STOP.
	* @throws IllegalArgumentException The str parameter was not of length 3.
	*/
	public static AminoAcid lookup( String str )
	{
		if (str.length() != 3)
		{
			throw new IllegalArgumentException( str );
		}
		
		AminoAcid s= (AminoAcid) mHash.get( toCanonical( str ) );
		
		if( s == null )
		{
			s= AminoAcid.UNDEFINED;
		}
		
		return s;
	}
	
	/**
	 * convert to uppercase, use T not U
	 */
	 private static String toCanonical( String str )
	 {
	 	return str.toUpperCase().replace( 'U', 'T' );
	 }
	 
	 		 	
	/**
	*
	*/
	private void register( String lookups, AminoAcid s )
	{
		String[] array= lookups.split( ":" );
		
		for( int i= 0; i < array.length ; i++ )
		{
			if( array[i] != null )
				mHash.put( array[i], s );		//key,value
		}
		
	}	
				
	/**
	* Get the IUPAC single letter representation of this AminoAcid.
	* 
	* @return The uppercase single letter representation. "X" is returned
	* for both STOP and UNDEFINED.
	*/
	public String toString()
	{
		return mS;
	}
	
	/**
	* Check if a three-letter Dna/Rna codon codes for this AminoAcid.
	*
	* @param str A string of exactly three letters representing a codon. 
	* Letters can be upper or lower case.
	* Both U and T are accepted and considered equivalent.
	* @return True if str codes for this AminoAcid.
	*/
	public boolean isSynonym( String str )
	{
		String s= ":" + toCanonical( str ) + ":";
		return ( 0 <= mLookups.indexOf( s ));
	}
	
	
	/**
	 * Translate from IUPAC nucleotide string to IUPAC protein string. The nucleotides are considered in
	 * groups of three, starting at the beginning of the string, and translated using the
	 * lookup() and toString() methods. If the string
	 * length is not divisible by three, the final one or two letters are ignored (no error indication).
	 *
	 * @param str A Dna/Rna string.
	 * @return A string of the IUPAC amino acid single letter abbreviations for the translated codons.
	 */
	static public String translate( String str )
	{
		StringBuffer p= new StringBuffer( str.length()/3 );
		for( int i=0; i < str.length() - 2; i+=3 )
		{	
			p.append( lookup( str.substring( i, i + 3 ) ));
		}
		return p.toString();
	}
	
	
	/**
	 * Get BLOSUM62 distance between two AminoAcids. Distances between
	 * AminoAcid.STOP, AminoAcid.UNDEFINED and other AminoAcids are not defined.
	 *
	 * @param a An AminoAcid.
	 * @param b Another AminoAcid.
	 * @return An integer distance based on the Blosum62 table.
	 * @throws IllegalArgumentException A parameter specified
	 * either a STOP or UNDEFINED AminoAcid.
	 */
	public static int blosum62( AminoAcid a, AminoAcid b )
	{
		if( a != AminoAcid.STOP && b != AminoAcid.STOP &&
		    a != AminoAcid.UNDEFINED && b != AminoAcid.UNDEFINED )
		{
			//stay within bounds of triangular matrix
			if( a.mIndex <= b.mIndex )
				return mBlosum62[ b.mIndex ][ a.mIndex ];
			else
				return mBlosum62[ a.mIndex ][ b.mIndex ];
			
		}
		else	
		{
			throw new IllegalArgumentException( a.toString() + " " + b.toString() );
		}
	}
	
	/**
	* Some simple tests of this class.
	* @deprecated
	*/
	static public void main( String args[] )
	{
		
		AminoAcid s= AminoAcid.lookup( "GCA" );
		
		assert s.isSynonym( "GCC" );
		assert AminoAcid.ALA.isSynonym( "GCG" );
		assert !AminoAcid.ALA.isSynonym( "AAA" );
				
        System.out.println( "JavaGene.AminoAcid tests succeeded." );
	}
	
}


