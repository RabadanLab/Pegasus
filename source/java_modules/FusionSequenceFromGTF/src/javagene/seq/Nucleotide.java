package javagene.seq;

import java.util.*;

/**
* A nucleotide object, with Dna/Rna utility methods. Understands  
* both the standard A C T G U and the "ambiguous" IUPAC nucleotides, ie, R Y B D H V N.
*
* @author Hanno Hinsch
*/
public class Nucleotide
{
	/**
	 * Sequence type: DNA or RNA.
	 */
	public enum Type { DNA, RNA };
	
	private String mS;			//the print representation
	private String mLookups;	//all lookup representations for this symbol
	private String mMatches;	//all print representations that are "equal" to this symbol
	
	private static HashMap<String,Nucleotide> mHash;
	
	public static final Nucleotide G;
	public static final Nucleotide C;
	public static final Nucleotide T;
	public static final Nucleotide A;
	public static final Nucleotide U;
	/**
	 * A or G.
	 */
	public static final Nucleotide R;
	/**
	 * C or T/U.
	 */
	public static final Nucleotide Y;
	/**
	 * G or C or T/U (Anything but A).
	 */
	public static final Nucleotide B;
	/**
	 * A or G or T/U (Anything but C).
	 */
	public static final Nucleotide D;
	/**
	 * A or C or T/U (Anything but G).
	 */
	public static final Nucleotide H;
	/**
	 * A or C or G (Anything but T/U).
	 */
	public static final Nucleotide V;
	/**
	 * Any nucleotide.
	 */
	public static final Nucleotide N;

	static
	{
		mHash= new HashMap<String,Nucleotide>();
		
		G= new Nucleotide( "G", ":G:", ":G:R:B:D:V:N:" );
		C= new Nucleotide( "C", ":C:", ":C:Y:B:H:V:N:" );
		T= new Nucleotide( "T", ":T:", ":T:U:Y:B:D:H:N:" );
		A= new Nucleotide( "A", ":A:", ":A:R:D:H:V:N:" );
		U= new Nucleotide( "U", ":U:", ":U:T:Y:B:D:H:N:" );
		R= new Nucleotide( "R", ":R:", ":R:A:G:N:" );
		Y= new Nucleotide( "Y", ":Y:", ":Y:T:C:U:N:" );
		B= new Nucleotide( "B", ":B:", ":B:G:C:T:U:N:" );
		D= new Nucleotide( "D", ":D:", ":D:A:G:T:U:N:" );
		H= new Nucleotide( "H", ":H:", ":H:A:C:T:U:N:" );
		V= new Nucleotide( "V", ":V:", ":V:A:C:G:N:" );
		N= new Nucleotide( "N", ":N:", ":N:A:C:G:T:U:R:Y:B:D:H:V:" );
	}
		
	
					
	/**
	* Find the Nucleotide object corresponding to an IUPAC symbol.
	*
	* @param str One-letter IUPAC nucleotide symbol, either upper or lower case.
	* @return The corresponding Nucleotide object.
	* @throws IllegalArgumentException The parameter did not match a known nucleotide.
	*/	
	public static Nucleotide lookup( String str )
	throws IllegalArgumentException
	{
		Nucleotide s= (Nucleotide) mHash.get( str.toUpperCase() );
		
		if( s == null )
		{
			throw new IllegalArgumentException("Nucleotide <" + str + "> not found." );
		}
		
		return s;
	}
	
	/**
	*
	*/
	private void register( String lookups, Nucleotide s )
	{
		String[] array= lookups.split( ":" );
		
		for( int i= 0; i < array.length ; i++ )
		{
			if( array[i] != null )
				mHash.put( array[i], s );		//key,value
		}
		
	}
	
	/**
	 *
	 */	
	private Nucleotide( String s, String lookups, String matches )
	{
		mS= s;
		mLookups= lookups;
		mMatches= matches;
		
		register( lookups, this );
	}
	
	/**
	 * Get IUPAC representation.
	 *
	 * @return One-letter uppercase IUPAC represenation.
	 */
	public String toString()
	{
		return mS;
	}
	
	/**
	* Check if a Nucleotide object matches this object, allowing ambiguous
	* matches.
	* The nucleotides match if their symbol sets overlap.
	* For example Nucleotide.Y and Nucleotide.T match because Y means "C or T",
	* so "T" is in both sets.
	*
	* @param n The other nucleotide object.
	* @return True if the two symbol sets overlap. 
	*
	*/	
	public boolean matches( Nucleotide n )
	{
		String str= ":" + n.toString().toUpperCase() + ":";
		return ( this == n  || 0 <= mMatches.indexOf( str ));
	}
	
	/**
	* Check if a IUPAC nucleotide symbol matches this object, allowing ambiguous matches.
	* The nucleotides match if their symbol sets overlap.
	* For example Nucleotide.Y and symbol "T" match because Y means "C or T",
	* so "T" is in both sets.
	*
	* @param s The one-letter symbol, upper or lower case.
	* @return True if the two symbol sets overlap. 
	*
	*/	
	public boolean matches( String s )
	{
		String str= ":" + s.toUpperCase() + ":";
		return ( 0 <= mMatches.indexOf( str ));
	}

	/**
	 * Complement a string of IUPAC nucleotides. 
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to complement.
	 * @param type Nucleotide.Type.DNA or Nucleotide.Type.RNA.
	 * @return A complemented string. Case is preserved.
	 */
	public static String complement( String s, Type type  )
	{
			switch( type )
			{
				case RNA:
					return rnaComplement( s );
					
				case DNA:
				default:
					return dnaComplement( s );
					
			}
		
	}
	
	/**
	 * Complement a string of IUPAC DNA nucleotides (output A C T G only).
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to complement.
	 * @return A complemented string, ie A->T, T->A, C->G, G->C, U->A. Case is preserved.
	 */
	public static String dnaComplement( String s )
	{
		char cgene[]= s.toCharArray();	
		
		//complement
		int i;
		for( i=0; i < cgene.length; i ++ )
		{
			switch( cgene[ i ])
			{
				case 'A': cgene[ i ]= 'T'; break;
				case 'T': cgene[ i ]= 'A'; break;
				case 'U': cgene[ i ]= 'A'; break;
				case 'C': cgene[ i ]= 'G'; break;
				case 'G': cgene[ i ]= 'C'; break;
				case 'a': cgene[ i ]= 't'; break;
				case 't': cgene[ i ]= 'a'; break;
				case 'u': cgene[ i ]= 'a'; break;
				case 'c': cgene[ i ]= 'g'; break;
				case 'g': cgene[ i ]= 'c'; break;
			}
				
		}
		
		return new String( cgene );	
	}
	
		
	/**
	 * Complement a string of IUPAC RNA nucleotides (output A C U G only).
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to complement.
	 * @return A complemented string, ie A->U, U->A, C->G, G->C, T->A. Case is preserved.
	 */
	public static String rnaComplement( String s )
	{
		char cgene[]= s.toCharArray();	
		
		//complement
		int i;
		for( i=0; i < cgene.length; i ++ )
		{
			switch( cgene[ i ])
			{
				case 'A': cgene[ i ]= 'U'; break;
				case 'U': cgene[ i ]= 'A'; break;
				case 'T': cgene[ i ]= 'A'; break;
				case 'C': cgene[ i ]= 'G'; break;
				case 'G': cgene[ i ]= 'C'; break;
				case 'a': cgene[ i ]= 'u'; break;
				case 'u': cgene[ i ]= 'a'; break;
				case 't': cgene[ i ]= 'a'; break;
				case 'c': cgene[ i ]= 'g'; break;
				case 'g': cgene[ i ]= 'c'; break;
			}
				
		}
		
		return new String( cgene );	
	}

	/**
	 * Reverse and complement a string of IUPAC DNA nucleotides.
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to reverse complement.
	 * @param type Nucleotide.Type.DNA or Nucleotide.Type.RNA.
	 * @return A complemented string. Case is preserved.
	 */
	public static String reverseComplement( String s, Type type  )
	{
			switch( type )
			{
				case RNA:
					return rnaReverseComplement( s );
					
				case DNA:
				default:
					return dnaReverseComplement( s );
					
			}
		
	}
	
	
	/**
	 * Reverse and complement a string of IUPAC DNA nucleotides (A C T G only).
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to reverse and complement.
	 * @return The reverse-complemented string, ie A->T, T->A, C->G, G->C. Case is preserved.
	 */
	public static String dnaReverseComplement( String s )
	{
		StringBuffer r= new StringBuffer( dnaComplement( s ) );
		return r.reverse().toString();
	}
	
	/**
	 * Reverse and complement a string of IUPAC RNA nucleotides (A C U G only).
	 *
	 * @param s The string of one-letter upper or lower case IUPAC nucleotides to reverse and complement.
	 * @return The reverse-complemented string, ie A->U, U->A, C->G, G->C. Case is preserved.
	 */
	public static String rnaReverseComplement( String s )
	{
		StringBuffer r= new StringBuffer( rnaComplement( s ) );
		return r.reverse().toString();
	}

	/**
	 * Check if two symbols complement one another, accepting wobble pairs (G-U, G-T) in addition to 
	 * standard (A-T, A-U, G-C) pairs.
	 *
	 * @param a A nucleotide symbol, either upper or lower case.
	 * @param b A nucleotide symbol, either upper or lower case.
	 * @return True if the nucleotides complement one another.
	 */
	public static boolean isWobbleComplement( char a, char b )
	{
		a= Character.toUpperCase( a );
		b= Character.toUpperCase( b );
		
		return ( a == 'A' && b == 'T' ) ||
		       ( a == 'A' && b == 'U' ) ||
		       ( a == 'T' && b == 'A' ) ||
		       ( a == 'T' && b == 'G' ) ||
		       ( a == 'U' && b == 'A' ) ||
		       ( a == 'U' && b == 'G' ) ||
		       ( a == 'G' && b == 'C' ) ||
		       ( a == 'G' && b == 'U' ) ||
		       ( a == 'G' && b == 'T' ) ||
		       ( a == 'C' && b == 'G' );
	}
	
	
	/**
	 * Check if two symbols complement one another. Accepts only
	 * standard (A-T, A-U, G-C) pairs.
	 *
	 * @param a A nucleotide symbol, either upper or lower case.
	 * @param b A nucleotide symbol, either upper or lower case.
	 * @return True if the nucleotides complement one another.
	 */
	public static boolean isComplement( char a, char b )
	{
		a= Character.toUpperCase( a );
		b= Character.toUpperCase( b );
		
		return ( a == 'A' && b == 'T' ) ||
		       ( a == 'T' && b == 'A' ) ||
		       ( a == 'U' && b == 'A' ) ||
		       ( a == 'A' && b == 'U' ) ||
		       ( a == 'G' && b == 'C' ) ||
		       ( a == 'C' && b == 'G' );
	}


	/**
	 * @deprecated
	 */
	static public void main( String args[] )
	throws Exception
	{
		
		Nucleotide A= Nucleotide.lookup( "A" );
		
		assert Nucleotide.lookup("A").matches( Nucleotide.lookup("a"));
		assert A == Nucleotide.A;
		assert A.matches( "a" );
		assert A.matches( Nucleotide.A );
		assert A.matches( Nucleotide.D );
		assert !A.matches( Nucleotide.B );
		assert A.matches( Nucleotide.N );
		assert A.matches( Nucleotide.R );
		
		System.out.println( "JavaGene.Nucleotide tests succeeded." );
	}
	
}


