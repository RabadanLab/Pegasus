package javagene.seq;
 
import java.io.*;
import java.util.*;

/**
 * A SeqI implementation that keeps the entire sequence in memory.
 * The sequence must be ungapped (ie contiguous) and is 
 * untyped (ie might be dna, rna, codons, ambiguous, etc.).
 *
 *
 * @author Hanno Hinsch
 */
public class SeqFast implements SeqI
{
	String mData="";
	String mId="";
	String mDesc="";	
	
	private SeqFast() {};
	
	/**
	 * Construct SeqFast object from string data.
	 *
	 * @param seqData
	 * @param id
	 * @param description
	 */
	public SeqFast( String seqData, String id, String description )
	{
		mData= seqData;
		mId= id;
		mDesc= description;
	}


	
	/**
	 * Construct SeqFast object from an arbitrary SeqI object.
	 * This can be useful when working with subsequences
	 * of SeqBig objects, as SeqFast supports additional methods
	 * that SeqBig does not.
	 *
	 * @param sequence A SeqI object.
	 */
	public SeqFast( SeqI sequence )
	{
		mData= sequence.toString();
		mId= sequence.id();
		mDesc= sequence.description();
		
		assert sequence.bounds().length() == mData.length(): "S=" + sequence.bounds().length() +"; D=" + mData.length();
	}	
	
	
	public String id()
	{
		return mId;
	}
	
	
	public String description()
	{
		return mDesc;
	}

	/**
	 * Get the bounding location of this sequence.
	 * SeqFast objects always start at index 0.
	 * To use a coordinate
	 * system that starts at a non-zero
	 * genomic location, wrap
	 * the SeqFast object in a SeqStranded object.
	 * 
	 *
	 * @return The bounding location.
	 */
	public Location bounds()
	{
		return new Location( 0, mData.length());
	}
	
	/**
     * Create a new sequence containing only the portion specified by the location parameter.
     * May make a copy of the underlying data.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The new sequence.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */
	public SeqFast subseq( Location location )
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqFast handles only positive strand locations." );
        }
        
        return new SeqFast( mData.substring( location.start(), location.end()), mId, mDesc );
	}
	
    /**
     * Create a new sequence containing only the portion specified by the location parameter,
     * with the specified id and description.
     *
     * @param location The portion of the sequence to get.
     * @param id The new single-word id.
     * @param description The new description line.
     * @return The subsequence.
     * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
     * of the sequence.
     */
    public SeqFast subseq( Location location, String id, String description )
    {
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqFast handles only positive strand locations." );
        }
        
        return new SeqFast( mData.substring( location.start(), location.end()), id, description );
    }
    
    public String toString()
	{
		return mData;
	}
	
	public String toString( Location location ) 
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqFast handles only positive strand locations." );
        }
        
        return mData.substring( location.start(), location.end());
	}
	
	/**
	 * Get the IUPAC representation of the entire sequence, bracketing the portion
	 * specified by the location parameter with the specified begin and end tags. Use this
	 * to create HTML-ready output of sequences with a specific location identified.
	 * For example, try using "&lt;span style='color:red'&gt;" and "&lt;/span&gt;" as tags.
     * 
     * @param location Location of the interesting portion of the sequence.
     * @param beginTag The string to put before the portion of interest.
     * @param endTag The string to put after the portion of interest.
     * @return The formatted string.
     * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
     * of the sequence.
     */
	public String toTaggedString( Location location, String beginTag, String endTag ) 
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqFast handles only positive strand locations." );
        }
        
        return toString( bounds().prefix( location )) + beginTag + toString( location ) + endTag + toString( bounds().suffix( location ));
	}
	
	
	/**
	 * Get the IUPAC character string representation of the subsequence specified by 
	 * the location parameter, in reverse order. Note: just reversed, not complemented.
	 * See the Nucleotide class for complementation methods.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The reversed string.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */	
	public String toReverseString( Location location ) 
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqFast handles only positive strand locations." );
        }
        
        StringBuffer s= new StringBuffer( mData.substring( location.start(), location.end()));
		return s.reverse().toString();
	}
				
	/**
	 * Return the location of the parameter string in this sequence.
	 * Returns null if the searchText is not found.
	 *
	 * @param searchText The string to search for.
	 * @return The first location of the searchText in the string,
	 * or null if the text is not found.
	 */
	 public Location locationOf( String searchText )
	 {
	 	int i= mData.indexOf( searchText );
	 	if( i >= 0 )
	 	{
		 	return new Location( i, searchText.length() );
	 	}
	 	else
	 	{
	 		return null;
	 	}
	 }
	 	
	 		
	/**
	* @deprecated
	*/  	
	public static void main( String args[] )
	{
		//Sequence seq= new Sequence( new Location( 12, 19) , "NATGAUGUUGUUGUUAAUC" );
		
		//assert seq.get( 15, 3 ).equals( "GAU" );
		//assert seq.indexOf( "AUC", 0 ) == seq.location().end() - 3;
		
		//Sequence sub= seq.subseq( new Location( 12, 21 ));
		//assert sub.get( 12, 3 ).equals( "NAT" );
		//assert sub.get( 14 ) == 'T';
		
		//assert (seq.location().end() - 3) == seq.framedIndexOf( "AUC", 13 );
		//assert -1 == seq.framedIndexOf( "AUC", 12 );
		//Log.log( " " + seq.framedLastIndexOf( "GUU", new Location( seq.location().start(), seq.location().end() ) ));
		//assert 24 == seq.framedLastIndexOf( "GUU", new Location( seq.location().start(), seq.location().end() )); 
	}
}