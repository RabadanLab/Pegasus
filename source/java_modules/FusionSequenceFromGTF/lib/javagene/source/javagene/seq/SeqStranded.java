package javagene.seq;

import java.util.*;

/**
 * A SeqI wrapper that presents a stranded coordinate system 
 * mapped to an underlying SeqI object.
 * <br><br>
 * The bounds of a SeqStranded object can be set to an arbitrary location; 
 * they do not necessarily start at 0, and they may be on the reverse
 * strand.
 * <br><br>
 * In JavaGene's coordinate system a numerical increase
 * is always upstream. The direction (strand) of the underlying sequence
 * is indicated by the strand of the bounding location.
 *<br><br>
 * The sequence type (RNA or DNA) is always set in constructor.
 *
 */
public class SeqStranded implements SeqI
{
	
		
	private SeqI mSeq;
	private Location mLocation;
	private Nucleotide.Type mType;
	
	//The genomic position of the very first symbol in base sequence data
	private int mBaseSeqOffset;
		
	private SeqStranded(){};
	
	/**
	 * Construct from an arbitrary SeqI object. Shares the underlying data.
	 * The bounding location starts at 0.
	 *
	 * @param baseSequence A SeqI object.
	 * @param type Nucleotide.Type.DNA or Nucleotide.Type.RNA.
	 */
	public SeqStranded( SeqI baseSequence, Nucleotide.Type type )
	{
		//FIXME what if SeqStranded is passed in? bounds() won't be right!
		mSeq= baseSequence;
		mLocation= baseSequence.bounds();
		mType= type;
		
		mBaseSeqOffset= baseSequence.bounds().start();
	}

	
	/**
	 * Construct a clone. Shares the underlying data.
	 *
	 * @param baseSequence A SeqI object.
	 */
	public SeqStranded( SeqStranded baseSequence )
	{
		mSeq= baseSequence.mSeq;
		mLocation= baseSequence.mLocation;;		
		mBaseSeqOffset= baseSequence.mBaseSeqOffset;
		mType= baseSequence.mType;
	}
		
	
			
	/**
	 * Construct from SeqI using specified location as the bounding
	 * location.
	 * Shares the underlying data.
	 * The start index of the bounding location may be non-zero.
	 * The bounding location may be on either strand.
	 * Note: the length of the sequence must equal
	 * the length of the location.
	 *
	 * @param baseSequence A SeqI object.
	 * @param location The bounding location.
	 * @param type Nucleotide.Type.DNA or Nucleotide.Type.RNA
	 * @throws IllegalArgumentException The sequence length and location length
	 * were not equal.
	 */
	public SeqStranded( SeqI baseSequence, Location location, Nucleotide.Type type )
	{
		if( baseSequence.bounds().length() != location.length() )
		{
			throw new IllegalArgumentException( "Sequence length and location length must match." );
		}
		
		mSeq= baseSequence;
		mLocation= location;	
		mBaseSeqOffset= location.start();
		mType= type;
	}
	
	
					
	public String id() { return mSeq.id(); };
	
	public String description(){ return mSeq.description(); };
	
	/**
	 * Get the bounding location of this sequence.
	 * The location can start at any index, on either
	 * strand.
	 *
	 * @return The bounding location.
	 */
	public Location bounds() { return mLocation; };
				
		
	/**
	 * Get the base sequence that was specified in the constructor.
	 *
	 * @return The base sequence.
	 */	
	public SeqI baseSeq()
	{
		return mSeq;
	}
	
	/**
     * Create a new sequence containing only the portion specified by the location parameter.
     * Shares the underlying data (does not make a copy).
	 *
	 * The specified location must be contained by the bounding location,
	 * and on the same strand. Consider using the toString( Location ) method
	 * if the location might be on either strand.
	 * 
	 * @param location The location of the subsequence.
	 * @return The subsequence.
	 * @throws IndexOutOfBoundsException The bounding location did not
	 * contain the specified location.
	 */
	public SeqStranded subseq( Location location )
	{
		if( ! mLocation.contains( location ))
		{
			throw new IndexOutOfBoundsException( "Specified location not within bounding location." );
		}
		
		SeqStranded g= new SeqStranded( this );
		
		g.mLocation= location;
		
		return g;
	}

    /**
     * Create a new sequence containing only the portion specified by the location parameter,
     * with the specified id and description.
     * Shares the underlying data (does not make a copy).
     *
     * The specified location must be contained by the bounding location,
     * and on the same strand. Consider using the toString( Location ) method
     * if the location might be on either strand.
     * 
     * @param location The location of the subsequence.
     * @param id The new single-word id.
     * @param description The new description line.
     * @return The subsequence.
     * @throws IndexOutOfBoundsException The bounding location did not
     * contain the specified location.
     */
    public SeqStranded subseq( Location location, String id, String description )
    {
        if( ! mLocation.contains( location ))
        {
            throw new IndexOutOfBoundsException( "Specified location not within bounding location." );
        }
        
        SeqStranded g= new SeqStranded( mSeq.subseq( mapToBase( location ), id, description ), location, mType );
        
        return g;
    }
    
	/**
	 * Get the IUPAC character string representation of the subsequence specified by 
	 * the location parameter. Returns the reverse complement if the specified location
	 * and bounding locations are on opposite strands. Uses the sequence type (Rna or Dna) 
	 * to choose whether T or U should be used when complementing (if necessary). 
	 * The specified location must be 
	 * within the bounding location.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The string.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */
	public String toString( Location location )
	{
		if( location.isSameStrand( bounds() ))
		{
			return mSeq.toString( mapToBase( location ));
		}
		else
		{
			return Nucleotide.reverseComplement( mSeq.toString( mapToBase( location )), mType);
		}
	}


	
	/**
	 * Get the IUPAC character string representation of this entire sequence,
	 * as defined by the bounding location.
	 *
	 * @return The text string.
	 */
	public String toString()
	{
		return mSeq.toString( mapToBase( mLocation ));
	}
	
	/**
	 * Map a genomic location to the equivalent location
	 * on the base sequence. If the specified location
	 * is on the opposite strand, this method first converts
	 * the specified location to the other strand, and then
	 * maps the result to the base sequence.
	 *
	 * @param genomicLocation The location to map. Must be within the location bounds,
	 * though both strands are allowed.
	 * @return The equivalent location on the base sequence.
	 * @throws IndexOutOfBoundsException The genomicLocation did not map
	 * within the base sequence bounds.
	 */
	public Location mapToBase( Location genomicLocation )
	{
		if( ! bounds().isSameStrand( genomicLocation ))
		{
			genomicLocation= genomicLocation.opposite();
		}
		
		if( ! mLocation.contains( genomicLocation ))
		{
			throw new IndexOutOfBoundsException( "Specified location not within bounding location." );
		}
				
		return new Location( genomicLocation.start() - mBaseSeqOffset, genomicLocation.end() - mBaseSeqOffset );
	}
	
	/**
	 * Map a location on the base sequence to the equivalent
	 * genomic location. The returned genomic location is 
	 * always on the same strand as the location bounds.
	 * 
	 * @param baseLocation The location on the base string.
	 * @return A location within the location bounds.
	 * @throws IndexOutOfBoundsException The specified baseLocation could not
	 * be mapped within the SeqStranded object's location bounds.
	 */
	public Location mapFromBase( Location baseLocation )
	{
		if( ! mSeq.bounds().contains( baseLocation ) )
		{
			throw new IndexOutOfBoundsException( "Specified location not within bounding location of base sequence." );			
		}
		
		return new Location( baseLocation.start() + mBaseSeqOffset, baseLocation.end() + mBaseSeqOffset );
	}
	


}