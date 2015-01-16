package javagene.seq;

import java.util.*;

/**
 * A sequence of either nucleotides or amino acids.
 * Sequences are in IUPAC single letter representation, 
 * ie A,C,G,T,U etc. Sequences are not mutable (once constructed,
 * cannot be changed).
 * <br><br>
 * The sequence objects derived from the SeqI interface have modest functionality; 
 * generally speaking, the good stuff is found in the Location, LocIterator, AminoAcid, and Nucleotide classes.
 */
public interface SeqI
{
	/**
	 * Get the (single-word) identifier of the sequence. This is the first identifier
	 * in a Fasta file.
	 *
	 * @return The single-word identifier.
	 */
	public String id();
	
	/**
	 * Get the text description of the sequence. This is all the text on the description line
	 * of a Fasta file after the id.
	 *
	 * @return The text description.
	 */
	public String description();
	
	/**
	 * Get the location that specifies the bounds of this sequence.
	 *
	 * @return The bounding location.
	 */
	public Location bounds();
	
	/**
	 * Create a new sequence containing only the portion specified by the location parameter.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The subsequence.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */
	public SeqI subseq( Location location );
	
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
    public SeqI subseq( Location location, String id, String description );
    
    /**
	 * Get the IUPAC character string representation of the subsequence specified by 
	 * the location parameter.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The string.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */
	public String toString( Location location );
	
	/**
	 * Get the IUPAC character string representation of this entire sequence.
	 *
	 * @return The text string.
	 */
	public String toString();
}