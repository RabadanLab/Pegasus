package javagene.seq;

import javagene.util.*;

import java.util.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.regex.*;


/**
 * A SeqI implementation capable of handling large sequences that do
 * not fit into memory. The sequence data is kept on disk and accessed
 * as a memory-mapped file.
 * The sequence file must be in Fasta format.
 * Assumes all data lines (except last) are of equal length.
 * (This is checked when the file is opened.)
 *
 */
 
public class SeqBig implements SeqI
{
	String mId;
	String mDesc;
		
	ByteBuffer mBuffer;
	
	int mSeqStart= -1;
	int mSeqEnd= -1;
	
	int mFirstDataLine;
	int mLineLength;
	int mLinePad;
	int mFullLineCount;
	int mLastLineLength;
		

	private SeqBig() {};
	
	/**
	 * Construct SeqBig object given Fasta filename. This is
	 * usually done through the Fasta class.
	 *
	 * @param filename The path and name of the Fasta file.
	 * @throws IOException - Check the exception text for details.
	 */
	public SeqBig( String filename )
	throws IOException
	{
		Log.log( "SeqBig: Reading " + filename );
		
		RandomAccessFile raf= new RandomAccessFile( filename, "r" );
		long size= raf.length();
		FileChannel fc= raf.getChannel();
		
		mBuffer= fc.map( FileChannel.MapMode.READ_ONLY, 0, size );	
		
		fc.close();
		raf.close();
		
		//parse first line
		String firstLine= toString( 0, nextEOL( 0 ));
		Pattern pattern= Pattern.compile( "^>(\\S+)\\s*(.*)$" );
		Matcher matcher= pattern.matcher( firstLine );
		if( !matcher.find())
		{
			Log.log( firstLine );
			throw new IOException( "Failed to parse description line in Fasta file." );
		}
		
		mId= matcher.group( 1 );
		mDesc= matcher.group( 2 );
			
		//deduce padlength
		int n= nextEOL( 0 );
		mLinePad= findPadLength( n );
		mFirstDataLine= n + mLinePad;
		
		mLineLength= nextEOL( mFirstDataLine ) - mFirstDataLine;
		
		//count number of lines and verify format
		mFullLineCount= 0;
		int pos= mFirstDataLine;
		while( true )
		{
			int eol= nextEOL( pos );
							
			if( eol - pos < mLineLength || pos > mBuffer.limit())
			{
				//check that there really is no more sequence data
				for( int i= eol; i < mBuffer.limit(); i++ )
				{
					char c= Character.toChars( (int) mBuffer.get( i ) )[0];
					
					if( Character.isLetterOrDigit( c ))
					{
						Log.log( "pos=" + pos + "; eol=" + eol );
						Log.log( toString( pos, eol ));
						throw new IOException( "BigSequence.java cannot handle Fasta files with uneven line lengths." );
					}
				}
				
				mLastLineLength= eol - pos;			
				break;
			}
			else if( eol - pos > mLineLength )
			{
				Log.log( "pos=" + pos + "; eol=" + eol );
				Log.log( toString( pos, eol ));
				throw new IOException( "BigSequence.java cannot handle Fasta files with uneven line lengths." );
			}
			else
			{
				mFullLineCount++;
				pos= eol + mLinePad;
			}
			
		}
		
		mSeqStart= 0;
		mSeqEnd= mLineLength * mFullLineCount + mLastLineLength;
		
	}
		
	/**
	 *
	 */
	private SeqBig( SeqBig seq, Location location, String id, String description )
	{
		mId= id;
		mDesc= description;
		
		mBuffer= seq.mBuffer;
		
		mFirstDataLine= seq.mFirstDataLine;
		mLineLength= seq.mLineLength;
		mLinePad= seq.mLinePad;
		mFullLineCount= seq.mFullLineCount;
		mLastLineLength= seq.mLastLineLength;
		
		mSeqStart= location.start();
		mSeqEnd= location.end();
		
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
     * Create a new sequence containing only the portion specified by the location parameter.
     * Shares the underlying data (does not make a copy).
     *
     * @param location The portion of the sequence to get.
     * @return The subsequence.
     * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
     * of the sequence.
     */
    public SeqBig subseq( Location location )
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqBig handles only positive strand locations." );
        }
        
        return new SeqBig( this, location, mId, mDesc );
	}
		
    /**
     * Create a new sequence containing only the portion specified by the location parameter,
     * with the specified id and description. Shares the underlying data (does not make a copy).
     *
     * @param location The portion of the sequence to get.
     * @param id The new single-word id.
     * @param description The new description line.
     * @return The subsequence.
     * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
     * of the sequence.
     */
    public SeqBig subseq( Location location, String id, String description )
    {
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqBig handles only positive strand locations." );
        }
        
        return new SeqBig( this, location, id, description );
    }
        
    /**
	 * Get the bounding location of this sequence.
	 * SeqBig objects always start at index 0.
	 * To use a coordinate
	 * system that starts at a non-zero
	 * genomic location, wrap
	 * the SeqBig object in a SeqStranded object.
	 *
	 * @return The bounding location.
	 */
	public Location bounds()
	{
		return new Location( mSeqStart, mSeqEnd );	
	}
	
	/**
	 * Get the IUPAC character string representation of this entire sequence.
	 * Makes a copy of the underlying data.
	 *
	 * @return The text string.
	 */
	public String toString()
	{
		return toString( bufferPos( mSeqStart ), bufferPos( mSeqEnd ) );
	}
	
	/**
	 * Get the IUPAC character string representation of the subsequence specified by 
	 * the location parameter. Makes a copy of the underlying data.
	 *
	 * @param location The portion of the sequence to get.
	 * @return The string.
	 * @throws IndexOutOfBoundsException The location parameter was not within the bounding location
	 * of the sequence.
	 */
	public String toString( Location location)
	{
        if( location.isNegative())
        {
            throw new IndexOutOfBoundsException( "Location is on negative strand; SeqBig handles only positive strand locations." );
        }
        
        return toString( bufferPos( location.start() ), bufferPos( location.end() ) );
	}
	
	/**
	 * Convert genomic index to buffer index
	 */
	private int bufferPos( int n )
	{
		//FIXME check bounds
		return mFirstDataLine + (n % mLineLength) + (n / mLineLength) * (mLineLength + mLinePad);
	}	
	
	/**
	 * Find index of next EOL after pos.
	 */
	private int nextEOL( int pos )
	{
		for( int i= pos; i < mBuffer.limit(); i++ )
		{
			char c= Character.toChars( (int) mBuffer.get( i ) )[0];
			
			if( c == '\r' || c == '\n' )
			{
				return i;
			} 
		}
		
		return mBuffer.limit();
	}
		
	
	/**
	 * given pos at EOL, how many chars til next line?
	 */
	private int findPadLength( int pos )
	{
		for( int i= 0; i < mBuffer.limit(); i++ )
		{
			char c= Character.toChars( (int) mBuffer.get( pos + i ) )[0];
			
			if( c != '\r' && c != '\n' )
			{
				return i;
			} 
		}
		
		return mBuffer.limit() - pos;
	}
		
		
	/**
	 * NB Uses buffer coordinates. Strips cr/lf's.
	 */
	private String toString( int start, int end )
	{
		StringBuffer sb= new StringBuffer(  end - start );
		
		for( int i= start; i < end; i++ )
		{
			char c= Character.toChars( (int) mBuffer.get( i ) )[0];
			
			if( c != '\r' && c != '\n' )
			{
				sb.append( c );
			} 
		}
			
		return sb.toString();
	}
	
		
	/**
	 * @deprecated
	 */
	public static void main( String args[] )
	throws Exception
	{
		
		SeqBig bs= new SeqBig( "..\\test\\test.fa" );
		//BigSequence bs= new BigSequence( "..\\ucsc-hg16-mm3\\mm3chr5.fa" );
		
		Log.log( "LineCount=" + bs.mFullLineCount );
		Log.log( "FirstLine=" + bs.mFirstDataLine);
		Log.log( "LineLength=" + bs.mLineLength);
		Log.log( "LinePad=" + bs.mLinePad);
		Log.log( "LastLineLen=" + bs.mLastLineLength);

		Log.log( "pos 779-784: " + bs.subseq( new Location( 779, 784 )).toString());
	}
	
		
}