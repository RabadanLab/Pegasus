package javagene.io;

import javagene.seq.*;
import javagene.util.*;

import java.io.*;
import java.util.*;

/**
* Read and write sequences as Fasta files. A Fasta file has a single header line
* of the following format:
*<pre>
    >id description
</pre>
* where the > symbol is required, the id is a single word, and the description can
* be many words, extending to the end of line. The proper interpretation of the id and
* description fields is not fixed. Data lines immediately follow the header line. Multiple
* sequences in a single file are separated from one another by a blank line.
*
* @author Hanno Hinsch
*/
public class Fasta
{
	
	private Fasta() {};
	
	
	/**
	 * Read a Fasta file containing a single sequence into a SeqBig object.
	 *
	 * @param filename The path to the file.
	 * @return The SeqBig object.
	 * @throws IOException Something went wrong -- check exception detail message.
	 */
	public static SeqBig readBig( String filename )
	throws IOException
	{
		return new SeqBig( filename );
	}
	
	/**
	 * Read a file containing (possibly) multiple sequences.
	 *
	 * @param filename The path to the file.
	 * @param bufferSize An estimate of the size of the largest sequence (hint: unless
	 * you're short on RAM, generously round up. ) 
	 * @return An ArrayList of SeqMem objects.
	 * @throws IOException Something went wrong -- check exception detail message.
	 */
	public static ArrayList<SeqI> readMany( String filename, int bufferSize )
	throws IOException
	{
		Log.log( "Fasta: Reading " + filename );
		
		BufferedReader br = new BufferedReader(new FileReader( filename ));
		
		ArrayList<SeqI> sequences= new ArrayList<SeqI>();
		
		StringBuffer temp= new StringBuffer( bufferSize );
		
		SeqFast seq= null;
		String id= null;
		String description= null;
		
		for( String s= br.readLine(); null != s; s=br.readLine() )
		{
    		s= s.trim();
    		if( 0 < s.length() )	//skip blank lines
    		{
	    		if( s.charAt(0) == '>' )
	    		{
	    			//save current sequence, if any
	    			if( id != null )
	    			{
						seq= new SeqFast(temp.toString(), id, description );	
						
						sequences.add( seq );
						Log.log( "Fasta: " + seq.id() + " loaded." );
						id= null;
						temp.setLength( 0 );	//empty the buffer
					}
	    			
	    			//id, desc line
	    			int i= s.indexOf( ' ' );
	    			if( i < 0 )
	    			{
	    				id= s.substring( 1 );	//whole line after >
	    				description="";
	    			}
	    			else
	    			{
		    			id= s.substring( 1, i );
		    			description= s.substring( i );
		    		}
		    	}
	    		else
	    		{
					//sequence data
					if( id != null )
					{
						temp.append( s );
					}
			    	else
			    	{
			    		throw new IOException( "Unexpected char before >" );
			    	}
				}
		    	
		    }
		
		}
		
		if( id != null )	//grab final sequence
		{
			seq= new SeqFast(temp.toString(), id, description );	
			
			sequences.add( seq );
			Log.log( "Fasta: " + seq.id() + " loaded." );
		}
		
		br.close();
				
		return sequences;
	}
	
	/**
	 * Read a file containing a single sequence.
	 *
	 * @param filename The path to the file.
	 * @param bufferSize An estimate of the size of the sequence (hint: unless
	 * you're short on RAM, generously round up. ) 
	 * @return A SeqMem object.
	 * @throws IOException Something went wrong -- check exception detail message.
	 */
	public static SeqI read( String filename, int bufferSize )
	throws IOException
	{
		return ( readMany( filename, bufferSize )).get(0);
	}
	
	
	/**
	 * Write an ArrayList of sequences to file.
	 *
	 * @param sequences The ArrayList of SeqI sequences.
	 * @param filename The path to the file.
	 * @throws IOException Something went wrong -- check exception detail message.
	 */
	public static void write( ArrayList<SeqI> sequences, String filename )
	throws IOException
	{
		Log.log( "Fasta: Writing " + filename );
		
		BufferedWriter bw= new BufferedWriter( new FileWriter( filename ));
				
		for( SeqI seq: sequences )
		{
			
			bw.write( ">" + seq.id() + " " + seq.description() + " " + (seq.bounds().bioStart()) + "-" + seq.bounds().bioEnd() );
			bw.newLine();
			
			LocIterator li= seq.bounds().iterator( 60, 60 );
			while( li.hasNext() )
			{
				bw.write( seq.toString( li.next()) );
				bw.newLine();
			}
			
            bw.write( seq.toString( li.remainder()) );
            bw.newLine();
                    
			bw.newLine();
		}
		
		bw.close();		
	}
	
	/**
	 * Write a sequence to file.
	 *
	 * @param sequence The SeqI object to write.
	 * @param filename The path to the file.
	 * @throws IOException Something went wrong -- check exception detail message.
	 */
	public static void write( SeqI sequence, String filename )
	throws IOException
	{
		ArrayList<SeqI> list= new ArrayList<SeqI>();		
		list.add( sequence );
		
		write( list, filename );
	}
	
	
	/**
	 * @deprecated
	 */
	public static void main( String args[] )
	throws IOException
	{
		ArrayList<SeqI> seqs= Fasta.readMany( "../rfam/rfam.fasta", 500 );
		
		int max= (seqs.size() < 10)?seqs.size():10;
		for( int i=0; i < max; i++ )
		{
			Log.log( seqs.get(i).id() + seqs.get(i).description() + "\n" );
			Log.log( seqs.get(i).toString() );
		}
		
	}
}