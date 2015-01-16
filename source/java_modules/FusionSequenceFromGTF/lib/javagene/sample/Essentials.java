import java.util.*;

import javagene.seq.*;
import javagene.io.*;
import javagene.util.*;

/**
 *
 * Quick introduction to the Sequence and Location objects.
 *
 */
public class Essentials
{
    
    public static void main( String[] args )
    throws Exception
    {
 
        //
        // A minimal DNA sequence -- our "fakegene"
        SeqI sequence= new SeqFast( 
                                    "GGGATGAAAGGGCCCTAGCCC",        //the sequence data
                                    "fakegene",                     //the one-word id
                                    "A short fake gene" );          //the text description
         //
         // About sequence objects:
         //
         // 1. They implement the SeqI interface.
         // 2. They don't do much.
         // 3. They live in Fasta files.
         // 4. They live either in memory (SeqFast) or live on disk (SeqBig).
         // 5. If it's nucleotide and you need to work with both strands,
         //    or you're using genomic coordinates, you'll need SeqStranded.
         //
         
         
         //
         // A Location object refers to a location on a sequence.
         // In this case, the Start (ATG) codon on our fake gene.
         Location locStart= new Location( 3, 6 );
         
         // About Location objects:
         //
         // 1. They represent a range (a start and an end ).
         // 2. You need them all the time.
         // 3. Nearly every method in JavaGene uses a Location for 
         //    some purpose or another.
         // 4. The first index is 0, not 1 (like Java strings)
         // 5. The Location class is designed to let you do sequence
         //    arithmetic without doing any arithmetic (you'll see).
         
  
         
         // get the sequence data at a location (our start codon)
         assert "ATG".equals( sequence.toString( locStart ));
          
                 
         // Get everything before the Start codon.
         // Explanation: sequence.bounds() gets the location of the entire sequence,
         // prefix( locStart ) gets everything in bounds() before locStart.
         Location locUTR= sequence.bounds().prefix( locStart );       
         assert "GGG".equals( sequence.toString( locUTR ));
         
                  
         // fancier: find the STOP codon in our sequence using
         // an iterator (created by the window() method).
         Location locStop= null;
         for( Location loc: sequence.bounds().window( 3, 3 ) )
         {
            if( AminoAcid.STOP.isSynonym( sequence.toString( loc )))
            {
                locStop= loc;   //found it
                break;
            }          
         }
         
         
         // Use prefix() and suffix() methods to pick out portions of a location.
         // Explanation: suffix( - 6 ) gets the last 6 symbols in location,
         // prefix( 3 ) gets the first three,
         // so here we get the first 3 of the last 6...
         assert locStop.equals( sequence.bounds().suffix( - 6 ).prefix( 3 ));
        
         
         // ...final note: the Location class has many other useful methods 
         // that are not shown in this sample. See the JavaDoc.
    }
    
}
