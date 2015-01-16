import java.util.*;

import javagene.seq.*;
import javagene.io.*;
import javagene.util.*;

/**
 * 
 * JavaGene makes it easy to write programs that work with
 * both strands of a nucleotide sequence.
 *
 */
public class Strands
{
    
    public static void main( String[] args )
    throws Exception
    {
 
        //
        // Location objects are strand-aware
        //
        
        // location on positive strand
        Location locPos= new Location( 12, 93 );
        
        // same place on negative strand
        Location locNeg= new Location( -93, -12 );
        
        // check...
        assert ! locPos.isSameStrand( locNeg );
        assert locNeg.isNegative();
        
        // opposite() maps from one to the other...
        assert locPos.equals( locNeg.opposite() );
        
        // plus() maps any location to its positive strand image
        assert locPos.equals( locNeg.plus() );
        
        
        
        //
        // SeqStranded is a sequence object that understands strands
        // and genomic locations. Ordinary sequence objects have a bounding
        // location that starts at 0 and runs along the positive strand.
        // SeqStranded objects can start anywhere (any position on positive 
        // or negative strand).
        //
        // You can turn any sequence into a SeqStranded object if you know the
        // bounding coordinates.
        //
        
        // get sample sequence from file
        SeqI sequence= Fasta.readBig( "./sample/sample.fa" );
        
        // genomic location of sequence in sample fasta file
        Location locGenomicBounds= new Location( 36457708, 39608421 ); 
        
        // wrap it
        SeqStranded seqWrapped= new SeqStranded( sequence, locGenomicBounds, Nucleotide.Type.DNA );
        
        // now genomic locations work properly
        Location locExonInGenome= new Location( 37000000, 37000100 );    //hypothetical exon       
        assert seqWrapped.bounds().contains( locExonInGenome );
        
        // we can get the data on the positive strand
        String sPlus= seqWrapped.toString(  locExonInGenome );
        
        // and OPPOSITE STRAND LOCATIONS AUTOMATICALLY GET THE
        // REVERSE-COMPLEMENT when using toString( location ). 
        // Very handy!
        String sReverseComp= seqWrapped.toString( locExonInGenome.opposite() );
        
        // check my honesty...
        assert sReverseComp.equals( Nucleotide.dnaReverseComplement( sPlus ));
        
        
        
        
        //
        // All the Location methods work on both positive and negative
        // strand locations. 
        //
        // On both strands, end >= start (always).
        //
        // So,
        // upstream() is always upstream 
        // before() is always before
        // ... etc, regardless of strand
        //
        // This greatly simplifies code that needs to work on both strands.
        
        
        
        
        
        
        //
        // Comparing locations on opposite strands:
        //
        // Q: What should the overlaps() method do when one location is on
        // the positive strand, the other on the negative?
        // A: Because one could argue about it, overlaps() just throws an exception.
        //
        // Most methods that consider two locations, ie, contains(), overlaps(), 
        // startsBefore(), union(), intersection(), and so on, throw an exception
        // if the two locations are on opposite strands because the proper behavior is
        // debatable.
        //
        // When the default behavior isn't what I need, I just write a little wrapper method
        // that does exactly what I want. For example, to create a strand-indifferent 
        // version of the overlaps() method, I would map both locations
        // to the positive strand image using the plus() method before calling the
        // overlaps() method.
        //
        
        
        }
}        
