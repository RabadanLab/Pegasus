import java.util.*;

import javagene.seq.*;
import javagene.io.*;
import javagene.util.*;

/**
 * Illustrate use of the AminoAcid and Nucleotide classes.
 */
public class Bio
{
    
    public static void main( String[] args )
    throws Exception
    {
         
        // create a sample sequence
        SeqI seq= new SeqFast( "atatatatacgcgcg", "sample id", "sample desc" );
        
        //
        // The Nucleotide class
        //
         
        // get complements
        String s4= Nucleotide.dnaComplement( seq.toString() );
        String s5= Nucleotide.dnaReverseComplement( seq.toString() );
        
        // Using "ambiguous" nucleotides: 
        // Here we count the purines (IUPAC "R" matches both "A" and "G" )
        int count= 0;
        for( Location loc: seq.bounds() )
        {
            String s= seq.toString( loc );
            if( Nucleotide.R.matches( s ) ) count++;
        }
        
        
        // JavaGene is fluent in RNA
        String s8= Nucleotide.rnaReverseComplement( seq.toString() );
        
        assert Nucleotide.isWobbleComplement( 'G','U' );
        assert Nucleotide.isWobbleComplement( 'G','C' );
        assert ! Nucleotide.isWobbleComplement( 'G','A' );
        
        //
        // The AminoAcid class
        //
         
        // translate nucleotide string to amino acid string
        String aa= AminoAcid.translate( seq.toString() );
        
        // check if codons are synonymous (two ways)
        assert AminoAcid.lookup( "att").isSynonym( "ATC" );      
        assert AminoAcid.CYS == AminoAcid.lookup( "TGC" );
        
        //note that the STOP codon is treated like an AA        
        assert AminoAcid.STOP.isSynonym( "TAA" );
               

    }
    
}
