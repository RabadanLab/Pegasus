import java.util.*;

import javagene.seq.*;
import javagene.io.*;
import javagene.util.*;

/**
 *
 * Get the translated protein sequence for each of the
 * genes in a GFF file.
 *
 * Illustrates reading of Fasta and GFF files, selection
 * of features that have specific attributes, splicing,
 * translation into protein sequence.
 *
 */
public class GenesToProteins
{
    
    public static void main( String[] args )
    throws Exception
    {
 
        // 
        // First read the sample Fasta file into a sequence object.
        //      
        // The 20000 parameter is a hint to the read() method to 
        // expect a sequence of around this size.
        // 
        // This sequence is a portion of human chr21
        //
        // NOTE - you may need to edit this filepath, depending on where
        // you installed the javagene.jar file and samples!
        //     
        SeqI seqRaw= Fasta.read( "./sample/sample.fa", 20000 );

        //
        // Since the genomic location and type of the sequence are
        // not in the Fasta file, we have to add this information here 
        // and create a SeqStranded object.
        // (yes, sorry, I know this step is confusing...)
        SeqI seqChr21Portion= new SeqStranded( seqRaw, 
                                               Location.fromBio( 36457709, 39608421, '+'  ), 
                                               Nucleotide.Type.DNA );
        
        //
        // Now read in a GTF file full of "features".
        //        
        // This sample file contains features describing known genes in our sequence.
        //
        FeatureList listChr21Genes= Gff.read( "./sample/sample.gtf" );
        
        //
        // Get a list of the genes in the GTF file:
        // attributeValues() returns a collection of all the values found for 
        // the specified attribute key; in this case, the "gene_id" key of a GTF feature.
        Collection<String> geneIds= listChr21Genes.attributeValues( "gene_id" );
        
        //
        // We want to show the proteins that are coded
        // so iterate over this collection of gene ids. 
        for( String gene: geneIds )
        {
            // we want the translated protein residues for each gene...
            
            // get all the individual features (exons, CDS regions, etc.) of this gene
            FeatureList listGene= listChr21Genes.selectByAttribute( "gene_id", gene ); 
            
            // now select only the coding regions of this gene
            listGene= listGene.selectByType( "CDS" );
            
            // sort them
            listGene= listGene.sortByStart();
            
            // splice them together
            String codingRegion= listGene.splice( seqChr21Portion );
            
            // all of the above statements could have been chained, like so:
            // codingRegion= listChr21Genes.selectByAttribute( "gene_id", gene ).selectByType( "CDS" ).sortByStart().splice( seqChr21Portion );
            
            // translate them
            String protein= AminoAcid.translate( codingRegion );
            
            // show translations on screen 
            Log.log(  gene + "= " + protein );
            Log.log();
        }
        
    }
    
}
