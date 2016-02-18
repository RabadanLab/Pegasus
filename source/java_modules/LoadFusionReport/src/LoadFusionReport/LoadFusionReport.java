package LoadFusionReport;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
//import java.util.ArrayList;

import org.hsqldb.Server;

import gnu.getopt.Getopt;

public class LoadFusionReport {
	private Connection conn;
	private Server server;

	public LoadFusionReport()
	{
		
	}

	public void run(
			String fusion_input_report_file, 
			String sample_name,
			String sample_type,
			String fusion_program, 
			String db_filename_prefix,
			String fusion_table) throws ClassNotFoundException, SQLException
	{

		//Class.forName("org.hsqldb.jdbcDriver");
		
		//conn = DriverManager.getConnection("jdbc:hsqldb:hsql://localhost/xdb", "SA", "");
		
        try {
        	
        	if(IsSampleAlreadyPresent(sample_name, fusion_program, fusion_table))
        	{
        		System.err.println("Sample " + sample_name + " already loaded in the database....");
        		System.exit(1);
        	}

            // add some rows - will create duplicates if run more then once
            // the id column is automatically generated
        	String fusion_report_line;
			File fusion_input_report_fr = new File(fusion_input_report_file);
			BufferedReader fusion_input_report_br = new BufferedReader(new FileReader(fusion_input_report_fr));
			int count_fusion = 0;
			while((fusion_report_line = fusion_input_report_br.readLine())!=null)
			{
				if(fusion_report_line.matches("^#.*$"))
				{
					continue;
				}
				
				if(fusion_report_line.split("\t")[0].compareTo("cluster_id")==0)
				{
					continue;
				}
				
				if(fusion_program.compareTo("defuse")==0)
				{
					insertDefuseRecord(fusion_report_line, sample_name, sample_type, fusion_table);
				}
				
				if(fusion_program.compareTo("chimerascan")==0)
				{
					insertChimerascanRecord(fusion_report_line, sample_name, sample_type, fusion_table);
				}
				
				if(fusion_program.compareTo("bellerophontes")==0)
				{
					insertBellerophontesRecord(fusion_report_line, sample_name, sample_type, fusion_table);
				}
				
				if(fusion_program.compareTo("general")==0)
				{
					insertGeneralFusionRecord(fusion_report_line, sample_name, sample_type, fusion_table);
				}
				
				count_fusion++;
			}
			int loaded_fusion = CountLoadedFusions(fusion_program, sample_name, sample_type, fusion_table);
			boolean exit_for_failure = false;
			while (loaded_fusion != count_fusion && !exit_for_failure)
			{
				System.err.print("Error: Fusion not completely loaded...\nRetrying:\n");
				Thread.sleep(5000);
				loaded_fusion = CountLoadedFusions(fusion_program, sample_name, sample_type, fusion_table);
				if(loaded_fusion != count_fusion)
				{
					exit_for_failure = true;
					System.err.print("Error: Fusion not completely loaded...\nExiting:\n");

				}
			}
			
			System.err.print("Loaded " + loaded_fusion + " fusions\n");
			
        	fusion_input_report_br.close();
        	//conn.close();
        	
            // at end of program
            //this.shutdown();
//        } catch (SQLException ex3) {
//            ex3.printStackTrace();
        } catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private int CountLoadedFusions(
		String fusion_program,
		String sample_name,
		String sample_type,
		String fusion_table
	)
	{
		int count = 0;
        Statement st = null;
        ResultSet rs = null;
        String expression = "";
        try {
        	
	        st = conn.createStatement();         // statement objects can be reused with
	
	        expression = 
	        	"SELECT COUNT(*)  " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE sample_name='"+ sample_name + "' " + 
	        	"AND sample_type='"+ sample_type + "'" +	        	
	        	"AND program='"+ fusion_program + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	
	        st.close();    

	        // if ResultSet is nulll returns -1
	        if(rs == null) 
	        {
	        	return count;
	        }
	        
			rs.next();
			count = rs.getInt(1);
			rs.close();
	    }
        catch (SQLException e) {
        	// TODO Auto-generated catch block
        	System.err.println("Error executing query:\n" + expression);
        	e.printStackTrace();
        }
		
		return count;
	}
	
	
	public void insertDefuseRecord(
			String defuse_record, 
			String sample_name, 
			String sample_type,
			String fusion_table
	)
	{
		String[] fusion_report_fields = defuse_record.split("\t");
		
		String program = "defuse";
		String tot_count = fusion_report_fields[58];
		String splitr_count = fusion_report_fields[2];
		String gene_chromosome1 = fusion_report_fields[26];
		String gene_chromosome2 = fusion_report_fields[27];
		String gene_start1 = fusion_report_fields[34];
		String gene_end1 = fusion_report_fields[28];
		String gene_start2 = fusion_report_fields[35];
		String gene_end2 = fusion_report_fields[29];
		String gene_strand1 = fusion_report_fields[36];
		String gene_strand2 = fusion_report_fields[37];
		String gene_name1 = fusion_report_fields[32];
		String gene_name2 = fusion_report_fields[33];
		String gene_id1 = fusion_report_fields[22];
		String gene_id2 = fusion_report_fields[23];
		String genomic_break_pos1 = fusion_report_fields[39];
		String genomic_break_pos2 = fusion_report_fields[40];
		
		String query = 
	        "INSERT INTO " + fusion_table + "(" +
	        "Sample_name," +
	        "program," +
	        "tot_count," +
	        "splitr_count," +
	        "gene_chromosome1," +
	        "gene_chromosome2," +
	        "gene_start1," +
	        "gene_end1," +
	        "gene_start2," +
	        "gene_end2," +
	        "gene_strand1," +
	        "gene_strand2," +
	        "gene_name1," +
	        "gene_name2," +
	        "gene_id1," +
	        "gene_id2," +
	        "genomic_break_pos1," +
	        "genomic_break_pos2," +
	        "sample_type" +
	        ") VALUES(" +
	        "'" + sample_name + "', " +
	        "'" + program + "', " +
	        "'" + tot_count + "', " +
	        "'" + splitr_count + "', " +
	        "'" + gene_chromosome1 + "', " +
	        "'" + gene_chromosome2 + "', " +
	        "'" + gene_start1 + "', " +
	        "'" + gene_end1 + "', " +
	        "'" + gene_start2 + "', " +
	        "'" + gene_end2 + "', " +
	        "'" + gene_strand1 + "', " +
	        "'" + gene_strand2 + "', " +
	        "'" + gene_name1 + "', " +
	        "'" + gene_name2 + "', " +
	        "'" + gene_id1 + "|', " +
	        "'" + gene_id2 + "|', " +
	        "'" + genomic_break_pos1 + "', " +
	        "'" + genomic_break_pos2 + "', " +
	        "'" + sample_type + "'" +
	        ")";
		
        try {
			this.update(query);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			System.err.println("Error executing query:\n" + query);
			e.printStackTrace();
		}
		
	}
	
	public void insertChimerascanRecord(
			String defuse_record, 
			String sample_name,
			String sample_type, 
			String fusion_table
	)
	{
		String[] fusion_report_fields = defuse_record.split("\t");
		
		String program = "chimerascan";
		String tot_count = fusion_report_fields[16];
		String splitr_count = fusion_report_fields[17];
		String gene_chromosome1 = fusion_report_fields[0].replaceAll("chr", "");
		//String tmp = gene_chromosome1.replaceAll("chr", "");
		String gene_chromosome2 = fusion_report_fields[3].replaceAll("chr", "");		
		//gene_chromosome2.replaceAll("/chr/", "");
		String gene_start1 = fusion_report_fields[1];
		String gene_end1 = fusion_report_fields[2];
		String gene_start2 = fusion_report_fields[4];
		String gene_end2 = fusion_report_fields[5];
		String gene_strand1 = fusion_report_fields[8];
		String gene_strand2 = fusion_report_fields[9];
		String gene_name1 = fusion_report_fields[12];
		String gene_name2 = fusion_report_fields[13];
		String gene_id1_field = fusion_report_fields[10];
		String gene_id2_field = fusion_report_fields[11];
		String genomic_break_pos1 = gene_strand1.compareTo("+")==0 ? gene_end1 : gene_start1;
		String genomic_break_pos2 = gene_strand2.compareTo("+")==0 ? gene_start2 : gene_end2;
		String gene_ids1 = "";
		String gene_ids2 = "";
		
		String[] gene_ids1_field = gene_id1_field.split(",");
		for (String gene_ids1_comma : gene_ids1_field)
		{
			gene_ids1 += gene_ids1_comma.split(":")[0] + "|";
		}

		String[] gene_ids2_field = gene_id2_field.split(",");
		for (String gene_ids2_comma : gene_ids2_field)
		{
			gene_ids2 += gene_ids2_comma.split(":")[0] + "|";
		}

		String query = 
	        "INSERT INTO " + fusion_table + "(" +
	        "Sample_name," +
	        "program," +
	        "tot_count," +
	        "splitr_count," +
	        "gene_chromosome1," +
	        "gene_chromosome2," +
	        "gene_start1," +
	        "gene_end1," +
	        "gene_start2," +
	        "gene_end2," +
	        "gene_strand1," +
	        "gene_strand2," +
	        "gene_name1," +
	        "gene_name2," +
	        "gene_id1," +
	        "gene_id2," +
	        "genomic_break_pos1," +
	        "genomic_break_pos2," +
	        "sample_type" +
	        ") VALUES(" +
	        "'" + sample_name + "', " +
	        "'" + program + "', " +
	        "'" + tot_count + "', " +
	        "'" + splitr_count + "', " +
	        "'" + gene_chromosome1 + "', " +
	        "'" + gene_chromosome2 + "', " +
	        "'" + gene_start1 + "', " +
	        "'" + gene_end1 + "', " +
	        "'" + gene_start2 + "', " +
	        "'" + gene_end2 + "', " +
	        "'" + gene_strand1 + "', " +
	        "'" + gene_strand2 + "', " +
	        "'" + gene_name1 + "', " +
	        "'" + gene_name2 + "', " +
	        "'" + gene_ids1.toString() + "', " +
	        "'" + gene_ids2.toString() + "', " +
	        "'" + genomic_break_pos1 + "', " +
	        "'" + genomic_break_pos2 + "', " +
	        "'" + sample_type + "'" +
	        ")";
		
        try {
			this.update(query);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			System.err.println("Error executing query:\n" + query);
			e.printStackTrace();
		}		
	}
	
	public void insertGeneralFusionRecord(
			String record, 
			String sample_name,
			String sample_type, 
			String fusion_table
	)
	{
		String[] fusion_report_fields = record.split("\t");
		
		String program = "general";
		String tot_count = fusion_report_fields[13];
		String splitr_count = fusion_report_fields[12];
		String gene_chromosome1 = fusion_report_fields[3].replaceAll("chr", "");
		String gene_chromosome2 = fusion_report_fields[7].replaceAll("chr", "");		
		String gene_start1 = fusion_report_fields[8];
		String gene_end1 = fusion_report_fields[9];
		String gene_start2 = fusion_report_fields[10];
		String gene_end2 = fusion_report_fields[11];
		String gene_strand1 = fusion_report_fields[2];
		String gene_strand2 = fusion_report_fields[6];
		String gene_name1 = fusion_report_fields[0];
		String gene_name2 = fusion_report_fields[4];
		String gene_id1_field = fusion_report_fields[1];
		String gene_id2_field = fusion_report_fields[5];
		String genomic_break_pos1 = gene_strand1.compareTo("+")==0 ? gene_end1 : gene_start1;
		String genomic_break_pos2 = gene_strand2.compareTo("+")==0 ? gene_start2 : gene_end2;
		String gene_ids1 = "";
		String gene_ids2 = "";
		
		String[] gene_ids1_field = gene_id1_field.split(",");
		for (String gene_ids1_comma : gene_ids1_field)
		{
			gene_ids1 += gene_ids1_comma.split(":")[0] + "|";
		}

		String[] gene_ids2_field = gene_id2_field.split(",");
		for (String gene_ids2_comma : gene_ids2_field)
		{
			gene_ids2 += gene_ids2_comma.split(":")[0] + "|";
		}

		String query = 
	        "INSERT INTO " + fusion_table + "(" +
	        "Sample_name," +
	        "program," +
	        "tot_count," +
	        "splitr_count," +
	        "gene_chromosome1," +
	        "gene_chromosome2," +
	        "gene_start1," +
	        "gene_end1," +
	        "gene_start2," +
	        "gene_end2," +
	        "gene_strand1," +
	        "gene_strand2," +
	        "gene_name1," +
	        "gene_name2," +
	        "gene_id1," +
	        "gene_id2," +
	        "genomic_break_pos1," +
	        "genomic_break_pos2," +
	        "sample_type" +
	        ") VALUES(" +
	        "'" + sample_name + "', " +
	        "'" + program + "', " +
	        "'" + tot_count + "', " +
	        "'" + splitr_count + "', " +
	        "'" + gene_chromosome1 + "', " +
	        "'" + gene_chromosome2 + "', " +
	        "'" + gene_start1 + "', " +
	        "'" + gene_end1 + "', " +
	        "'" + gene_start2 + "', " +
	        "'" + gene_end2 + "', " +
	        "'" + gene_strand1 + "', " +
	        "'" + gene_strand2 + "', " +
	        "'" + gene_name1 + "', " +
	        "'" + gene_name2 + "', " +
	        "'" + gene_ids1.toString() + "', " +
	        "'" + gene_ids2.toString() + "', " +
	        "'" + genomic_break_pos1 + "', " +
	        "'" + genomic_break_pos2 + "', " +
	        "'" + sample_type + "'" +
	        ")";
		
        try {
			this.update(query);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			System.err.println("Error executing query:\n" + query);
			e.printStackTrace();
		}		
	}
	
	public void insertBellerophontesRecord(
			String defuse_record, 
			String sample_name,
			String sample_type, 
			String fusion_table
	)
	{
		String[] fusion_report_fields = defuse_record.split("\t");
		
		String program = "bellerophontes";
		String tot_count = fusion_report_fields[14].replaceAll(" encompassing: ", "");
		String splitr_count = fusion_report_fields[13].replaceAll(" spanning: ", "");;
		String gene_chromosome1 = fusion_report_fields[3].replaceAll("chr", "");
		//String tmp = gene_chromosome1.replaceAll("chr", "");
		String gene_chromosome2 = fusion_report_fields[7].replaceAll("chr", "");		
		//gene_chromosome2.replaceAll("/chr/", "");
		String gene_start1 = fusion_report_fields[8];
		String gene_end1 = fusion_report_fields[9];
		String gene_start2 = fusion_report_fields[11];
		String gene_end2 = fusion_report_fields[12];
		String gene_strand1 = fusion_report_fields[2];
		String gene_strand2 = fusion_report_fields[6];
		String gene_name1 = fusion_report_fields[0];
		String gene_name2 = fusion_report_fields[4];
		String gene_id1_field = fusion_report_fields[1];
		String gene_id2_field = fusion_report_fields[5];
		String genomic_break_pos1 = gene_strand1.compareTo("+")==0 ? gene_end1 : gene_start1;
		String genomic_break_pos2 = gene_strand2.compareTo("+")==0 ? gene_start2 : gene_end2;
		String gene_ids1 = "";
		String gene_ids2 = "";
		
		String[] gene_ids1_field = gene_id1_field.split(",");
		for (String gene_ids1_comma : gene_ids1_field)
		{
			gene_ids1 += gene_ids1_comma.split(":")[0] + "|";
		}

		String[] gene_ids2_field = gene_id2_field.split(",");
		for (String gene_ids2_comma : gene_ids2_field)
		{
			gene_ids2 += gene_ids2_comma.split(":")[0] + "|";
		}

		String query = 
	        "INSERT INTO " + fusion_table + "(" +
	        "Sample_name," +
	        "program," +
	        "tot_count," +
	        "splitr_count," +
	        "gene_chromosome1," +
	        "gene_chromosome2," +
	        "gene_start1," +
	        "gene_end1," +
	        "gene_start2," +
	        "gene_end2," +
	        "gene_strand1," +
	        "gene_strand2," +
	        "gene_name1," +
	        "gene_name2," +
	        "gene_id1," +
	        "gene_id2," +
	        "genomic_break_pos1," +
	        "genomic_break_pos2," +
	        "sample_type" +
	        ") VALUES(" +
	        "'" + sample_name + "', " +
	        "'" + program + "', " +
	        "'" + tot_count + "', " +
	        "'" + splitr_count + "', " +
	        "'" + gene_chromosome1 + "', " +
	        "'" + gene_chromosome2 + "', " +
	        "'" + gene_start1 + "', " +
	        "'" + gene_end1 + "', " +
	        "'" + gene_start2 + "', " +
	        "'" + gene_end2 + "', " +
	        "'" + gene_strand1 + "', " +
	        "'" + gene_strand2 + "', " +
	        "'" + gene_name1 + "', " +
	        "'" + gene_name2 + "', " +
	        "'" + gene_ids1.toString() + "', " +
	        "'" + gene_ids2.toString() + "', " +
	        "'" + genomic_break_pos1 + "', " +
	        "'" + genomic_break_pos2 + "', " +
	        "'" + sample_type + "'" +
	        ")";
		
        try {
			this.update(query);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			System.err.println("Error executing query:\n" + query);
			e.printStackTrace();
		}		
	}
	
	public boolean IsSampleAlreadyPresent(String sample_name, String program, String fusion_table)
	{
        Statement st = null;
        ResultSet rs = null;
        
        boolean results = true;

        try {
        	
	        st = conn.createStatement();         // statement objects can be reused with
	
	        String expression = 
	        	"SELECT DISTINCT(sample_name) " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE sample_name='"+ sample_name + "' " + 
	        	"AND program='"+ program + "'";
	        
	        // repeated calls to execute but we
	        // choose to make a new one each time
	        rs = st.executeQuery(expression);    // run the query
	
	        // do something with the result set.
 			if(!rs.next())
			{
				results = false;
			}
 			else
 			{
 				results = true;
 			}
	        st.close();    
        // NOTE!! if you close a statement the associated ResultSet is
        // closed too
        // so you should copy the contents to some other object.
        // the result set is invalidated also  if you recycle an Statement
        // and try to execute some other query before the result set has been
        // completely examined.
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return results;
	}
	
    public void shutdown() throws SQLException {

        Statement st = conn.createStatement();

        // db writes out to files and performs clean shuts down
        // otherwise there will be an unclean shutdown
        // when program ends
        st.execute("SHUTDOWN");
        conn.close();    // if there are no other open connection
    }
    
    public synchronized void update(String expression) throws SQLException {

        Statement st = null;

        st = conn.createStatement();    // statements

        int i = st.executeUpdate(expression);    // run the query

        if (i == -1) {
            System.out.println("db error : " + expression);
        }

        st.close();
    }

	
	public void FusionDB_start(String database_path, String fusion_db)
	{
		server = new Server();
		server.setLogWriter(null);
		server.setSilent(true);
		server.setDatabasePath(0, "file:" + database_path);
		server.setDatabaseName(0, fusion_db);
		server.start();
	}
	
    
    public void ConnectFusionDB(String DBAlias)
    {
		try 
		{
			Class.forName("org.hsqldb.jdbcDriver");
			conn = DriverManager.getConnection("jdbc:hsqldb:hsql://localhost/" + DBAlias, "SA", "");
		} 
		catch (ClassNotFoundException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
    }
    
    public void CloseConnectionFusionDB()
    {
		try 
		{
			//conn.prepareStatement("shutdown").execute();
			String sql = "SHUTDOWN";
			Statement stmt = conn.createStatement();
			stmt.executeUpdate(sql);
			stmt.close();
	    	conn.close();
		} 
		catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
    }

	public void FusionDB_stop()
	{
		server.stop();
	}

    
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String fusion_input_report_file_path = "";
		String fusion_program_type = "";
		String db_filename_prefix = "";
		String fusion_table = "";
		String sample_info = "";
		String db_info = "";
        Getopt g = new Getopt("LoadFusionReport", args, "i:p:d:s:t:f:");
        //
        int c;

        while ((c = g.getopt()) != -1)
        {
            switch(c)
            {
            	case 'i':
            		fusion_input_report_file_path = g.getOptarg();
            		break;

            	case 'f':
            		db_info = g.getOptarg();
            		break;

                 case 'p':
                	 fusion_program_type = g.getOptarg();
                	 break;

                 case 'd':
                	 db_filename_prefix = g.getOptarg();
                	 break;

                 case 't':
                	 fusion_table = g.getOptarg();
                	 break;

                 case 's':
                	 sample_info = g.getOptarg();
                	 break;

                 case '?':
                	 System.exit(1);
                	 break;

                 default:
                	 System.exit(1);
                	 System.out.println("default");
            }
        }
		
        if (	
        		fusion_input_report_file_path.compareTo("")==0 || 
        		fusion_table.compareTo("")==0 || 
        		fusion_program_type.compareTo("")==0
        	)
        {
        	System.err.println("Error!\nusage: java -jar LoadFusionReport " +
        			"\n\t\t-i [ Fusion Report Input File ]" +
        			"\n\t\t-d [ DB filename prefix ]" +
        			"\n\t\t-f [ path to db ]" +
        			"\n\t\t-t [ Fusion Table ]" +
        			"\n\t\t-s [ Sample Info (sample_name|sample_type) ]" +
        			"\n\t\t-p [ Fusion Program ]");
        	System.exit(1);
        }
        
        String[] sample_info_fields = sample_info.split("\\|");
    	LoadFusionReport lfr = new LoadFusionReport();
        try
        {
        	lfr.FusionDB_start(db_info, "xdb");
        	lfr.ConnectFusionDB("xdb");
        	lfr.run(fusion_input_report_file_path, sample_info_fields[0], sample_info_fields[1], fusion_program_type, db_filename_prefix, fusion_table);
        	lfr.CloseConnectionFusionDB();
        	lfr.FusionDB_stop();
        }
        catch (Exception e) {
        	e.printStackTrace();
			// TODO: handle exception
		}
        finally
        {
        	lfr.FusionDB_stop();
        	System.exit(0);
        }

	}

}
