package QueryFusionDatabase;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;

import org.hsqldb.Server;

import gnu.getopt.Getopt;

public class QueryFusionDatabase {
	private Server server;
	private Connection conn;
	
	public QueryFusionDatabase()
	{
//		Class.forName("org.hsqldb.jdbcDriver");
//		
//		conn = DriverManager.getConnection("jdbc:hsqldb:hsql://localhost/xdb", "SA", "");

	}
	
	public void run(String command,
			String fusion_table,
			String sample_type
	) throws ClassNotFoundException, SQLException
	{
		if(command.equals("listsamples"))
		{
			ArrayList<String> samples = getSampleAvailable(fusion_table);
			
			PrintHeading();
			
			for (String sample : samples)
			{
				//getFusionListFromSample();
				
				ArrayList<String> fusion_ids = getFusionIDsBySample(sample, fusion_table);
				
				for(String fusion_id : fusion_ids)
				{
					Print2ProgOverlappingSampleListByFusionID(fusion_id, fusion_table);
				}
			}
		}
		
		else if(command.equals("deleteAll"))
		{
			deleteAllFusions(fusion_table);
		}
		
		else if(command.equals("kinases"))
		{
			if(sample_type.compareTo("")==0)
			{
				System.err.println("Error: with kinases command Sample_type is mandatory!");
			}
			

			ArrayList<String> fusion_ids = getFusionIDsBySample_Type(sample_type, fusion_table);
				
			for(String fusion_id : fusion_ids)
			{
				PrintKinasesFusion(fusion_id, fusion_table);
			}			
		}

		else if(command.equals("kinases_cf"))
		{
			if(sample_type.compareTo("")==0)
			{
				System.err.println("Error: with kinases_cf command Sample_type is mandatory!");
			}
			

			ArrayList<String> fusion_ids = getFusionIDsBySample_Type(sample_type, fusion_table);
				
			for(String fusion_id : fusion_ids)
			{
				PrintKinasesFullFormatNoProgramOverlap(fusion_id, fusion_table);
			}			
		}

		else if(command.equals("all_with_kinases"))
		{
			if(sample_type.compareTo("")==0)
			{
				System.err.println("Error: with kinases_cf command Sample_type is mandatory!");
			}
			

			ArrayList<String> fusion_ids = getFusionIDsBySample_Type(sample_type, fusion_table);
				
			for(String fusion_id : fusion_ids)
			{
				PrintAllWithInfoKinaseNoProgramOverlap(fusion_id, fusion_table);
			}			
		}

		else if(command.equals("kinases_cf_TP"))
		{
			if(sample_type.compareTo("")==0)
			{
				System.err.println("Error: with kinases_cf_TP command Sample_type is mandatory!");
			}
			

			ArrayList<String> fusion_ids = getFusionIDsBySample_Type(sample_type, fusion_table);
				
			for(String fusion_id : fusion_ids)
			{
				PrintThreePrimeKinasesFullFormatNoProgramOverlap(fusion_id, fusion_table);
			}			
		}
		else
		{
			System.err.println("Error: Unknown Command.");
		}


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
	
    public synchronized void update(String expression) throws SQLException {

        Statement st = null;

        st = conn.createStatement();    // statements

        int i = st.executeUpdate(expression);    // run the query

        if (i == -1) {
            System.out.println("db error : " + expression);
        }

        st.close();
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
		//server.shutdown();
	}
	
	public void PrintHeading()
	{
		StringBuffer sb = new StringBuffer();
		
		sb.append("FUSION_ID\t");
		sb.append("SAMPLE_NAME\t");
		sb.append("PROGRAM\t");
		sb.append("TOT_COUNT\t");
		sb.append("SPLITR_COUNT\t");
		sb.append("GENE_CHROMOSOME1\t");
		sb.append("GENE_CHROMOSOME2\t");
		sb.append("GENE_START1\t");
		sb.append("GENE_END1\t");
		sb.append("GENE_START2\t");
		sb.append("GENE_END2\t");
		sb.append("GENE_STRAND1\t");
		sb.append("GENE_STRAND2\t");
		sb.append("GENE_NAME1\t");
		sb.append("GENE_NAME2\t");
		sb.append("GENOMIC_BREAK_POS1\t");
		sb.append("GENOMIC_BREAK_POS2\t");
		sb.append("SAMPLE_TYPE\t");
		sb.append("SAMPLE_NAME_LIST\t");
		sb.append("SAMPLE_TYPE_LIST\t");
		sb.append("Kinase\n");
		
		System.out.print(sb.toString());
	}
	
	public void Print2ProgOverlappingSampleListByFusionID(String fusion_id, String fusion_table)
	{
		ArrayList<String> sample_types_defuse = new ArrayList<String>();
		ArrayList<String> sample_names_defuse = new ArrayList<String>();
		ArrayList<String> sample_types_chimerascan = new ArrayList<String>();
		ArrayList<String> sample_names_chimerascan = new ArrayList<String>();

        Statement st = null;
        ResultSet rs = null;
        
        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT * " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE fusion_id='"+ fusion_id + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	        
	        if(rs.next())
	        {
	        	//String program = rs.getString("program");
	        	String gene_chromosome1 = rs.getString("gene_chromosome1");
	        	String gene_chromosome2 = rs.getString("gene_chromosome2");
	        	int genomic_break_pos1 = Integer.parseInt(rs.getString("genomic_break_pos1"));
	        	int genomic_break_pos2 = Integer.parseInt(rs.getString("genomic_break_pos2"));	    
	        	
	        	ResultSetMetaData md = rs.getMetaData();
        		int col = md.getColumnCount();
        		
		        StringBuffer sb = new StringBuffer();

		        for (int i = 1; i <= col; i++){
		        	String col_name = md.getColumnName(i);
		        	sb.append(rs.getString(col_name) + "\t");
	        	}
	        	
		        String nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "defuse" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
		        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_defuse.add(rs.getString("sample_type"));
		        	sample_names_defuse.add(rs.getString("sample_name"));
		        }
		        
		        nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "chimerascan" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
			        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_chimerascan.add(rs.getString("sample_type"));
		        	sample_names_chimerascan.add(rs.getString("sample_name"));
		        } 
		        
		        if(
		        		!sample_names_chimerascan.isEmpty() && 
		        		!sample_names_defuse.isEmpty() && 
		        		!sample_types_chimerascan.isEmpty() && 
		        		!sample_types_defuse.isEmpty()
		        )
		        {
		        	TreeMap<String, String> uniq_types = new TreeMap<String, String>();
		        	TreeMap<String, String> uniq_sample_names = new TreeMap<String, String>();
		        	
		        	System.out.print(sb.toString());
		        	for(String sample_name : sample_names_chimerascan)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : sample_names_defuse)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : uniq_sample_names.keySet())
		        	{
		        		System.out.print(sample_name + ",");
		        	}
		        	System.out.print("\t");
		        	
		        	
		        	for(String sample_type : sample_types_chimerascan)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : sample_types_defuse)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : uniq_types.keySet())
		        	{
		        		System.out.print(sample_type + ",");
		        	}
		        	
		        	System.out.print("\t");
		        	System.out.print("Kinase_Unspecified\t");
		        	System.out.print("\n");
		        }

	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		//return samples;

	}
	
	public void PrintThreePrimeKinasesFullFormatNoProgramOverlap(String fusion_id, String fusion_table)
	{
		ArrayList<String> sample_types_defuse = new ArrayList<String>();
		ArrayList<String> sample_names_defuse = new ArrayList<String>();
		ArrayList<String> sample_types_chimerascan = new ArrayList<String>();
		ArrayList<String> sample_names_chimerascan = new ArrayList<String>();

        Statement st = null;
        ResultSet rs = null;
        
        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT * " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE fusion_id='"+ fusion_id + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	        
	        if(rs.next())
	        {
	        	String program = rs.getString("program");
	        	//String gene_id1 = rs.getString("GENE_ID1");
	        	String gene_id2 = rs.getString("GENE_ID2");
	        	String gene_chromosome1 = rs.getString("gene_chromosome1");
	        	String gene_chromosome2 = rs.getString("gene_chromosome2");
	        	int genomic_break_pos1 = Integer.parseInt(rs.getString("genomic_break_pos1"));
	        	int genomic_break_pos2 = Integer.parseInt(rs.getString("genomic_break_pos2"));	    
	        	
	        	ResultSetMetaData md = rs.getMetaData();
        		int col = md.getColumnCount();
        		
		        StringBuffer sb = new StringBuffer();

		        for (int i = 1; i <= col; i++){
		        	String col_name = md.getColumnName(i);
		        	sb.append(rs.getString(col_name) + "\t");
	        	}
	        	
		        String nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "defuse" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
		        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_defuse.add(rs.getString("sample_type"));
		        	sample_names_defuse.add(rs.getString("sample_name"));
		        }
		        
		        nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "chimerascan" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
			        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_chimerascan.add(rs.getString("sample_type"));
		        	sample_names_chimerascan.add(rs.getString("sample_name"));
		        } 
		        
		        //kinases check
	        	//HashMap<String, String> kinases_G1 = getKinasesList(program, gene_id1);
	        	HashMap<String, String> kinases_G2 = getKinasesList(program, gene_id2);
		        
		        if(
		        		kinases_G2.size()>0
		        )
		        {
		        	TreeMap<String, String> uniq_types = new TreeMap<String, String>();
		        	TreeMap<String, String> uniq_sample_names = new TreeMap<String, String>();
		        	
		        	System.out.print(sb.toString());
		        	for(String sample_name : sample_names_chimerascan)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : sample_names_defuse)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : uniq_sample_names.keySet())
		        	{
		        		System.out.print(sample_name + ",");
		        	}
		        	System.out.print("\t");
		        	
		        	
		        	for(String sample_type : sample_types_chimerascan)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : sample_types_defuse)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : uniq_types.keySet())
		        	{
		        		System.out.print(sample_type + ",");
		        	}
		        	
		        	
		        	System.out.print("\tKINASE\n");
		        }

	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		//return samples;

	}


	public void PrintKinasesFullFormatNoProgramOverlap(String fusion_id, String fusion_table)
	{
		ArrayList<String> sample_types_defuse = new ArrayList<String>();
		ArrayList<String> sample_names_defuse = new ArrayList<String>();
		ArrayList<String> sample_types_chimerascan = new ArrayList<String>();
		ArrayList<String> sample_names_chimerascan = new ArrayList<String>();

        Statement st = null;
        ResultSet rs = null;
        
        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT * " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE fusion_id='"+ fusion_id + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	        
	        if(rs.next())
	        {
	        	String program = rs.getString("program");
	        	String gene_id1 = rs.getString("GENE_ID1");
	        	String gene_id2 = rs.getString("GENE_ID2");
	        	String gene_chromosome1 = rs.getString("gene_chromosome1");
	        	String gene_chromosome2 = rs.getString("gene_chromosome2");
	        	int genomic_break_pos1 = Integer.parseInt(rs.getString("genomic_break_pos1"));
	        	int genomic_break_pos2 = Integer.parseInt(rs.getString("genomic_break_pos2"));	    
	        	
	        	ResultSetMetaData md = rs.getMetaData();
        		int col = md.getColumnCount();
        		
		        StringBuffer sb = new StringBuffer();

		        for (int i = 1; i <= col; i++){
		        	String col_name = md.getColumnName(i);
		        	sb.append(rs.getString(col_name) + "\t");
	        	}
	        	
		        String nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "defuse" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
		        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_defuse.add(rs.getString("sample_type"));
		        	sample_names_defuse.add(rs.getString("sample_name"));
		        }
		        
		        nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "chimerascan" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
			        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_chimerascan.add(rs.getString("sample_type"));
		        	sample_names_chimerascan.add(rs.getString("sample_name"));
		        } 
		        
		        //kinases check
	        	HashMap<String, String> kinases_G1 = getKinasesList(program, gene_id1);
	        	HashMap<String, String> kinases_G2 = getKinasesList(program, gene_id2);
		        
		        if(
		        		kinases_G1.size()>0 || 
		        		kinases_G2.size()>0
		        )
		        {
		        	TreeMap<String, String> uniq_types = new TreeMap<String, String>();
		        	TreeMap<String, String> uniq_sample_names = new TreeMap<String, String>();
		        	
		        	System.out.print(sb.toString());
		        	for(String sample_name : sample_names_chimerascan)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : sample_names_defuse)
		        	{
		        		uniq_sample_names.put(sample_name, "");
		        	}
		        	for(String sample_name : uniq_sample_names.keySet())
		        	{
		        		System.out.print(sample_name + ",");
		        	}
		        	System.out.print("\t");
		        	
		        	
		        	for(String sample_type : sample_types_chimerascan)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : sample_types_defuse)
		        	{
		        		uniq_types.put(sample_type, "");
		        	}
		        	for(String sample_type : uniq_types.keySet())
		        	{
		        		System.out.print(sample_type + ",");
		        	}
		        	
		        	
		        	System.out.print("\tKINASE\n");
		        }

	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		//return samples;

	}

	public void PrintAllWithInfoKinaseNoProgramOverlap(String fusion_id, String fusion_table)
	{
		ArrayList<String> sample_types_defuse = new ArrayList<String>();
		ArrayList<String> sample_names_defuse = new ArrayList<String>();
		ArrayList<String> sample_types_bellerophontes = new ArrayList<String>();
		ArrayList<String> sample_names_bellerophontes = new ArrayList<String>();
		ArrayList<String> sample_types_chimerascan = new ArrayList<String>();
		ArrayList<String> sample_names_chimerascan = new ArrayList<String>();

        Statement st = null;
        ResultSet rs = null;
        
        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT * " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE fusion_id='"+ fusion_id + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	        
	        if(rs.next())
	        {
	        	String program = rs.getString("program");
	        	String gene_id1 = rs.getString("GENE_ID1");
	        	String gene_id2 = rs.getString("GENE_ID2");
	        	String gene_chromosome1 = rs.getString("gene_chromosome1");
	        	String gene_chromosome2 = rs.getString("gene_chromosome2");
	        	int genomic_break_pos1 = Integer.parseInt(rs.getString("genomic_break_pos1"));
	        	int genomic_break_pos2 = Integer.parseInt(rs.getString("genomic_break_pos2"));	    
	        	
	        	ResultSetMetaData md = rs.getMetaData();
        		int col = md.getColumnCount();
        		
		        StringBuffer sb = new StringBuffer();

		        for (int i = 1; i <= col; i++){
		        	String col_name = md.getColumnName(i);
		        	sb.append(rs.getString(col_name) + "\t");
	        	}
	        	
		        String nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "defuse" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
		        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_defuse.add(rs.getString("sample_type"));
		        	sample_names_defuse.add(rs.getString("sample_name"));
		        }
		        
		        nested_expression = 
			        	"SELECT sample_name, sample_type " +
			        	"FROM " + fusion_table + " " +
			        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
			        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
			        	"AND program='"+ "bellerophontes" + "'" +
			        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
			        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
			        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
			        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
			        
		        rs = st.executeQuery(nested_expression);
			        
		        while(rs.next())
		        {
		        	sample_types_bellerophontes.add(rs.getString("sample_type"));
		        	sample_names_bellerophontes.add(rs.getString("sample_name"));
		        }
		        
		        nested_expression = 
		        	"SELECT sample_name, sample_type " +
		        	"FROM " + fusion_table + " " +
		        	"WHERE gene_chromosome1='"+ gene_chromosome1 + "'" +
		        	"AND gene_chromosome2='"+ gene_chromosome2 + "'" +
		        	"AND program='"+ "chimerascan" + "'" +
		        	"AND genomic_break_pos1>'"+ (genomic_break_pos1-5) + "'" +
		        	"AND genomic_break_pos1<'"+ (genomic_break_pos1+5) + "'" +
		        	"AND genomic_break_pos2>'"+ (genomic_break_pos2-5) + "'" +
		        	"AND genomic_break_pos2<'"+ (genomic_break_pos2+5) + "'";
			        
		        rs = st.executeQuery(nested_expression);
		        
		        while(rs.next())
		        {
		        	sample_types_chimerascan.add(rs.getString("sample_type"));
		        	sample_names_chimerascan.add(rs.getString("sample_name"));
		        } 
		        
		        //kinases check
	        	HashMap<String, String> kinases_G1 = getKinasesList(program, gene_id1);
	        	HashMap<String, String> kinases_G2 = getKinasesList(program, gene_id2);
		        

	        	TreeMap<String, String> uniq_types = new TreeMap<String, String>();
	        	TreeMap<String, String> uniq_sample_names = new TreeMap<String, String>();
	        	
	        	System.out.print(sb.toString());
	        	for(String sample_name : sample_names_chimerascan)
	        	{
	        		uniq_sample_names.put(sample_name, "");
	        	}
	        	for(String sample_name : sample_names_bellerophontes)
	        	{
	        		uniq_sample_names.put(sample_name, "");
	        	}
	        	for(String sample_name : sample_names_defuse)
	        	{
	        		uniq_sample_names.put(sample_name, "");
	        	}
	        	for(String sample_name : uniq_sample_names.keySet())
	        	{
	        		System.out.print(sample_name + ",");
	        	}
	        	System.out.print("\t");
	        	
	        	
	        	for(String sample_type : sample_types_chimerascan)
	        	{
	        		uniq_types.put(sample_type, "");
	        	}
	        	for(String sample_type : sample_types_bellerophontes)
	        	{
	        		uniq_types.put(sample_type, "");
	        	}
	        	for(String sample_type : sample_types_defuse)
	        	{
	        		uniq_types.put(sample_type, "");
	        	}
	        	for(String sample_type : uniq_types.keySet())
	        	{
	        		System.out.print(sample_type + ",");
	        	}
		        	
		        if(
		        		kinases_G1.size()>0 && 
		        		kinases_G2.size()>0
		        )
		        {
		        	System.out.print("\tBOTH_KINASE\n");
		        }
		        else if(kinases_G1.size()>0)
		        {
		        	System.out.print("\t5p_KINASE\n");
		        }
		        else if(kinases_G2.size()>0)
		        {
		        	System.out.print("\t3p_KINASE\n");
		        }
		        else
		        {
		        	System.out.print("\tNO_KINASE\n");
		        }

	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		//return samples;

	}

	
	public HashMap<String, String> getKinasesList(String program, String geneID)
	{
        Statement st = null;
        ResultSet rs = null;
        
    	HashMap<String, String> kinases = new HashMap<String, String>();
        
        try {
        	
	        st = conn.createStatement();
	        
        	String[] gene_id_list = geneID.split("\\|");
        	HashMap<String, String> gene_id_list_map = new HashMap<String, String>();
        	for (String s : gene_id_list)
        	{
        		gene_id_list_map.put(s.split("\\.")[0], null);
        	}
        	
        	String nested_expression;
        	if (program.compareTo("chimerascan")==0)
        	{
        		nested_expression = 
        			"SELECT * FROM ENSEMBL_UCSC_ID, KINASE_ENSEMBL_ID_SYM, ENSEMBL_GENE_TRANSCRIPT_ID " +
        			"WHERE KINASE_ENSEMBL_ID_SYM.ID=ENSEMBL_GENE_TRANSCRIPT_ID.GENE_ID " +
        			"AND ENSEMBL_GENE_TRANSCRIPT_ID.TRANSCRIPT_ID=ENSEMBL_UCSC_ID.ENSEMBL_ID AND (";
        	
        		Iterator<String> s_iter = gene_id_list_map.keySet().iterator();
        		if(s_iter.hasNext())
        		{
        			nested_expression += "ENSEMBL_UCSC_ID.UCSC_ID like '" + s_iter.next() + "%'";
        		}
        		while (s_iter.hasNext())
        		{
        			nested_expression += "OR ENSEMBL_UCSC_ID.UCSC_ID like '" + s_iter.next() + "%'";
        		}
        		nested_expression += ");";
        	}
        	else
        	{
        		nested_expression = "SELECT * FROM KINASE_ENSEMBL_ID_SYM ";
        		Iterator<String> s_iter = gene_id_list_map.keySet().iterator();
        		if(s_iter.hasNext())
        		{
        			nested_expression += "where ID='" + s_iter.next() + "'";
        		}
        		while (s_iter.hasNext())
        		{
        			nested_expression += "OR ID='" + s_iter.next() + "'";
        		}
        	}
        	rs = st.executeQuery(nested_expression);
	        
	        while(rs.next())
	        {
	        	kinases.put(rs.getString("SYMBOL"), null);
	        }
        	
	        
	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return kinases;
	}
	
	public void PrintKinasesFusion(String fusion_id, String fusion_table)
	{

        Statement st = null;
        ResultSet rs = null;
        
        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT GENE_ID1,GENE_ID2,GENE_NAME1,GENE_NAME2, PROGRAM " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE fusion_id='"+ fusion_id + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	        
	        while(rs.next())
	        {
	        	String program = rs.getString("program");
	        	String gene_id1 = rs.getString("GENE_ID1");
	        	String gene_id2 = rs.getString("GENE_ID2");
	        	String gene_name1 = rs.getString("GENE_NAME1");
	        	String gene_name2 = rs.getString("GENE_NAME2");	  
	        	HashMap<String, String> kinases_G1 = getKinasesList(program, gene_id1);
	        	HashMap<String, String> kinases_G2 = getKinasesList(program, gene_id2);
	        	
	        	
	        	if(kinases_G1.size()>0 || kinases_G2.size()>0)
	        	{
	        		System.out.print(fusion_id + "\t" + gene_name1 + "\t" + gene_name2 + "\t");
	        		if(kinases_G1.size()>0)
	        		{
	        			for (String g_id : kinases_G1.keySet())
	        			{
	        				System.out.print(g_id + ",");
	        			}
	        		}
	        		else
	        		{
	        			System.out.print("NO_KINASE");
	        		}
	        		System.out.print("\t");
	        	
	        		if(kinases_G2.size()>0)
	        		{
	        			for (String g_id : kinases_G2.keySet())
	        			{
	        				System.out.print(g_id + ",");
	        			}
	        		}
	        		else
	        		{
	        			System.out.print("NO_KINASE");
	        		}
	        		System.out.print("\n");
	        	}

	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		//return samples;

	}

	
	public ArrayList<String> getFusionIDsBySample_Type(String sample_type, String fusion_table)
	{
		ArrayList<String> fusion_ids = new ArrayList<String>();
		
        Statement st = null;
        ResultSet rs = null;

        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT fusion_id, " +
	        	"program, gene_chromosome1, " +
	        	"gene_chromosome2, genomic_break_pos1, " +
	        	"genomic_break_pos2 " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE sample_type='"+ sample_type + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	
	        while (rs.next())
	        {
	        	String fusion_id = rs.getString("fusion_id");
	        	
	        	fusion_ids.add(fusion_id);
	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		return fusion_ids;
	}
	
	
	
	
	
	public ArrayList<String> getFusionIDsBySample(String sample, String fusion_table)
	{
		ArrayList<String> fusion_ids = new ArrayList<String>();
		
        Statement st = null;
        ResultSet rs = null;

        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT fusion_id, " +
	        	"program, gene_chromosome1, " +
	        	"gene_chromosome2, genomic_break_pos1, " +
	        	"genomic_break_pos2 " +
	        	"FROM " + fusion_table + " " +
	        	"WHERE sample_name='"+ sample + "'";
	        
	        rs = st.executeQuery(expression);    // run the query
	
	        while (rs.next())
	        {
	        	String fusion_id = rs.getString("fusion_id");
	        	
	        	fusion_ids.add(fusion_id);
	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		return fusion_ids;
	}
	
	public ArrayList<String> deleteAllFusions(String fusion_table)
	{
		ArrayList<String> samples = new ArrayList<String>();
        Statement st = null;
        ResultSet rs = null;

        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"DELETE FROM " + fusion_table + " ";
	        
	        st.executeUpdate(expression);    // run the query
	
	        st.close();
	        
	        st = conn.createStatement();
	    	
	        expression = 
	        	"SELECT * FROM " + fusion_table + " ";
	        
	        rs = st.executeQuery(expression);    // run the query
	        //rs.last();
	        int count = rs.getRow();
	        if (count!=0)
	        {
	        		System.err.printf("Error: not all the fusions have been properly deleted\n");
	        }
	        else
	        {
	        	System.err.printf("All fusions have been properly deleted\n");
	        	
	        }

	        st.close();    
        
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		return samples;
	}
	
	public ArrayList<String> getSampleAvailable(String fusion_table)
	{
		ArrayList<String> samples = new ArrayList<String>();
        Statement st = null;
        ResultSet rs = null;

        try {
        	
	        st = conn.createStatement();
	
	        String expression = 
	        	"SELECT distinct(Sample_name)" +
	        	"FROM " + fusion_table + " ";
	        
	        rs = st.executeQuery(expression);    // run the query
	
	        while (rs.next())
	        {
	        	samples.add(rs.getString("SAMPLE_NAME"));
	        }

	        st.close();    
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		return samples;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String command = "";
		String fusion_table = "";
		String sample_type = "";
		String db_info = "";
		Getopt g = new Getopt("QueryFusionDatabase", args, "c:t:s:d:");
        //
        int c;

        while ((c = g.getopt()) != -1)
        {
            switch(c)
            {
	        	case 'd':
	        		db_info = g.getOptarg();
	        		break;
	        		
	        	case 'c':
	        		command = g.getOptarg();
	        		break;
	        		
            	case 's':
            		sample_type = g.getOptarg();
            		break;

            	case 't':
            		fusion_table = g.getOptarg();
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
        		command.compareTo("")==0 ||
        		fusion_table.compareTo("")==0
        	)
        {
        	System.err.println("Error!\nusage: java -jar QueryFusionDatabase " +
        			"\n\t\t-t [ Fusion Table ]" +
        			"\n\t\t-c [ Command Option ]" +
        			"\n\t\t-c [ deleteAll ]" +
        			"\n\t\t\t [ listsamples ]" +
        			"\n\t\t\t [ kinases ]" +
        			"\n\t\t\t [ kinases_cf ]" +
        			"\n\t\t\t [ all_with_kinases ]" +
        			"\n\t\t\t [ kinases_cf_TP ]" +
        			"");
        	System.exit(1);
        }
        
        QueryFusionDatabase qfd = new QueryFusionDatabase();
        
        try
        {
        	
        	
        	//qfd.FusionDB_start("db/mydb", "xdb");
        	qfd.FusionDB_start(db_info, "xdb");
        	qfd.ConnectFusionDB("xdb");
        	qfd.run(command, fusion_table, sample_type);
        	qfd.CloseConnectionFusionDB();
        	//qfd.FusionDB_stop();
        	
        	System.err.print("Completed.\n");
        }
        catch (Exception e) {
        	e.printStackTrace();
        	qfd.FusionDB_stop();
			// TODO: handle exception
		}
        finally
        {
        	qfd.FusionDB_stop();
        	System.exit(0);
        }
	}

}
