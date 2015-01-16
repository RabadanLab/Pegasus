package example;

import gnu.getopt.Getopt;

import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.logging.Logger;


public class UniprotWebServiceClient
{
    private static final String	UNIPROT_SERVER	= "http://www.uniprot.org/";
    private static final Logger	LOG		= Logger.getAnonymousLogger();

    private static void run(String tool, ParameterNameValue[] params)
        throws Exception
    {
        StringBuilder locationBuilder = new StringBuilder(UNIPROT_SERVER + tool + "/?");
        for (int i = 0; i < params.length; i++)
        {
            if (i > 0)
                locationBuilder.append('&');
            locationBuilder.append(params[i].name).append('=').append(params[i].value);
        }
        String location = locationBuilder.toString();
        URL url = new URL(location);
        LOG.info("Submitting...");
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        HttpURLConnection.setFollowRedirects(true);
        conn.setDoInput(true);
        conn.connect();

        int status = conn.getResponseCode();
        while (true)
        {
            int wait = 0;
            String header = conn.getHeaderField("Retry-After");
            if (header != null)
                wait = Integer.valueOf(header);
            if (wait == 0)
                break;
            LOG.info("Waiting (" + wait + ")...");
            conn.disconnect();
            Thread.sleep(wait * 1000);
            conn = (HttpURLConnection) new URL(location).openConnection();
            conn.setDoInput(true);
            conn.connect();
            status = conn.getResponseCode();
        }
        if (status == HttpURLConnection.HTTP_OK)
        {
            LOG.info("Got a OK reply");

            InputStream reader = conn.getInputStream();
            URLConnection.guessContentTypeFromStream(reader);
            StringBuilder builder = new StringBuilder();
            int a = 0;
            while ((a = reader.read()) != -1)
            {
                builder.append((char) a);
            }
            System.out.println(builder.toString());
        }
        else
            LOG.severe("Failed, got " + conn.getResponseMessage() + " for "
                + location);
        conn.disconnect();
    }

    public String run_public(String tool, ParameterNameValue[] params)
    throws Exception
{
    	StringBuffer sb = new StringBuffer();
    StringBuilder locationBuilder = new StringBuilder(UNIPROT_SERVER + tool + "/?");
    for (int i = 0; i < params.length; i++)
    {
        if (i > 0)
            locationBuilder.append('&');
        locationBuilder.append(params[i].name).append('=').append(params[i].value);
    }
    String location = locationBuilder.toString();
    URL url = new URL(location);
    LOG.info("Submitting...");
    HttpURLConnection conn = (HttpURLConnection) url.openConnection();
    HttpURLConnection.setFollowRedirects(true);
    conn.setDoInput(true);
    conn.connect();

    int status = conn.getResponseCode();
    while (true)
    {
        int wait = 0;
        String header = conn.getHeaderField("Retry-After");
        if (header != null)
            wait = Integer.valueOf(header);
        if (wait == 0)
            break;
        LOG.info("Waiting (" + wait + ")...");
        conn.disconnect();
        Thread.sleep(wait * 1000);
        conn = (HttpURLConnection) new URL(location).openConnection();
        conn.setDoInput(true);
        conn.connect();
        status = conn.getResponseCode();
    }
    if (status == HttpURLConnection.HTTP_OK)
    {
        LOG.info("Got a OK reply");

        InputStream reader = conn.getInputStream();
        URLConnection.guessContentTypeFromStream(reader);
        StringBuilder builder = new StringBuilder();
        int a = 0;
        while ((a = reader.read()) != -1)
        {
            builder.append((char) a);
        }
        //System.out.println(builder.toString());
        sb.append(builder.toString());
    }
    else
        LOG.severe("Failed, got " + conn.getResponseMessage() + " for "
            + location);
    conn.disconnect();
    return sb.toString();
}

    
    public static void main(String[] args)
        throws Exception
    {
    	
		String uniprot_acc_id = "";
		String command = "";
		String format = "";
		String transcript_id = "";

        Getopt g = new Getopt("UniprotWebServiceClient", args, "c:i:f:t:");
        //
        int c;

        while ((c = g.getopt()) != -1)
        {
            switch(c)
            {
	        	case 'c':
	        		command = g.getOptarg();                	 
	        		break;
	
	        	case 'i':
	        		uniprot_acc_id = g.getOptarg();                	 
	        		break;

	        	case 'f':
	        		format = g.getOptarg();                	 
	        		break;

	        	case 't':
	        		transcript_id = g.getOptarg();                	 
	        		break;

				case '?':
					System.exit(1);
					break;
				
				default:
					System.exit(1);
					System.out.println("default");
            }
        }
		
        if (command.compareTo("")==0)
        {
        	System.err.println("Error!\nusage: java -jar UniprotWebServiceClient" +
        			" -i [ Command type (mapping | batch) ] " + 
        			" -i [ Uniprot AccessionID file (batch) ] " + 
        			" -t [ Transcript id (mapping) ] " + 
        			" -f [ Format ] (batch)");
        	System.exit(1);
        }
    	if(command.compareTo("batch")==0)
    	{
            if (uniprot_acc_id.compareTo("")==0 || format.compareTo("")==0)
            {
            	System.err.println("Error!\nusage: java -jar UniprotWebServiceClient" +
            			" -i [ Command type (mapping | batch) ] " + 
            			" -i [ Uniprot AccessionID file (batch) ] " + 
            			" -t [ Transcript id (mapping) ] " + 
            			" -f [ Format ] (batch)");
            	System.exit(1);
            }

            run("batch", new ParameterNameValue[] {
                	  new ParameterNameValue("query", uniprot_acc_id),
                  	  new ParameterNameValue("columns", "features"),
                  	  new ParameterNameValue("columns", "families"),
                  	  new ParameterNameValue("columns", "comments"),
                  	  new ParameterNameValue("format", format),
            });
    	}
    	else if(command.compareTo("mapping")==0)
    	{
            if (transcript_id.compareTo("")==0)
            {
            	System.err.println("Error!\nusage: java -jar UniprotWebServiceClient" +
            			" -i [ Command type (mapping | batch) ] " + 
            			" -i [ Uniprot AccessionID file (batch) ] " + 
            			" -t [ Transcript id (mapping) ] " + 
            			" -f [ Format ] (batch)");
            	System.exit(1);
            }
            run("mapping", new ParameterNameValue[] {
            	new ParameterNameValue("from", "ENSEMBL_TRS_ID"),
            	new ParameterNameValue("to", "ACC"),
                new ParameterNameValue("format", "tab"),
                new ParameterNameValue("query", transcript_id),
          });
    	}
    	else
    	{
    		System.err.println("Invalid command.");
    	}
    }

    public static class ParameterNameValue
    {
            private final String	name;
            private final String	value;

        public ParameterNameValue(String name, String value)
            throws UnsupportedEncodingException
        {
            this.name = URLEncoder.encode(name, "UTF-8");
            this.value = URLEncoder.encode(value, "UTF-8");
        }
    }
}
