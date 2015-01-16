package javagene.util;

import java.io.*;

/**
* Write text messages to the console.
*
* This class is a placeholder; rewrite it to
* fit your needs.
*
*/
public class Log
{
	
	private Log(){};
	
	/**
	 * Write a blank line.
	 */
	public static void log()
	{
		log( "" );
	}
	
	/**
	 * Write a line to the console.
	 * @param text What to write.
	 */
	public static void log( String text )
	{
		System.out.println( text );
	}
	


}