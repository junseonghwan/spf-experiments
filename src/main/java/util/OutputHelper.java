package util;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import briefj.BriefIO;

public class OutputHelper 
{
	public static void writeLines(File file, List<String> contents)
	{
		PrintWriter writer = BriefIO.output(file);
		for (String line : contents)
			writer.println(line);
		writer.close();
	}
	
	/**
	 * Helper function to write a vector of data to file as a column vector
	 * 
	 * @param filepath
	 * @param data
	 * @return
	 */
	public static void writeVector(File file, List<Double> data)
	{
		PrintWriter writer = BriefIO.output(file);
		for (double val : data)
		{
			writer.println(val);
		}
		writer.close();
	}

	public static void writeVector(String filepath, double [] data)
	{
		PrintWriter writer = BriefIO.output(new File(filepath));
		for (double val : data)
		{
			writer.println(val);
		}
		writer.close();
	}

	public static void writeTableAsCSV(File file, String [] header, double []... data)
	{
		if (header.length != data.length)
			throw new RuntimeException("The number of headers != number of columns in the data");
		
		PrintWriter writer = BriefIO.output(file);
		for (int i = 0; i < header.length; i++) {
			writer.print(header[i]);
			if (i < header.length - 1) writer.print(",");
			else writer.println();
		}
		for (int i = 0; i < data[0].length; i++)
		{
			for (int j = 0; j < data.length; j++) 
			{
				writer.print(data[j][i]);
				if (j < data.length - 1)
					writer.print(",");
				else
					writer.println();
			}
		}
		writer.close();
	}

}
