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
	public static <T> void writeVector(File file, List<T> data)
	{
		PrintWriter writer = BriefIO.output(file);
		for (T val : data)
		{
			writer.println(val.toString());
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

	public static void writeVector(File file, double [] data)
	{
		PrintWriter writer = BriefIO.output(file);
		for (double val : data)
		{
			writer.println(val);
		}
		writer.close();
	}

	public static <T> void writeTableAsCSV(File file, int numColumns, List<List<T>> data)
	{
		PrintWriter writer = BriefIO.output(file);
		for (List<T> d : data)
		{
			if (d.size() != numColumns)
				throw new RuntimeException("number of columns != actual data size");
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < numColumns; i++)
			{
				if (i == numColumns - 1)
					sb.append(d.get(i).toString());
				else
					sb.append(d.get(i).toString() + ", ");
			}
			writer.println(sb.toString());
		}
		writer.close();
	}
	
	public static void writeDoubleArray2DAsCSV(File file, String [] header, double [][] data)
	{
		PrintWriter writer = BriefIO.output(file);
		for (int i = 0; i < header.length; i++) {
			writer.print(header[i]);
			if (i < header.length - 1) writer.print(",");
			else writer.println();
		}
		for (int i = 0; i < data.length; i++)
		{
			for (int j = 0; j < data[i].length; j++) 
			{
				writer.print(data[i][j]);
				if (j < data[i].length - 1)
					writer.print(",");
				else
					writer.println();
			}
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
