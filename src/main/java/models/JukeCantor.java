package models;

import java.util.Random;

import bayonet.distributions.Multinomial;
import briefj.Indexer;

public class JukeCantor 
{
	private static Indexer<Character> indexer = DNAIndexer.indexer;
	
	public static double [][] getTransitionProbabilities(double t, double pi)
	{
		double [][] P = new double[4][4];
		
		for (int i = 0; i < 4; i++)
		{
			P[i][i] = 0.25 + 0.75 * Math.exp(-t * pi);
			for (int j = i + 1; j < 4; j++)
			{
				P[i][j] = 0.25 - 0.25 * Math.exp(-t * pi);
				P[j][i] = 0.25 - 0.25 * Math.exp(-t * pi);
			}
		}
		return P;
	}

	public static String simulate(Random random, final String parentSequence, double t, double pi)
	{
		StringBuilder childSequence = new StringBuilder(parentSequence);
		double [][] P = getTransitionProbabilities(t, pi);
		for (int i = 0; i < parentSequence.length(); i++)
		{
			double [] p = P[indexer.o2i(parentSequence.charAt(i))];
			// sample
			int j = Multinomial.sampleMultinomial(random, p);
			childSequence.setCharAt(i, indexer.i2o(j));
		}
		return childSequence.toString();
	}

	public static void main(String [] args)
	{
		String parent = "AAAAACCCGGGCCCTACGGT";
		Random random = new Random(123);
		System.out.println(parent);
		double mu = 0.3;
		double t = 0.1;
		System.out.println(simulate(random, parent, t, mu));
	}
}
