package models;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Normal;

public class Kitagawa 
{
	public static Pair<double [], double []> simulate(Random random, double var_v, double var_w, int numPoints)
	{
		double [] x = new double[numPoints];
		double [] y = new double[numPoints];
		x[0] = Normal.generate(random, 0.0, 5.0);
		y[0] = Math.pow(x[0], 2) / 20.0 + Normal.generate(random, 0.0, var_w);
		for (int i = 1; i < numPoints; i++)
		{
			x[i] = x[i-1]/2 + 25 * x[i-1] / (1 + Math.pow(x[i-1], 2.0)) + 8 * Math.cos(1.2*i) + Normal.generate(random, 0.0, var_v);
			y[i] = Math.pow(x[i], 2.0) / 20.0 + Normal.generate(random, 0.0, var_w);
		}

		return Pair.of(x, y);
	}
	
	public static void main(String [] args)
	{
		Pair<double [], double []> ret = Kitagawa.simulate(new Random(1), 10, 10, 100);
		double [] x = ret.getLeft();
		double [] y = ret.getRight();
		for (int i = 0; i < x.length; i++)
		{
			System.out.println(x[i] + ", " + y[i]);
		}
	}

}
