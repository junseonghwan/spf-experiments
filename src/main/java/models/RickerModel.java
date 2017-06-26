package models;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Normal;
import bayonet.distributions.Poisson;

public class RickerModel 
{
	/**
	 * Return a pair of lists of length T + 1 (including N0 and the observation generated for time 0)
	 * 
	 * @param random
	 * @param T
	 * @param N0
	 * @param phi
	 * @param var
	 * @param r
	 * @return
	 */
	public static Pair<List<Double>, List<Integer>> simulate(Random random, int T, double N0, double phi, double var, double r)
	{
		List<Double> latent = new ArrayList<>();
		List<Integer> obs = new ArrayList<>();

		latent.add(N0);
		obs.add(Poisson.generate(random, phi*N0));
		
		for (int t = 1; t <= T; t++)
		{
			double et = Normal.generate(random, 0.0, var);
			double Nt = r * latent.get(t-1) * Math.exp(-latent.get(t-1) + et);
			int Yt = Poisson.generate(random, phi * Nt);

			obs.add(Yt);
			latent.add(Nt);
		}

		return Pair.of(latent, obs);
	}

	public static void main(String [] args)
	{
		Random random = new Random(123);
		int T = 10;
		Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(random, T, 7, 10.0, 0.3, 44.7);
		for (int i = 0; i <= T; i++)
		{
			System.out.println(ret.getLeft().get(i) + ", " + ret.getRight().get(i));
		}
	}
}
