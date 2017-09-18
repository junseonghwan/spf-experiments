package dynamic.models;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import pmcmc.Model;
import pmcmc.proposals.RealVectorParameters;
import bayonet.distributions.Normal;
import bayonet.distributions.Poisson;
import briefj.BriefIO;
import briefj.Indexer;

public class RickerModel implements Model<RealVectorParameters>
{
	private RealVectorParameters params;
	private RealVectorParameters old;
	private Indexer<String> paramIndexer;

	public static String [] paramNames1 = new String[]{"N0", "phi", "var", "r"};
	public static String [] paramNames2 = new String[]{"phi", "var", "r"};
	
	public RickerModel(double N0, double phi, double var, double r)
	{
		params = new RealVectorParameters(new double[]{N0, phi, var, r});
		this.paramIndexer = new Indexer<>();
		for (String paramName : paramNames1)
			this.paramIndexer.addToIndex(paramName);
	}

	public RickerModel(double phi, double var, double r)
	{
		params = new RealVectorParameters(new double[]{phi, var, r});
		this.paramIndexer = new Indexer<>();
		for (String paramName : paramNames2)
			this.paramIndexer.addToIndex(paramName);
	}

	public static List<Pair<List<Double>, List<Integer>>> generate(Random random, int numSimulations, int T, double N0, double phi, double var, double r)
	{
		// generate the latent variable and the data (all at once)
		List<Pair<List<Double>, List<Integer>>> data = new ArrayList<>();
		for (int i = 1; i <= numSimulations; i++)
		{
			Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(new Random(random.nextLong()), T, N0, phi, var, r);
			data.add(ret);
		}
		
		return data;
	}

	public static List<Pair<List<Double>, List<Integer>>> generate(Random random, int numSimulations, int T, int maxN0, double phi, double var, double r)
	{
		// generate the latent variable and the data (all at once)
		List<Pair<List<Double>, List<Integer>>> data = new ArrayList<>();
		for (int i = 1; i <= numSimulations; i++)
		{
			Pair<List<Double>, List<Integer>> ret = RickerModel.simulate(new Random(random.nextLong()), T, random.nextDouble()*maxN0, phi, var, r);
			data.add(ret);
		}
		
		return data;
	}

	public static Pair<List<Double>, List<Integer>> readFromFile(String path)
	{
		List<Double> Ns = new ArrayList<>();
		List<Integer> ys = new ArrayList<>();
		for (String line : BriefIO.readLines(new File(path)))
		{
			String [] row = line.split(",");
			try {
				int y = Integer.parseInt(row[0]);
				double N = Double.parseDouble(row[1]);
				Ns.add(N);
				ys.add(y);
			} catch (Exception ex) {
				
			}
		}
		return Pair.of(Ns, ys);
	}
	
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

	@Override
	public RealVectorParameters getModelParameters() {
		return params;
	}

	@Override
	public void updateModelParameters(RealVectorParameters p) {
		this.old = this.params;
		this.params = p;
	}

	@Override
	public void revert() {
		this.params = this.old;
		this.old = null;
	}
	
	public double getParamValue(String paramName)
	{
		if (!paramIndexer.containsObject(paramName))
			throw new RuntimeException("Parameter name " + paramName + " does not exist for Ricker model.");
		return this.params.getVector()[paramIndexer.o2i(paramName)];
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
