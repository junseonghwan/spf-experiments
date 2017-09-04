package dynamic.models;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import pmcmc.Model;
import pmcmc.proposals.RealVectorParameters;
import bayonet.distributions.Normal;

public class KitagawaModel implements Model<RealVectorParameters>
{
	private RealVectorParameters params;
	private RealVectorParameters old;
	
	public KitagawaModel(double var_v, double var_w)
	{
		this.params = new RealVectorParameters(new double[]{var_v, var_w});
	}
	
	public double getVarV()
	{
		return this.params.getVector()[0];
	}

	public double getVarW()
	{
		return this.params.getVector()[1];
	}

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

	@Override
	public RealVectorParameters getModelParameters() {
		return this.params;
	}

	@Override
	public void updateModelParameters(RealVectorParameters p) {
		this.old = params;
		this.params = p;
	}

	@Override
	public void revert() {
		this.params = this.old;
		this.old = null;
	}

	public static void main(String [] args)
	{
		Pair<double [], double []> ret = KitagawaModel.simulate(new Random(1), 10, 10, 100);
		double [] x = ret.getLeft();
		double [] y = ret.getRight();
		for (int i = 0; i < x.length; i++)
		{
			System.out.println(x[i] + ", " + y[i]);
		}
	}
}
