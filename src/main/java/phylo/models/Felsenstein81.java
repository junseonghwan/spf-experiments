package phylo.models;

import bayonet.math.NumericalUtils;
import briefj.Indexer;
import phylo.EvolutionaryModel;

public class Felsenstein81 implements EvolutionaryModel
{
	private double [] pi;
	private double [][] Q;
	public Felsenstein81(double [] pi)
	{
		this.pi = pi;

		double sum = 0.0;
		for (double w : pi)
			sum += w;
		NumericalUtils.checkIsClose(sum, 1.0, 1e-6);

		constructRateMatrix(pi);
	}
	
	private void constructRateMatrix(double [] pi)
	{
		double sum = 0.0;
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		int numChars = dnaIndexer.size();
		Q = new double[numChars][numChars]; 
		for (int i = 0; i < numChars; i++)
		{
			sum = 0.0;
			for (int j = 0; j < numChars; j++)
			{
				if (i != j) {
					Q[i][j] = pi[j];
					sum += pi[j];
				}
			}
			Q[i][i] = -sum;
		}		
	}

	@Override
	public double[][] getRateMatrix() 
	{
		if (Q == null)
			constructRateMatrix(this.pi);
		return Q;
	}

	@Override
	public double[] getStationaryDistribution() 
	{
		return pi;
	}

}
