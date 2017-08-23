package phylo.models;

import briefj.Indexer;
import phylo.EvolutionaryModel;

public class JukesCantor implements EvolutionaryModel
{
	private JukesCantorParam param = null;
	
	public JukesCantor(double mutationRate)
	{
		this.param = new JukesCantorParam(mutationRate);
	}

	@Override
	public double[][] getRateMatrix() 
	{
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		int numChars = dnaIndexer.size();
		double [][] Q = new double[numChars][numChars]; 
		for (int i = 0; i < numChars; i++)
		{
			Q[i][i] = -3*param.mutationRate/numChars;
			for (int j = i + 1; j < numChars; j++)
			{
				Q[i][j] = param.mutationRate/numChars;
				Q[j][i] = param.mutationRate/numChars;
			}
		}
		return Q;
	}

	@Override
	public double[] getStationaryDistribution() 
	{
		int b = DNAIndexer.indexer.size();
		double [] pi = new double[b];
		for (int i = 0; i < b; i++)
		{
			pi[i] = 0.25;
		}
		return pi;
	}	

	public static class JukesCantorParam
	{
		private double mutationRate;
		public JukesCantorParam(double mutationRate) {
			this.mutationRate = mutationRate;
		}
		public double getMutationRate() { return mutationRate; }
	}

}
