package phylo.models;

import org.jblas.DoubleMatrix;

import blang.variables.RealVariable;
import briefj.Indexer;
import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PhyloOptions;
import pmcmc.ModelParameters;

public class JukesCantor implements EvolutionaryModel
{
	private JukesCantorParam param = null;
	private DoubleMatrix Q;
	
	public JukesCantor(double mutationRate)
	{
		this.param = new JukesCantorParam(mutationRate);
		constructRateMatrix();
	}
	
	public void constructRateMatrix()
	{
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		int numChars = dnaIndexer.size();
		double [][] Q = new double[numChars][numChars]; 
		for (int i = 0; i < numChars; i++)
		{
			Q[i][i] = -3*param.getMutationRate()/numChars;
			for (int j = i + 1; j < numChars; j++)
			{
				Q[i][j] = param.getMutationRate()/numChars;
				Q[j][i] = param.getMutationRate()/numChars;
			}
		}
		this.Q = new DoubleMatrix(Q);
	}

	@Override
	public DoubleMatrix getRateMatrix() 
	{
		if (Q == null)
			constructRateMatrix();
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

	public static class JukesCantorParam implements ModelParameters
	{
		private RealVariable mutationRate;
		private RealVariable oldValue;
		public JukesCantorParam(double mutationRate) {
			this.mutationRate = new RealVariable(mutationRate);
		}
		public double getMutationRate() { return mutationRate.getValue(); }

		@Override
		public void update(ModelParameters p) {
			this.oldValue = this.mutationRate;
			this.mutationRate = ((JukesCantorParam)p).mutationRate;
		}

		@Override
		public void revert() {
			this.mutationRate = this.oldValue;
			this.oldValue = null;
			PhyloOptions.calc = new FelsensteinPruningAlgorithm(new JukesCantor(this.mutationRate.getValue()));
		}
		@Override
		public String asCommaSeparatedLine() {
			return mutationRate.getValue() + "";
		}
	}

}
