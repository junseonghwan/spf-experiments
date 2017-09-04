package phylo.models;

import org.jblas.DoubleMatrix;

import briefj.Indexer;

import phylo.EvolutionaryModel;
import pmcmc.proposals.RealVectorParameters;

public class JukesCantorModel implements EvolutionaryModel<RealVectorParameters>
{
	private RealVectorParameters param = null;
	private RealVectorParameters old = null;
	private DoubleMatrix Q;
	
	public JukesCantorModel(double mutationRate)
	{
		this.param = new RealVectorParameters(new double[]{mutationRate});
	}

	public JukesCantorModel(RealVectorParameters param)
	{
		if (param.getDim() != 1)
			throw new RuntimeException("The Jukes Cantor model parameter only has 1 parameter.");
		this.param = param;
		constructRateMatrix(param.getVector()[0]);
	}
	
	public void constructRateMatrix(double mutationRate)
	{
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		int numChars = dnaIndexer.size();
		double [][] Q = new double[numChars][numChars]; 
		for (int i = 0; i < numChars; i++)
		{
			Q[i][i] = -3*mutationRate/numChars;
			for (int j = i + 1; j < numChars; j++)
			{
				Q[i][j] = mutationRate/numChars;
				Q[j][i] = mutationRate/numChars;
			}
		}
		this.Q = new DoubleMatrix(Q);
	}

	@Override
	public DoubleMatrix getRateMatrix() 
	{
		if (Q == null)
			constructRateMatrix(this.param.getVector()[0]);
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

	@Override
	public RealVectorParameters getModelParameters() {
		return param;
	}

	@Override
	public void updateModelParameters(RealVectorParameters p) {
		this.old = this.param;
		this.param = p;
		constructRateMatrix(this.param.getVector()[0]);
	}

	@Override
	public void revert() {
		this.param = this.old;
		this.old = null;
		constructRateMatrix(this.param.getVector()[0]);
	}	

	/*
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
			PhyloOptions.calc = new FelsensteinPruningAlgorithm(new JukesCantorModel(this.mutationRate.getValue()));
		}
		@Override
		public String asCommaSeparatedLine() {
			return mutationRate.getValue() + "";
		}
	}

	@Override
	public RealVectorParameters getModelParameters() {
		// TODO Auto-generated method stub
		return null;
	}
	*/

}
