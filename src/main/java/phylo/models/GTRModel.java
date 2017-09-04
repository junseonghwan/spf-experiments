package phylo.models;

import org.jblas.DoubleMatrix;

import briefj.Indexer;
import phylo.EvolutionaryModel;
import pmcmc.proposals.RealVectorParameters;

public class GTRModel implements EvolutionaryModel<RealVectorParameters>
{
	private RealVectorParameters gtrParams;
	private RealVectorParameters oldParams;
	private double [] pi;
	private DoubleMatrix Q;
	private double [][] params;
	public GTRModel(double [] pi, double [] params)
	{
		this.pi = pi;
		this.gtrParams = new RealVectorParameters(params);
		constructRateMatrix();
	}
	
	private void constructRateMatrix()
	{
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		int numChars = dnaIndexer.size();
		double [][] Q = new double[numChars][numChars];
		double [] vec = gtrParams.getVector();
		params = new double[numChars][numChars];
		params[0][1] = vec[0];
		params[0][2] = vec[1];
		params[0][3] = vec[2];
		params[1][2] = vec[3];
		params[1][3] = vec[4];
		params[2][3] = vec[5];

		for (int i = 0; i < numChars; i++)
		{
			for (int j = i + 1; j < numChars; j++)
			{
				Q[i][j] = params[i][j] * pi[j];
				Q[j][i] = Q[i][j];
			}
		}
		for (int i = 0; i < numChars; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < numChars; j++)
			{
				if (i != j)
					sum += Q[i][j]; 
			}
			Q[i][i] = -sum;
		}
		this.Q = new DoubleMatrix(Q);
	}

	@Override
	public DoubleMatrix getRateMatrix() {
		return Q;
	}

	@Override
	public double[] getStationaryDistribution() {
		return this.pi;
	}
	
	public static class GTRModelParams
	{
		double alpha, beta, gamma, delta, epsilon, eta;
		public GTRModelParams(double a, double b, double g, double d, double e, double f)
		{
			this.alpha = a;
			this.beta = b;
			this.gamma = g;
			this.delta = d;
			this.epsilon = e;
			this.eta = f;
		}
		
		@Override
		public String toString()
		{
			return alpha + ", " + beta + ", " + gamma + ", " + delta + ", " + epsilon + ", " + eta;
		}
	}

	@Override
	public RealVectorParameters getModelParameters() {
		return gtrParams;
	}

	@Override
	public void updateModelParameters(RealVectorParameters p) {
		this.oldParams = this.gtrParams;
		this.gtrParams = p;
	}

	@Override
	public void revert() {
		this.gtrParams = this.oldParams;
		this.oldParams = null;
	}

}
