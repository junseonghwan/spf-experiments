package pmcmc.proposals;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.MultivariateSummaryStatistics;

import briefj.BriefIO;
import briefj.run.Results;
import pmcmc.MCMCProposal;

public class MultivariateGaussianRandomWalk implements MCMCProposal<RealVectorParameters> {

	private double [] mu;
	private double [][] covariances;
	private int dim;

	public MultivariateGaussianRandomWalk(int dim, double [] sd) 
	{
		this(new double[dim], sd);
	}
	
	public MultivariateGaussianRandomWalk(double [] mu, double [] sd) 
	{
		this.mu = mu;
		this.dim = mu.length;
		this.covariances = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			this.covariances[i][i] = sd[i];
		}
	}
	
	@Override
	public RealVectorParameters initialize(Random random) {
		return propose(random, new RealVectorParameters(mu));
	}

	@Override
	public RealVectorParameters propose(Random random, RealVectorParameters curr) {
		double [] currVec = curr.getVector();

		MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(currVec, covariances);
		double [] vec = mvn.sample();
		return new RealVectorParameters(vec);
	}

	@Override
	public double logProposalDensity(RealVectorParameters curr, RealVectorParameters prev) {
		double [] currVec = curr.getVector();
		double [] prevVec = prev.getVector();
		MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(prevVec, covariances);
		return Math.log(mvn.density(currVec));
	}

	@Override
	public void adapt(List<RealVectorParameters> params) {
		// compute the sample covariance matrix
		MultivariateSummaryStatistics summ = new MultivariateSummaryStatistics(dim, true);
		for (RealVectorParameters param : params)
		{
			summ.addValue(param.getVector());
		}
		RealMatrix cov = summ.getCovariance().scalarMultiply(2.38*2.38/dim);
		System.out.println("Adapting covariance matrix: ");
		this.covariances = cov.getData();
		File resultsDir = Results.getResultFolder();
		PrintWriter writer = BriefIO.output(new File(resultsDir, "covMatrix.csv"));
		for (int i = 0; i < dim; i++)
		{
			StringBuilder sb = new StringBuilder();
			for (int j = 0; j < dim; j++) {
				System.out.print(covariances[i][j] + " ");
				sb.append(covariances[i][j] + " ");
			}
			System.out.println();
			writer.println(sb.toString());
		}
		writer.close();
	}
}
