package phylo;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Multinomial;
import phylo.models.DNAIndexer;
import pmcmc.proposals.RealVectorParameters;

public class FelsensteinPruning implements LikelihoodCalculatorInterface 
{
	private EvolutionaryModel<RealVectorParameters> model;
	
	public FelsensteinPruning(EvolutionaryModel<RealVectorParameters> model)
	{
		this.model = model;
	}

	public EvolutionaryModel<RealVectorParameters> getModel() { return model; }
	
	@Override
	public Pair<Double, double [][]> computeLikelihoodTable(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2, boolean peek)
	{
		if (peek) {
			double logLik = computeLoglikInStream(t1, t2, b1, b2);
			return Pair.of(logLik, null);
		}
			
		double [][] transitionProbs1 = PhyloUtils.getTransitionMatrix(model, b1);
		double [][] transitionProbs2 = PhyloUtils.getTransitionMatrix(model, b2);

		double [][] table1 = t1.getTaxon().getLikelihoodTable();
		double [][] table2 = t2.getTaxon().getLikelihoodTable();

		//double [] pi = model.getStationaryDistribution();

		int numSites = table1.length;
		int numChars = table1[0].length;

		double [][] likelihoodTable = new double[numSites][numChars];

		double logRatio = 0.0;
		for (int s = 0; s < numSites; s++)
		{
			double currentNorm = 0.0;
			for (int x = 0; x < numChars; x++)
			{
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (int y = 0; y < numChars; y++)
				{
					sum1 += (transitionProbs1[x][y] * table1[s][y]);
					sum2 += (transitionProbs2[x][y] * table2[s][y]);
				}
				likelihoodTable[s][x] = sum1 * sum2;
				currentNorm += likelihoodTable[s][x];
			}
			Multinomial.normalize(likelihoodTable[s]);
			logRatio += Math.log(currentNorm);
		}
		double newLogLikelihood = logRatio + t1.logLikelihood() + t2.logLikelihood(); // loglikelihood for the merged tree

		return Pair.of(newLogLikelihood, likelihoodTable);
	}
	
	@Override
	public double computeLoglik(double [][] likelihoodTable)
	{
		double logLik = 0.0;
		int s = likelihoodTable.length;
		int b = DNAIndexer.indexer.size();
		double [] pi = model.getStationaryDistribution();
		for (int site = 0; site < s; site++)
		{
			double sum = 0.0;
			for (int x = 0; x < b; x++)
			{
				sum += likelihoodTable[site][x]*pi[x];
			}
			logLik += Math.log(sum);
		}
		return logLik;
	}
	
	@Override
	public double computeLoglikInStream(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2) 
	{
		double [][] transitionProbs1 = PhyloUtils.getTransitionMatrix(model, b1);
		double [][] transitionProbs2 = PhyloUtils.getTransitionMatrix(model, b2);

		double [][] table1 = t1.getTaxon().getLikelihoodTable();
		double [][] table2 = t2.getTaxon().getLikelihoodTable();

		double [] pi = model.getStationaryDistribution();

		int numSites = table1.length;
		int numChars = table1[0].length;

		double logRatio = 0.0;

		for (int s = 0; s < numSites; s++)
		{
			double C = 0.0;
			for (int x = 0; x < numChars; x++)
			{
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (int y = 0; y < numChars; y++)
				{
					sum1 += (transitionProbs1[x][y] * table1[s][y]);
					sum2 += (transitionProbs2[x][y] * table2[s][y]);
				}
				C += pi[x] * sum1 * sum2;
			}
			logRatio += Math.log(C);
		}
		double newLogLikelihood = logRatio + t1.logLikelihood() + t2.logLikelihood(); // loglikelihood for the merged tree
		//System.out.println("logRatio: " + logRatio + ", newLogLik: " + newLogLikelihood);
		return newLogLikelihood;
	}
}
