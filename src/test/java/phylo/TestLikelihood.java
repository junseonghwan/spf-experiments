package phylo;

import org.junit.Test;

import bayonet.math.NumericalUtils;

import java.util.Random;

import org.junit.Assert;
import phylo.models.DNAIndexer;
import phylo.models.JukesCantorModel;

/**
 * Test the likelihood computation.
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class TestLikelihood 
{
	public Random random = new Random(123);
	
	@Test
	public void test1()
	{
		String s1 = "ACGT";
		String s2 = "CCAA";
		double b1 = 1.5;
		double b2 = 2.2;
		double mu = 0.1;
		performTest(s1, s2, b1, b2, mu);
		
		int S = random.nextInt(1000);
		StringBuilder sb1 = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		for (int s = 0; s < S; s++)
		{
			sb1.append(DNAIndexer.indexer.i2o(random.nextInt(4)));
			sb2.append(DNAIndexer.indexer.i2o(random.nextInt(4)));
		}
		performTest(sb1.toString(), sb2.toString(), random.nextDouble()*10, random.nextDouble()*10, random.nextDouble()*10);
	}
	
	/**
	 * A simple test of likelihood computation on a cherry computed implemented manually versus the implementation in FelsensteinPruningAlgorithm 
	 */
	public void performTest(String s1, String s2, double b1, double b2, double mu)
	{
		Taxon taxon1 = new Taxon("t1", s1);
		Taxon taxon2 = new Taxon("t2", s2);
		int S = s1.length();

		EvolutionaryModel jcModel = new JukesCantorModel(mu);
		double [] pi = jcModel.getStationaryDistribution();
		double [][] P1 = PhyloUtils.getTransitionMatrix(jcModel, b1);
		double [][] P2 = PhyloUtils.getTransitionMatrix(jcModel, b2);
		double [][] likTableTruth = new double[s1.length()][4];
		for (int s = 0; s < s1.length(); s++)
		{
			for (String str : DNAIndexer.indexer.objectsList())
			{
				int i = DNAIndexer.indexer.o2i(str);
				String str1 = s1.charAt(s) + "";
				String str2 = s2.charAt(s) + "";
				int j1 = DNAIndexer.indexer.o2i(str1);
				int j2 = DNAIndexer.indexer.o2i(str2);
				likTableTruth[s][i] = P1[i][j1] * P2[i][j2];
			}
		}
		double trueMarginalLogLikelihood = 0.0;
		for (int s = 0; s < S; s++)
		{
			double sum = 0.0;
			for (String str1 : DNAIndexer.indexer.objectsList())
			{
				int i = DNAIndexer.indexer.o2i(str1);
				sum += pi[i] * likTableTruth[s][i];
			}
			trueMarginalLogLikelihood += Math.log(sum);
		}

		// construct likelihood table for the parent node and compare against the hand computation
		FelsensteinPruningAlgorithm peeling = new FelsensteinPruningAlgorithm(jcModel);
		RootedPhylogeny t1 = new RootedPhylogeny(taxon1);
		RootedPhylogeny t2 = new RootedPhylogeny(taxon2);
		double [][] likTable = peeling.computeLikelihoodTable(t1, t2, b1, b2, false).getRight();

		for (int s = 0; s < s1.length(); s++)
		{
			for (String str : DNAIndexer.indexer.objectsList())
			{
				int b = DNAIndexer.indexer.o2i(str);
				System.out.println(likTable[s][b] + ", " + likTableTruth[s][b]);
				Assert.assertTrue(NumericalUtils.isClose(likTable[s][b] - likTableTruth[s][b], 0.0, 1e-6));
			}
		}
		
		double marginalLogLikelihood = peeling.computeLoglik(likTable);
		System.out.println(trueMarginalLogLikelihood + ", " + marginalLogLikelihood);
		Assert.assertTrue(NumericalUtils.isClose(trueMarginalLogLikelihood - marginalLogLikelihood, 0.0, 1e-6));
	}

	
}
