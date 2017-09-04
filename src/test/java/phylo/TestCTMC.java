package phylo;

import org.junit.Assert;
import org.junit.Test;

import bayonet.math.NumericalUtils;
import phylo.models.Felsenstein81;
import phylo.models.JukesCantorModel;

public class TestCTMC 
{
	@Test
	public void testMatrixExponentiation()
	{
		double t = 1.2;
		double mu = 1.2;
		testJC(t, mu);
		testJC(0.0, 1.0);
		testJC(0.0, 0.0);
		testJC(1.0, 0.0);
		testJC(1.0, 1.0);
		
		double [] pi = new double[]{0.1, 0.3, 0.2, 0.4};
		testF81(pi, 32);
	}
	
	public void testF81(double [] pi, double t)
	{
		EvolutionaryModel f81 = new Felsenstein81(pi);
		double [][] P = PhyloUtils.getTransitionMatrix(f81, t);
		double sumSq = 0.0;
		for (int i = 0; i < 4; i++)
			sumSq += Math.pow(pi[i], 2.0);
		double beta = 1.0/(1 - sumSq);
		double expBetaT = Math.exp(-beta * t);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				double expected = 0.0;
				if (i == j) {
					expected = expBetaT + pi[j] * (1 - expBetaT);
				} else {
					expected = pi[j] * (1 - expBetaT);
				}
				System.out.println("P[" + i + "][" + j + "]: " + expected + ", " + P[i][j]);
				Assert.assertTrue(NumericalUtils.isClose(expected, P[i][j], 1e-6));
			}
		}		
	}
	
	public void testJC(double t, double mu)
	{
		EvolutionaryModel jcModel = new JukesCantorModel(mu);
		double [][] P = PhyloUtils.getTransitionMatrix(jcModel, t);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				double expected = 0.0;
				if (i == j) {
					expected = 0.25 + 0.75 * Math.exp(-t * mu);
				} else {
					expected = 0.25 - 0.25 * Math.exp(-t * mu);
				}
				System.out.println("P[" + i + "][" + j + "]: " + expected + ", " + P[i][j]);
				Assert.assertTrue(NumericalUtils.isClose(expected, P[i][j], 1e-6));
			}
		}		
	}

}
