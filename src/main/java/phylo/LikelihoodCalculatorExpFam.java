package phylo;

import java.util.Collections;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import briefj.Indexer;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCExpFam.LearnedReversibleModel;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import phylo.RootedPhylogeny;

public class LikelihoodCalculatorExpFam implements LikelihoodCalculatorInterface
{

	private CTMCExpFam<String> model;
	private LearnedReversibleModel learnedModel;
	
	public LikelihoodCalculatorExpFam(CTMCExpFam<String> model, LearnedReversibleModel learnedModel)
	{
		this.model = model;
		this.learnedModel = learnedModel;
	}
	
	public DoubleMatrix getTransitionMatrix(double t)
	{
		DoubleMatrix Q = new DoubleMatrix(learnedModel.getRateMatrix());
		DoubleMatrix Qt = Q.mul(t);
		// exponentiate the rate matrix
		DoubleMatrix P = MatrixFunctions.expm(Qt);
		
		// check P is proper transition matrix
		if (!check(P))
			throw new RuntimeException("Not a propoer transition matrix");
		
		return P;
	}
	
	public static boolean check(DoubleMatrix P)
	{
		int R = P.rows;
		int C = P.columns;
		for (int r = 0; r < R; r++)
		{
			double sum = P.getRow(r).sum();
			if (Math.abs(sum - 1.0) > 1e-4)
			{
				System.out.println("r: " + sum);
				return false;
			}
		}
		return true;
	}
	
	public Indexer<String> getStateIndexer()
	{
		return model.stateIndexer;
	}
	
	/*
	 * compute the likelihood table
	 */
	@Override
	public double [][] computeLikelihoodTable(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2)
	{
		double [][] P1 = getTransitionMatrix(b1).toArray2();
		double [][] P2 = getTransitionMatrix(b2).toArray2();
		
		Indexer<String> stateIndexer = getStateIndexer();
		double [][] l1 = t1.getTaxon().getLikelihoodTable();
		double [][] l2 = t2.getTaxon().getLikelihoodTable();
		int s = l1.length;
		int b = stateIndexer.size();

		double [][] likelihoodTable = new double[s][b];

		// compute the dynamic programming table for the parent node
		for (int site = 0; site < s; site++)
		{
			for (int x = 0; x < b; x++)
			{
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (int y = 0; y < b; y++)
				{
					sum1 += (P1[x][y]*l1[site][y]);
					sum2 += (P2[x][y]*l2[site][y]);
				}
				likelihoodTable[site][x] = sum1*sum2;
			}
		}

		return likelihoodTable;
	}

	@Override
	public double computeLoglik(double [][] likelihoodTable)
	{
		double logLik = 0.0;
		int s = likelihoodTable.length;
		int b = this.getStateIndexer().size();
		for (int site = 0; site < s; site++)
		{
			double sum = 0.0;
			for (int x = 0; x < b; x++)
			{
				sum += likelihoodTable[site][x]*learnedModel.pi[x];
			}
			logLik += Math.log(sum);
		}
		return logLik;
	}
	
	public double computeLogLik(RootedPhylogeny phylogeny)
	{
		double logLik = 0.0;
		
		List<RootedPhylogeny> nodes = new ArrayList<>();
		nodes.add(phylogeny);
		
		while (!nodes.isEmpty())
		{
			RootedPhylogeny node = nodes.remove(0);
			if (!node.isLeaf()) {
				RootedPhylogeny left = node.getLeftChild();
				RootedPhylogeny right = node.getRightChild();
				nodes.add(left);
				nodes.add(right);
				logLik += computeBranchLogLik(node, left, node.getLeftBranchLength());
				logLik += computeBranchLogLik(node, right, node.getRightBranchLength());
			}
		}

		return logLik;
	}
	
	public double computeBranchLogLik(RootedPhylogeny parent, RootedPhylogeny child, double bl)
	{
		double [][] P = getTransitionMatrix(bl).toArray2();
		
		String parentSeq = parent.getTaxon().getSequence();
		String childSeq = child.getTaxon().getSequence();
		
		double logLik = 0.0;
		for (int s = 0; s < parentSeq.length(); s++)
		{
			int i = model.stateIndexer.o2i(parentSeq.charAt(s) + "");
			int j = model.stateIndexer.o2i(childSeq.charAt(s) + "");
			logLik += Math.log(P[i][j]);
		}
		return logLik;
	}

	public static void main(String [] args)
	{
		Random rand = new Random(1);
		
		Indexer<String> indexer = new Indexer<>();
		indexer.addToIndex("A", "C", "G", "T");
		CTMCExpFam<String> model = CTMCExpFam.createModelWithFullSupport(indexer, true);		
	    model.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<String>()));
	    model.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<String>()));

		int p = model.nFeatures();
		System.out.println("nFeatures: " + p);
	
	    // draw w ~ MVN(0, I)
		double [] mu = new double[p];
		double [][] I = new double[p][p];
		for (int i = 0; i < p; i++)
		{
			I[i][i] = 1.0;
		}

		MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new CustomRandomGenerator(rand.nextLong()), mu, I);
		double [] w = mvn.sample();

		LearnedReversibleModel learnedModel = model.reversibleModelWithParameters(w);
		LikelihoodCalculatorExpFam calc = new LikelihoodCalculatorExpFam(model, learnedModel);

		Taxon t1 = new Taxon("T1");
		Taxon t2 = new Taxon("T2");
		Taxon t3 = new Taxon("T3");
		
		double b1 = 1.2;
		double b2 = 0.3;
		double totalProb = 0.0;
		for (int x = 0; x < 4; x++)
		{
			String seq1 = model.stateIndexer.i2o(x);
			t1.setSequence(seq1);
			t1.initLikelihoodTable(model.stateIndexer);
			RootedPhylogeny phylo1 = new RootedPhylogeny(t1);
			for (int y = 0; y < 4; y++)
			{
				String seq2 = model.stateIndexer.i2o(y);
				t2.setSequence(seq2);
				t2.initLikelihoodTable(model.stateIndexer);
				RootedPhylogeny phylo2 = new RootedPhylogeny(t2);
				RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t0"), phylo1, phylo2, b1, b1, b1);
				
				double [][] likelihoodTable = calc.computeLikelihoodTable(phylo1, phylo2, b1, b1);
				double logLik = calc.computeLoglik(likelihoodTable);
				double lik = Math.exp(logLik);
				totalProb += lik;
			}
		}
		
		//DoubleMatrix P = calc.getTransitionMatrix(bl);
		//System.out.println(P);
		//System.out.println(new DoubleMatrix(likelihoodTable));
		
		System.out.println(totalProb); // should be 1.0
		
	}
}
