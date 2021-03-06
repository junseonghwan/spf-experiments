package phylo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import phylo.models.Coalescent;
import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCExpFam.LearnedReversibleModel;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import conifer.io.Indexers;

public class GenerateSequenceDataFromExpFam {
	
	// generate the data along the phylogeny
	@SuppressWarnings("rawtypes")
	public static void generateSequencesCTMC(Random rand, CTMCExpFam<String> model, LearnedReversibleModel learnedModel, RootedPhylogeny phylogeny, int numSites)
	{
		// get the transition probabilities and the rate matrix from the model
		//double [][] transProbs = getTransitionProbs(learnedModel.getRateMatrix());
		//double [][] Q = learnedModel.getRateMatrix();

		// generate the data at the root from the stationary distribution
		String sequence = generateSequenceStationary(rand, model, numSites, learnedModel.pi);
		phylogeny.getTaxon().setSequence(sequence);

		// traverse the phylogeny in preorder (do it by recursion)
		generateSequence(rand, model, phylogeny, learnedModel); // right subtree
	}
	
	// helper function for recursively generating the data
	// generate the data for the children
	//private static void generateSequence(Random rand, CTMCExpFam<String> model, RootedPhylogeny tree, double [][] probs, double [][] Q)
	private static void generateSequence(Random rand, CTMCExpFam<String> model, RootedPhylogeny tree, LearnedReversibleModel learnedReversibleModel)
	{
		String sequence = tree.getTaxon().getSequence();
		RootedPhylogeny t1 = tree.getLeftChild();
		RootedPhylogeny t2 = tree.getRightChild();
		double b1 = tree.getLeftBranchLength();
		double b2 = tree.getRightBranchLength();

		String leftSequence = simulateSequence(rand, model, sequence, b1, learnedReversibleModel);
		String rightSequence = simulateSequence(rand, model, sequence, b2, learnedReversibleModel);
		
		t1.getTaxon().setSequence(leftSequence);
		t2.getTaxon().setSequence(rightSequence);
		
		if (!t1.isLeaf())
		{
			generateSequence(rand, model, t1, learnedReversibleModel);
		}
		if (!t2.isLeaf())
		{
			generateSequence(rand, model, t2, learnedReversibleModel);
		}
	}
	
	private static String simulateSequence(Random rand, CTMCExpFam<String> model, String sequence, double branchLength, LearnedReversibleModel learnedReversibleModel)
	{
		double [][] probs = getTransitionProbs(learnedReversibleModel, branchLength);
		StringBuilder newSequence = new StringBuilder();
		for (int i = 0; i < sequence.length(); i++)
		{
			String ch = sequence.charAt(i) + "";
			int state = model.stateIndexer.o2i(ch);
			int newState = Multinomial.sampleMultinomial(rand, probs[state]);
			String newCh = model.stateIndexer.i2o(newState);
			newSequence.append(newCh);
		}
		return newSequence.toString();
	}
	
	/*
	private static String simulateSequence(Random rand, CTMCExpFam<String> model, String sequence, double [][] probs, double [][] Q, double branchLength)
	{
		String newSequence = "";
		for (int i = 0; i < sequence.length(); i++)
		{
			double t = 0.0;
			String ch = sequence.charAt(i) + "";
			int state = model.stateIndexer.o2i(ch);
			while (true)
			{
	    		double holdingTime = Exponential.generate(rand, -1.0/Q[state][state]);
	    		t += holdingTime;

	    		if (t > branchLength)
	    			break;

	    		state = Multinomial.sampleMultinomial(rand, probs[state]);
	  			ch = model.stateIndexer.i2o(state);
			}
			newSequence += ch;
		}
		return newSequence;
	}
	*/
	
	public static String generateSequenceStationary(Random rand, CTMCExpFam<String> model, int numSites, double [] pi)
	{
		String seq = "";
		for (int s = 0; s < numSites; s++)
		{
	  		int state = Multinomial.sampleMultinomial(rand, pi);
	  		String ch = model.stateIndexer.i2o(state);
	  		seq += ch;
		}

		return seq;
	}

	/*
	public static double [][] getTransitionProbs(double [][] rateMatrix)
	{
		int numStates = rateMatrix[0].length;
		double [][] transProbs = new double[numStates][numStates];
		for (int i = 0; i < numStates; i++)
		{
			double norm = -rateMatrix[i][i];
  		for (int j = 0; j < numStates; j++)
  		{
  			if (i != j) 
  				transProbs[i][j] = rateMatrix[i][j]/norm;
  		}
		}
		
		return transProbs;
	}
	*/
	
	public static double[][] getTransitionProbs(LearnedReversibleModel learnedModel, double t)
	{
		DoubleMatrix Q = new DoubleMatrix(learnedModel.getRateMatrix());
		DoubleMatrix Qt = Q.mul(t);
		// exponentiate the rate matrix
		DoubleMatrix P = MatrixFunctions.expm(Qt);
		
		// check P is proper transition matrix
		if (!LikelihoodCalculatorExpFam.check(P))
			throw new RuntimeException("Not a propoer transition matrix");
		
		return P.toArray2();
	}


	public static void main(String [] args)
	{
		// test the data generation code
		Random rand = new Random(102);
		List<Taxon> taxa = new ArrayList<Taxon>();
		taxa.add(new Taxon("T0"));
		taxa.add(new Taxon("T1"));
		taxa.add(new Taxon("T2"));
		taxa.add(new Taxon("T3"));

		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(rand, taxa);

		int numSites = 10;

	    CTMCExpFam<String> model = CTMCExpFam.createModelWithFullSupport(Indexers.dnaIndexer(), true);
	    model.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<String>()));
	    model.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<String>()));

		int p = model.nFeatures();
		System.out.println("nFeatures: " + p);

		// generate w ~ mvn(0, I)
		double [] mu = new double[p];
		double [][] I = new double[p][p];
		for (int i = 0; i < p; i++)
		{
			I[i][i] = 1.0;
		}

		MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new CustomRandomGenerator(rand.nextLong()), mu, I);
		double [] w = mvn.sample();

		@SuppressWarnings("rawtypes")
		LearnedReversibleModel learnedModel = model.reversibleModelWithParameters(w);

		GenerateSequenceDataFromExpFam.generateSequencesCTMC(rand, model, learnedModel, phylogeny, numSites);

		System.out.println(phylogeny.getTreeString());
		System.out.println(phylogeny.getDataString());

	}

}
