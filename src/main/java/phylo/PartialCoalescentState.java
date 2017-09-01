package phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;

public class PartialCoalescentState {

	private List<RootedPhylogeny> trees = new ArrayList<RootedPhylogeny>();
	private int id = 0;
	private double height = 0.0;
	
	private static PartialCoalescentState initial = null;

	private PartialCoalescentState()
	{
		// prevent initialization outside of this class
	}
	
	public static PartialCoalescentState copy(PartialCoalescentState src)
	{
		PartialCoalescentState newState = new PartialCoalescentState();
		for (int i = 0; i < src.trees.size(); i++)
		{
			newState.trees.add(src.trees.get(i));
		}
		newState.id = src.id;
		newState.height = src.height;
		return newState;
	}
	
	// get initial state -- can only instantiate through this method
	public static PartialCoalescentState getInitial(List<Taxon> taxa)
	{
		if (initial == null)
		{
  		initial = new PartialCoalescentState();
  		for (Taxon taxon : taxa)
  		{
  			RootedPhylogeny phylo = new RootedPhylogeny(taxon);
  			//phylo.setLogLikelihood(PhyloOptions.calc.computeLoglik(taxon.getLikelihoodTable()));
  			initial.trees.add(phylo);
  		}
		}
		return initial;
	}
	
	// Leaves the current state un-touched, including the trees (RootedPhylogeny) and everything else related to it
	public Pair<Double, PartialCoalescentState> coalesce(Random rand, boolean isPeek)
	{
		PartialCoalescentState newState = new PartialCoalescentState();
		for (int i = 0; i < trees.size(); i++)
		{
			newState.trees.add(trees.get(i));
		}
		newState.id = this.id;

		// sample from exponential distribution
		double branchLength = Exponential.generate(rand, nChoose2(newState.trees.size()));
		// select two trees to merge at random
		RootedPhylogeny t1 = newState.trees.remove(rand.nextInt(newState.trees.size()));
		RootedPhylogeny t2 = newState.trees.remove(rand.nextInt(newState.trees.size()));
		
		newState.height = this.height + branchLength;
		
		// determine branch length for each of the subtrees to the parent
		double b1 = newState.height - t1.getHeight();
		double b2 = newState.height - t2.getHeight();
		
		// compute the likelihood
		double logw = 0.0;
		if (isPeek) {
			logw = PhyloOptions.calc.computeLoglikInStream(t1, t2, b1, b2);
		} else {
			RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + newState.id++), t1, t2, b1, b2, newState.height, true);
			newState.trees.add(parent);

			double [][] likelihoodTable = PhyloOptions.calc.computeLikelihoodTable(t1, t2, b1, b2);
			parent.getTaxon().setLikelihoodTable(likelihoodTable);
			logw = PhyloOptions.calc.computeLoglik(likelihoodTable);
			parent.setLogLikelihood(logw + t1.logLikelihood() + t2.logLikelihood()); // this is the likelihood for the tree
		}

		//System.out.println(newState.trees.size() + ", " + t1.getTaxon().toString() + ", " + t2.getTaxon().toString() + ", " + branchLength + ", " + newState.height + ", " + logLik);

		//return Pair.of(logLik + t1.logLikelihood() + t2.logLikelihood(), newState);
		return Pair.of(logw, newState); // returns the weight
	}
	
	public Pair<Double, PartialCoalescentState> priorPost(Random random, boolean isPeek)
	{
		PartialCoalescentState newState = copy(this);
		double [] logProbs = new double[(int)nChoose2(newState.numTrees())];
		List<Pair<Integer, Integer>> indices = new ArrayList<>();
		int k = 0;
		double delta = Exponential.generate(random, nChoose2(newState.numTrees()));
		newState.height = newState.height + delta;
		for (int i = 0; i < newState.trees.size(); i++)
		{
			RootedPhylogeny t1 = newState.trees.get(i);
			for (int j = i + 1; j < newState.trees.size(); j++)
			{
				RootedPhylogeny t2 = newState.trees.get(j);

				// determine branch length for each of the subtrees to the parent
				double b1 = newState.height - t1.getHeight();
				double b2 = newState.height - t2.getHeight();

				// compute the log-likelihood
				double logLik = PhyloOptions.calc.computeLoglikInStream(t1, t2, b1, b2);
				logProbs[k] = (logLik - (t1.logLikelihood() + t2.logLikelihood()));
				indices.add(Pair.of(i, j));
				k++;
			}
		}
		Multinomial.expNormalize(logProbs);
		int idx = Multinomial.sampleMultinomial(random, logProbs);
		double logw = Math.log(logProbs[idx]);

		if (isPeek) {
			return Pair.of(logw, null); 
		} else {
  		// now construct the likelihood table for the chosen root
  		// construct the likelihood table
  		int i = indices.get(idx).getLeft();
  		int j = indices.get(idx).getRight();
  		RootedPhylogeny t1 = newState.trees.get(i);
  		RootedPhylogeny t2 = newState.trees.get(j);
  		double b1 = newState.height - t1.getHeight();
  		double b2 = newState.height - t2.getHeight();
  		RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + newState.id++), t1, t2, b1, b2, newState.height, true);
  
  		// compute the log-likelihood
  		double [][] likelihoodTable = PhyloOptions.calc.computeLikelihoodTable(t1, t2, b1, b2);
  		parent.getTaxon().setLikelihoodTable(likelihoodTable);
  		parent.setLogLikelihood(PhyloOptions.calc.computeLoglik(likelihoodTable));
  
  		// remove the subtrees and add the new root
  		newState.trees.remove(t1);
  		newState.trees.remove(t2);
  		newState.trees.add(parent);
  		return Pair.of(logw, newState);
		}
	}

  public static double nChoose2(double n) { return n*(n-1)/2; }
		
	public int numTrees()
	{
		return trees.size();
	}
	
	public RootedPhylogeny getCoalescent()
	{
		if (numTrees() > 1)
			throw new RuntimeException("Coalescent has not been reached yet!");
		
		return trees.get(0);
	}
	
	public double logLikelihood()
	{
		double logLik = 0.0;
		for (RootedPhylogeny tree : trees)
		{
			logLik += tree.logLikelihood();
		}
		return logLik;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append(height + "\n");
		for (RootedPhylogeny tree : trees)
		{
			sb.append(tree.getTreeString());
		}
		return sb.toString();
	}
	
	public static void main(String [] args)
	{
		Random rand = new Random(2);
		for (int i = 0; i < 50; i++)
		{
			double param = 1.0/nChoose2(2.0);
  		double t1 = Exponential.generate(rand, param);
  		param = 1.0/nChoose2(3.0);
  		double t2 = Exponential.generate(rand, param);
  		System.out.println(t1 + ", " + t2);
		}
	}
}
