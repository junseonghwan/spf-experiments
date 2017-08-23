package phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import bayonet.distributions.Exponential;

public class PartialCoalescentState {

	private List<RootedPhylogeny> trees = new ArrayList<RootedPhylogeny>();
	private int id = 0;
	private double height = 0.0;
	
	private static PartialCoalescentState initial = null;

	private PartialCoalescentState()
	{
		// prevent initialization outside of this class
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
  			initial.trees.add(phylo);
  		}
		}
		return initial;
	}
	
	// Leaves the current state un-touched, including the trees (RootedPhylogeny) and everything else related to it
	public PartialCoalescentState coalesce(Random rand)
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
		
		RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + newState.id++), t1, t2, b1, b2, newState.height);
		newState.trees.add(parent);

		// compute the likelihood
		double [][] likelihoodTable = PhyloOptions.calc.computeLikelihoodTable(t1, t2, b1, b2);
		parent.getTaxon().setLikelihoodTable(likelihoodTable);
		double logLik = PhyloOptions.calc.computeLoglik(likelihoodTable);
		parent.setLogLikelihood(logLik);

		//System.out.println(newState.trees.size() + ", " + t1.getTaxon().toString() + ", " + t2.getTaxon().toString() + ", " + branchLength + ", " + newState.height + ", " + logLik);

		return newState;
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
		/*
		StringBuilder sb = new StringBuilder();
		for (RootedPhylogeny tree : trees)
		{
			sb.append(tree.getTreeString());
		}
		return sb.toString();
		*/
		return height + "";
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
