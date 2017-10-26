package phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import bayonet.math.NumericalUtils;

public class PartialCoalescentState {

	private List<RootedPhylogeny> trees = new ArrayList<RootedPhylogeny>();
	private int id = 0;
	private double height = 0.0;
	private double oldLogLik = 0.0;
	private double logLik = 0.0;
	
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
		newState.oldLogLik = src.logLik;
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
	  			phylo.setLogLikelihood(PhyloOptions.calc.computeLoglik(taxon.getLikelihoodTable()));
	  			initial.trees.add(phylo);
	  		}
		}
		return initial;
	}
	
	// Leaves the current state un-touched, including the trees (RootedPhylogeny) and everything else related to it
	public PartialCoalescentState coalesce(Random rand, boolean peek)
	{
		PartialCoalescentState newState = copy(this);

		// sample from exponential distribution
		double branchLength = Exponential.generate(rand, nChoose2(newState.trees.size()));
		//System.out.println("delta: " + branchLength);
		// select two trees to merge at random
		RootedPhylogeny t1 = newState.trees.remove(rand.nextInt(newState.trees.size()));
		RootedPhylogeny t2 = newState.trees.remove(rand.nextInt(newState.trees.size()));

		newState.height = this.height + branchLength;

		// determine branch length for each of the subtrees to the parent
		double b1 = newState.height - t1.getHeight();
		double b2 = newState.height - t2.getHeight();

		// compute the likelihood
		RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + newState.id++), t1, t2, b1, b2, newState.height, true);
		newState.trees.add(parent);

		Pair<Double, double [][]> likelihoodTable = PhyloOptions.calc.computeLikelihoodTable(t1, t2, b1, b2, peek);
		if (!peek)
			parent.getTaxon().setLikelihoodTable(likelihoodTable.getRight());
		parent.setLogLikelihood(likelihoodTable.getLeft()); // this is the likelihood for the tree

		return newState; // returns the weight
	}
	
	private PartialCoalescentState coalesce(int i, int j, double delta, boolean peek)
	{
		PartialCoalescentState newState = copy(this);
		newState.height = this.height + delta;
		RootedPhylogeny t1 = newState.trees.get(i);
		RootedPhylogeny t2 = newState.trees.get(j);

		double b1 = newState.height - t1.getHeight();
		double b2 = newState.height - t2.getHeight();

		RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + newState.id++), t1, t2, b1, b2, newState.height, true);

		newState.trees.remove(t1);
  		newState.trees.remove(t2);
  		newState.trees.add(parent);
  		
  		Pair<Double, double [][]> likelihoodTable = PhyloOptions.calc.computeLikelihoodTable(t1, t2, b1, b2, peek);
		if (!peek)
			parent.getTaxon().setLikelihoodTable(likelihoodTable.getRight());
		parent.setLogLikelihood(likelihoodTable.getLeft()); // this is the likelihood for the tree
		return newState;
	}
	
	public Pair<Double, PartialCoalescentState> priorPost(Random random, boolean peek)
	{
		double [] logProbs = new double[(int)nChoose2(this.numTrees())];
		List<Pair<Integer, Integer>> indices = new ArrayList<>();
		int k = 0;
		double delta = Exponential.generate(random, nChoose2(this.numTrees()));
		for (int i = 0; i < this.trees.size(); i++)
		{
			for (int j = i + 1; j < this.trees.size(); j++)
			{
				PartialCoalescentState tempState = this.coalesce(i, j, delta, true);
				logProbs[k] = tempState.logRatio();
				indices.add(Pair.of(i, j));
				k++;
			}
		}
		double logw = NumericalUtils.logAdd(logProbs);
		Multinomial.expNormalize(logProbs);
		int idx = Multinomial.sampleMultinomial(random, logProbs);

		if (peek) {
			return Pair.of(logw, null); 
		} else {
	  		// now construct the likelihood table for the chosen root
	  		// construct the likelihood table
	  		int i = indices.get(idx).getLeft();
	  		int j = indices.get(idx).getRight();
	  		PartialCoalescentState newState = coalesce(i, j, delta, false);
	  
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
	
	public double logRatio()
	{
		return logLik() - oldLogLik;
	}
	
	public double logLik()
	{
		this.logLik = 0.0;
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
