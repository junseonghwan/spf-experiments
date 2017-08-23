package phylo.models;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import phylo.RootedPhylogeny;
import phylo.Taxon;
import bayonet.distributions.Exponential;

public class Coalescent 
{
	/*
	 * Sample a phylogenetic X-tree with the labels given by the leaves parameter and the branch length drawn from exp(rate)
	 * The tree is sampled from the coalescent, i.e., randomly choose a pair to merge until onely one tree remains
	 * rand: randomness used in generating the data
	 * leaves: taxa at the leaves of the tree
	 * rate: branch length parameter for the exponential distribution
	 */
	public static RootedPhylogeny sampleFromCoalescent(Random rand, List<Taxon> leaves)
	{
		List<RootedPhylogeny> trees = new ArrayList<>();

		// create a single node trees out of the leaves
		for (Taxon taxon : leaves)
		{
			trees.add(new RootedPhylogeny(taxon));
		}
		
		int id = 0;
		double height = 0.0;
		while (trees.size() > 1)
		{
			int n = trees.size();
			// select trees to merge
			RootedPhylogeny t1 = trees.remove(rand.nextInt(trees.size()));
			RootedPhylogeny t2 = trees.remove(rand.nextInt(trees.size()));

			// sample from exp(rate)
			double branchLength = Exponential.generate(rand, nChoose2(n));
			height += branchLength;

			// determine branch length for each of the subtrees to its parent
			double b1 = height - t1.getHeight();
			double b2 = height - t2.getHeight();
			RootedPhylogeny parent = new RootedPhylogeny(new Taxon("t" + id++), t1, t2, b1, b2, height);
			trees.add(parent);
		}
		
		return trees.get(0);
	}
	
	public static double nChoose2(int n)
	{
		return n*(n-1)/2.0;
	}

}