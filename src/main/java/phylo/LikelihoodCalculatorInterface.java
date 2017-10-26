package phylo;

import org.apache.commons.lang3.tuple.Pair;

public interface LikelihoodCalculatorInterface 
{
	// returns the likelihood table and the loglikelihood for the merged tree with t1 and t2 as children and b1 b2 as branch lengths 
	public Pair<Double, double [][]> computeLikelihoodTable(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2, boolean peek);
	public double computeLoglik(double [][] likelihoodTable);
	// same as computeLikelihoodTable except that it does not return the likelihood table and hence, it does not allocate memory for the table
	public double computeLoglikInStream(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2);
	
	public static void computeDataLogLikTable(LikelihoodCalculatorInterface peeling, RootedPhylogeny root)
	{
		if (!root.isLeaf())
		{
			// recurse to ensure children have likelihood table constructed 
			computeDataLogLikTable(peeling, root.getLeftChild());
			computeDataLogLikTable(peeling, root.getRightChild());

			// construct the likelihood table
			double [][] likTbl = peeling.computeLikelihoodTable(root.getLeftChild(), root.getRightChild(), root.getLeftBranchLength(), root.getRightBranchLength(), false).getRight();
			root.getTaxon().setLikelihoodTable(likTbl);
		}
	}

}
