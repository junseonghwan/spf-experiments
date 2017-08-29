package phylo;

public interface LikelihoodCalculatorInterface 
{
	public double [][] computeLikelihoodTable(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2);
	public double computeLoglik(double [][] likelihoodTable);
	public double computeLoglikInStream(RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2);
}
