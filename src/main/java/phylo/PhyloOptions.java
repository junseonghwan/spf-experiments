package phylo;

import phylo.models.DNAIndexer;
import briefj.Indexer;

public class PhyloOptions {

	public static LikelihoodCalculatorInterface calc;
	public static Indexer<String> stateIndexer = DNAIndexer.indexer;
}
