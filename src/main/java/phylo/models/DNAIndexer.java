package phylo.models;

import briefj.Indexer;

public class DNAIndexer 
{
	public static Indexer<String> indexer;
	
	static {
		indexer = new Indexer<>();
		indexer.addToIndex("A");
		indexer.addToIndex("C");
		indexer.addToIndex("G");
		indexer.addToIndex("T");
	}

}
