package models;

import briefj.Indexer;

public class DNAIndexer 
{
	public static Indexer<Character> indexer;
	static {
		indexer = new Indexer<>();
		indexer.addToIndex('A');
		indexer.addToIndex('C');
		indexer.addToIndex('G');
		indexer.addToIndex('T');
	}

}
