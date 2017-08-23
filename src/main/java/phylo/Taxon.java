package phylo;

import org.jblas.DoubleMatrix;

import briefj.Indexer;

public class Taxon {
	
	private String name;
	private String sequence; // *aligned
	
	private double [][] likelihoodTable;
	
	public Taxon(String name)
	{
		this.name = name;
	}

	public Taxon(String name, String sequence)
	{
		this.name = name;
		this.setSequence(sequence);
	}
	
	public String getName()
	{
		return name;
	}
	
	public void setSequence(String sequence)
	{
		this.sequence = sequence;
		initLikelihoodTable(PhyloOptions.stateIndexer);
	}
	
	public String getSequence()
	{
		return this.sequence;
	}
	
	public void setLikelihoodTable(double [][] likelihoodTable)
	{
		this.likelihoodTable = likelihoodTable;
	}
	
	public void initLikelihoodTable(Indexer<String> indexer)
	{
		if (sequence != null)
		{
			likelihoodTable = new double[sequence.length()][indexer.size()];
			String ch = "";
			for (int site = 0; site < sequence.length(); site++)
			{
				ch = sequence.charAt(site) + "";
				int index = indexer.o2i(ch);
				likelihoodTable[site][index] = 1.0;
			}
		}
		else
			throw new RuntimeException("Invalid initialization of the DP table from Taxon class");
	}
	
	public double [][] getLikelihoodTable()
	{
		/*
		if (likelihoodTable == null)
		{
			initLikelihoodTable(PhyloOptions.stateIndexer);
		}
		*/
		return this.likelihoodTable;
	}
	
	public String toString()
	{
		return this.name;
	}
	
	@Override
	public int hashCode()
	{
		return this.name.hashCode();
	}
	
	public static void main(String [] args)
	{
		Taxon t1 = new Taxon("t1");
		t1.setSequence("ACCGT");
		Indexer<String> indexer = new Indexer<String>();
		indexer.addToIndex("A", "C", "G", "T");
		t1.initLikelihoodTable(indexer);
		DoubleMatrix likelihoodTable = new DoubleMatrix(t1.getLikelihoodTable());
		System.out.println(likelihoodTable);
	}
}
