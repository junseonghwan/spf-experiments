package phylo.models;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jblas.MatrixFunctions;

import bayonet.distributions.Multinomial;
import briefj.Indexer;
import phylo.EvolutionaryModel;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;

public class GenerateSequences 
{
	public static void generateSequencesFromModel(Random random, EvolutionaryModel model, RootedPhylogeny phylogeny, int numSites)
	{
		String rootSeq = generateSequenceFromStationaryDistribution(random, model, numSites);
		phylogeny.getTaxon().setSequence(rootSeq);
		dataGenerationRecursion(random, model, phylogeny);
	}
	
	private static void dataGenerationRecursion(Random random, EvolutionaryModel model, RootedPhylogeny tree)
	{
		String sequence = tree.getTaxon().getSequence();
		RootedPhylogeny t1 = tree.getLeftChild();
		RootedPhylogeny t2 = tree.getRightChild();
		double b1 = tree.getLeftBranchLength();
		double b2 = tree.getRightBranchLength();

		String leftSequence = generateSequence(random, model, sequence, b1);
		String rightSequence = generateSequence(random, model, sequence, b2);
		
		t1.getTaxon().setSequence(leftSequence);
		t2.getTaxon().setSequence(rightSequence);
		//t1.setLogLikelihood(PhyloOptions.calc.computeLoglik(t1.getTaxon().getLikelihoodTable()));
		//t2.setLogLikelihood(PhyloOptions.calc.computeLoglik(t2.getTaxon().getLikelihoodTable()));

		if (!t1.isLeaf())
		{
			dataGenerationRecursion(random, model, t1);
		}
		if (!t2.isLeaf())
		{
			dataGenerationRecursion(random, model, t2);
		}
	}
	
	private static String generateSequence(Random random, EvolutionaryModel model, String parentSeq, double bl)
	{
		Indexer<String> dnaIndexer = DNAIndexer.indexer;
		StringBuilder childSeq = new StringBuilder();
		double [][] P = MatrixFunctions.expm(model.getRateMatrix().mul(bl)).toArray2();
		double [] probs;
		for (int i = 0; i < parentSeq.length(); i++)
		{
			probs = P[dnaIndexer.o2i(parentSeq.charAt(i) + "")];
			int idx = Multinomial.sampleMultinomial(random, probs);
			childSeq.append(dnaIndexer.i2o(idx));
		}
		return childSeq.toString();
	}
	
	public static String generateSequenceFromStationaryDistribution(Random random, EvolutionaryModel model, int numSites)
	{
		Indexer<String> indexer = DNAIndexer.indexer;
		double [] pi = model.getStationaryDistribution();
		StringBuilder seq = new StringBuilder();
		for (int i = 0; i < numSites; i++)
		{
			int idx = Multinomial.sampleMultinomial(random, pi);
			seq.append(indexer.i2o(idx));
		}
		return seq.toString();
	}

	public static void main(String [] args)
	{
		Random rand = new Random(1022);
		List<Taxon> taxa = new ArrayList<Taxon>();
		taxa.add(new Taxon("T0"));
		taxa.add(new Taxon("T1"));
		taxa.add(new Taxon("T2"));
		taxa.add(new Taxon("T3"));
		taxa.add(new Taxon("T4"));
		taxa.add(new Taxon("T5"));

		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(rand, taxa);

		int numSites = 10;
		
		EvolutionaryModel jcModel = new JukesCantorModel(1.2);
		GenerateSequences.generateSequencesFromModel(rand, jcModel, phylogeny, numSites);

		System.out.println(phylogeny.getTreeString());
		System.out.println(phylogeny.getDataString());
	}
}
