package experiments.phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import bayonet.smc.ParticlePopulation;
import briefj.Indexer;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningAlgorithm;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantorModel;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;

public class JukesCantorSMC implements Runnable 
{
	@Option(required=true)
	public Random rand = new Random(1);
	@Option(required=true)
	public int numTaxa = 20;
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numSimulations = 5;
	@Option(required=true)
	public double mutationRate = 1.0;
	@Option(required=true)
	public String proposalType = "priorprior";
	
	@OptionSet(name="smc")
	public SMCOptions options = new SMCOptions();

	private Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	private List<Taxon> leaves = new ArrayList<Taxon>();
	private RootedPhylogeny phylogeny = null;

	@Override
	public void run() 
	{
		// generate the data once, run simulations multiple times	
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(mutationRate);
		PhyloOptions.calc = new FelsensteinPruningAlgorithm(model);

		// generate the data and the tree
		long seed = rand.nextLong();
		System.out.println("seed: " + seed);
		Random random = new Random(seed);
		phylogeny = Coalescent.sampleFromCoalescent(random, leaves);
		GenerateSequences.generateSequencesFromModel(random, model, phylogeny, numSites);
		
		SMCProblemSpecification<PartialCoalescentState> proposal = null;
		if (proposalType.equalsIgnoreCase("priorprior"))
			proposal = new PriorPriorProblemSpecification(leaves);
		else if (proposalType.equalsIgnoreCase("priorpost"))
			proposal = new PriorPostProblemSpecification(leaves);

		for (int i = 0; i < numSimulations; i++)
		{
			random = new Random(rand.nextLong());
			System.out.println("Simulation: " + (i+1));
			simulation(random, proposal, i+1);
		}
	}
	
	public void simulation(Random random, SMCProblemSpecification<PartialCoalescentState> proposal, int simulNo)
	{
		options.random = new Random(random.nextLong());
		options.verbose = true;

		SMCAlgorithm<PartialCoalescentState> smc = new SMCAlgorithm<>(proposal, options);
		ParticlePopulation<PartialCoalescentState> pop = smc.sample();
		PartialCoalescentStateProcessorUtil.generateOutputSMC(simulNo, taxonIndexer, pop, smc, phylogeny);
	}


	public static void main(String[] args) 
	{
		Mains.instrumentedRun(args, new JukesCantorSMC());
	}

}
