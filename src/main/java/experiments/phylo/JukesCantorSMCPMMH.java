package experiments.phylo;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import briefj.Indexer;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import phylo.EvolutionaryModel;
import phylo.FelsensteinPruningSystBiol2012;
import phylo.PartialCoalescentState;
import phylo.PhyloOptions;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.models.GenerateSequences;
import phylo.models.JukesCantorModel;
import pmcmc.LogZProcessor;
import pmcmc.PMCMCOptions;
import pmcmc.PMMHAlgorithm;
import pmcmc.prior.MultivariateUniformPrior;
import pmcmc.proposals.MultivariateIndependentGaussianRandomWalk;
import pmcmc.proposals.RealVectorParameters;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;

public class JukesCantorSMCPMMH implements Runnable {

	@Option(required=true)
	public Random rand = new Random(1);
	@Option(required=true)
	public int numSites = 1000;
	@Option(required=true)
	public int numTaxa = 4;
	@OptionSet(name="pmcmc")
	public PMCMCOptions pmcmcOptions = new PMCMCOptions();
	@OptionSet(name="smc")
	public SMCOptions smcOptions = new SMCOptions();
	@Option(required=true)
	public double trueMutationRate = 0.74;
	@Option(required=true)
	public double proposalVar = 0.1;
	@Option(required=true)
	public String proposalType = "priorprior";

	Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	List<Taxon> leaves = new ArrayList<Taxon>();

	private void generateData() {
		for (int n = 0; n < numTaxa; n++)
		{
			Taxon T = new Taxon("T" + n);
			leaves.add(T);
			taxonIndexer.addToIndex(T);
		}

		// instantiate the model to generate the data and the tree
		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(trueMutationRate);
		PhyloOptions.calc = new FelsensteinPruningSystBiol2012(model);
		RootedPhylogeny phylogeny = Coalescent.sampleFromCoalescent(rand, leaves);
		GenerateSequences.generateSequencesFromModel(rand, model, phylogeny, numSites);
	}
	
	@Override
	public void run()
	{
		generateData();

		EvolutionaryModel<RealVectorParameters> model = new JukesCantorModel(rand.nextDouble()*5);
		PhyloOptions.calc = new FelsensteinPruningSystBiol2012(model);

		// specify SMC algorithm
		smcOptions.random = new Random(rand.nextLong());
		smcOptions.verbose = false;
		
		SMCProblemSpecification<PartialCoalescentState> proposal = null;
		if (proposalType.equalsIgnoreCase("priorprior"))
			proposal = new PriorPriorProblemSpecification(leaves);
		else if (proposalType.equalsIgnoreCase("priorpost"))
			proposal = new PriorPostProblemSpecification(leaves);

		LogZProcessor<RealVectorParameters> logZProcessor = new LogZProcessor<>("smc");

		SMCAlgorithm<PartialCoalescentState> smc = new SMCAlgorithm<>(proposal, smcOptions);

		MultivariateIndependentGaussianRandomWalk igrwProposal = new MultivariateIndependentGaussianRandomWalk(new double[]{2.5}, new double[]{proposalVar});
		MultivariateUniformPrior uniformPrior = new MultivariateUniformPrior(new double[]{0.0, 5.0}, 1, false, true);
		
		pmcmcOptions.random = new Random(rand.nextLong());
		PMMHAlgorithm<RealVectorParameters, PartialCoalescentState> pmmh = new PMMHAlgorithm<>(model, smc, igrwProposal, uniformPrior, pmcmcOptions, null, logZProcessor, true);
		pmmh.sample();
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new JukesCantorSMCPMMH());
	}
}
