package experiments.phylo;

import java.io.File;
import java.util.List;

import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.run.Results;
import phylo.PartialCoalescentState;
import phylo.RootedPhylogeny;
import phylo.RootedPhylogenyProcessor;
import phylo.Taxon;
import spf.StreamingParticleFilter;
import util.OutputHelper;

public class PartialCoalescentStateProcessorUtil 
{
	public static double [][] computePairwiseDistances(RootedPhylogeny phylogeny, Indexer<Taxon> taxonIndexer)
	{
		int numTaxa = taxonIndexer.size();
		double [][] dd = new double[numTaxa][numTaxa];
		Counter<UnorderedPair<Taxon, Taxon>> truth = phylogeny.getPairwiseDistances(taxonIndexer);
		for (UnorderedPair<Taxon, Taxon> pair : truth)
		{
			int i = taxonIndexer.o2i(pair.getFirst());
			int j = taxonIndexer.o2i(pair.getSecond());
			double dist = truth.getCount(pair);
			dd[i][j] = dist;
			dd[j][i] = dist;
		}
		return dd;
	}
	
	private static String [] constructDistanceHeader(Indexer<Taxon> taxonIndexer)
	{
		int nTaxa = taxonIndexer.size();
		String [] header = new String[nTaxa];
		for (int i = 0; i < nTaxa; i++)
		{
			header[i] = taxonIndexer.i2o(i).getName();
		}
		return header;
	}

	public static void generateOutputSPF(int simulNo, Indexer<Taxon> taxonIndexer, StreamingParticleFilter<PartialCoalescentState> spf, RootedPhylogeny truePhylogeny)
	{
		File resultsDir = Results.getResultFolder();
		List<PartialCoalescentState> states = spf.sample().particles;
		
		String [] distanceMatrixHeader = constructDistanceHeader(taxonIndexer);
		
		// process the particle population
		RootedPhylogenyProcessor processor = new RootedPhylogenyProcessor(truePhylogeny, taxonIndexer);
		int N = states.size();
		for (PartialCoalescentState state : states)
		{
			processor.process(state, 1.0/N);
		}
		double [][] distMatrix = computePairwiseDistances(truePhylogeny, taxonIndexer);
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/pairwise-dist-spf.csv"), distanceMatrixHeader, distMatrix);
		OutputHelper.writeTableAsCSV(new File(resultsDir, "output" + simulNo + "/pairwise-dist-truth.csv"), distanceMatrixHeader, processor.getMeanPairwiseDistancesArray2D());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/height-truth.csv"), new double[]{truePhylogeny.getHeight()});
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/heights-spf.csv"), processor.getHeights());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/logZ-spf.csv"), new double[]{spf.logNormEstimate()});
		
		// output SPF statistics
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/rel-ess-spf.csv"), spf.relESS());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/ess-spf.csv"), spf.ess());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/timing-results-spf.csv"), spf.getExecutionTimes());
		OutputHelper.writeVector(new File(resultsDir, "output" + simulNo + "/num-implicit-particles-spf.csv"), spf.nImplicitParticles());
	}
}
