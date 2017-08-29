package phylo;

import org.apache.commons.math3.util.Pair;

import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import util.PartitionMetric;

public class RootedPhylogenyProcessor
{

	//private RootedPhylogeny truePhylogeny;

	private Pair<Integer, Double> minDist = new Pair<Integer, Double>(Integer.MAX_VALUE, Double.MAX_VALUE);
	private RootedPhylogeny minPhylogeny;

	private Pair<Integer, Double> maxDist = new Pair<Integer, Double>(Integer.MIN_VALUE, Double.MIN_VALUE);
	private RootedPhylogeny maxPhylogeny;

	private PartitionMetric pm;
	private int unweightedDistance = 0;
	private double weightedDistance = 0.0;
	private int numSamples = 0;
	private double normalizer = 0.0;
	
	private Indexer<Taxon> taxonIndexer;
	
	private double meanHeight = 0.0;;
	
	public RootedPhylogenyProcessor(RootedPhylogeny truePhylogeny, Indexer<Taxon> taxonIndexer)
	{
		this.taxonIndexer = taxonIndexer;
		//this.truePhylogeny = truePhylogeny;
		this.pm = new PartitionMetric(truePhylogeny);
	}
	
	Counter<UnorderedPair<Taxon, Taxon>> meanDistances = new Counter<UnorderedPair<Taxon,Taxon>>();
  public void process(PartialCoalescentState state, double weight) 
	{
		double height = state.getCoalescent().getHeight();
		//LogInfo.logs(height);
		meanHeight += weight*height;
		
		Pair<Integer, Double> distances = this.pm.computePartitionMetric(state.getCoalescent());
		unweightedDistance += weight*distances.getFirst().intValue();
		weightedDistance += weight*distances.getSecond().doubleValue();
		numSamples += 1;
		normalizer += weight;

		// compute the distance against the true phylogeny
		//LogInfo.logs(state.getCoalescent().getTreeString() + "\t weights=" + weight)
		if (distances.getSecond().doubleValue()  < minDist.getSecond().doubleValue())
		{
			minDist = distances;
			minPhylogeny = state.getCoalescent();
		}
		if (distances.getSecond().doubleValue() > maxDist.getSecond().doubleValue())
		{
			maxDist = distances;
			maxPhylogeny = state.getCoalescent();
		}
		
		// compute the pairwise distances
		Counter<UnorderedPair<Taxon, Taxon>> pairwise = state.getCoalescent().getPairwiseDistances(taxonIndexer);
		for (UnorderedPair<Taxon, Taxon> pair : pairwise.keySet())
		{
			meanDistances.incrementCount(pair, pairwise.getCount(pair)*weight);
		}
  }
	
	public int getNumSamples()
	{
		return numSamples;
	}
	
	public double getNorm()
	{
		return normalizer;
	}
	
	public Counter<UnorderedPair<Taxon, Taxon>> getMeanDistances()
	{
		return meanDistances;
	}
	
	public double getHeight()
	{
		return meanHeight/normalizer;
	}
	
	public double getAvgWeightedDistance()
	{
		return weightedDistance/normalizer;
	}
	
	public Pair<Integer, Double> getDistance(boolean min)
	{
		if (min)
			return minDist;
		return maxDist;
	}
	
	public RootedPhylogeny getMinPhylogeny()
	{
		return minPhylogeny;
	}
	
	public RootedPhylogeny getMaxPhylogeny()
	{
		return maxPhylogeny;
	}
	
	public String getInfo()
	{
		double avgUnweightedDist = (double)unweightedDistance / normalizer;
		double avgWeightedDist = weightedDistance / normalizer;
		return (avgUnweightedDist + ", " + avgWeightedDist);
	}

}
