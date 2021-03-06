package phylo;

import java.util.ArrayList;
import java.util.List;

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
	private List<Double> heights = new ArrayList<>();
	private List<Double> weights = new ArrayList<>();

	double [][] dd = null;
	Counter<UnorderedPair<Taxon, Taxon>> meanPairwiseDistances = new Counter<UnorderedPair<Taxon,Taxon>>();

	public RootedPhylogenyProcessor(RootedPhylogeny truePhylogeny, Indexer<Taxon> taxonIndexer)
	{
		this.taxonIndexer = taxonIndexer;
		//this.truePhylogeny = truePhylogeny;
		this.pm = new PartitionMetric(truePhylogeny);
	}

	public void process(PartialCoalescentState state, double weight) 
	{
		double height = state.getCoalescent().getHeight();
		heights.add(height);
		weights.add(weight);
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
			meanPairwiseDistances.incrementCount(pair, pairwise.getCount(pair)*weight);
		}
	}
	
	public List<Double> getHeights() { return heights; }
	public List<Double> getWeights() { return weights; }
	
	public int getNumSamples()
	{
		return numSamples;
	}
	
	public double getNorm()
	{
		return normalizer;
	}
	
	public Counter<UnorderedPair<Taxon, Taxon>> getMeanPairwiseDistances()
	{
		return meanPairwiseDistances;
	}
	
	private void constructMeanPairwiseDistancesAsArray()
	{
		int nTaxa = taxonIndexer.size();
		dd = new double[nTaxa][nTaxa];
		for (UnorderedPair<Taxon, Taxon> key : meanPairwiseDistances.keySet())
		{
			int i = taxonIndexer.o2i(key.getFirst());
			int j = taxonIndexer.o2i(key.getSecond());
			dd[i][j] = meanPairwiseDistances.getCount(key);
			dd[j][i] = dd[i][j];
		}
	}
	
	public double [][] getMeanPairwiseDistancesArray2D()
	{
		if (dd == null)
			constructMeanPairwiseDistancesAsArray();
		return dd;
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
