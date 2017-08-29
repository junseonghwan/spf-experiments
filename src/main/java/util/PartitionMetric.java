package util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Stack;

import org.apache.commons.math3.util.Pair;

import briefj.Indexer;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import phylo.models.Coalescent;
import phylo.RootedPhylogeny.PostorderObject;

public class PartitionMetric {

	private HashMap<Taxon, ClusterTable> X = new HashMap<Taxon, PartitionMetric.ClusterTable>();
	private Indexer<Taxon> taxonIndexer = new Indexer<Taxon>();
	
	private int numInternalEdges = 0;
	private double initialCost = 0.0;

	public PartitionMetric(RootedPhylogeny originalPhylo)
	{
		// re-label the taxa in post-order sequence and compute the X-table as in (Day 1985)
		List<PostorderObject> postorder = originalPhylo.getPostorderTraversal();
		
		int leafcode = 0, R = 0;
		for (int i = 0; i < postorder.size(); i++)
		{
			PostorderObject po = postorder.get(i);
			RootedPhylogeny tree = po.getRootedPhylogeny();

			// compute the X-table
			if (tree.isLeaf())
			{
				// re-labelling achieved through the Indexer
				Taxon taxon = tree.getTaxon();
				taxonIndexer.addToIndex(taxon);
				
				if (!X.containsKey(taxon))
				{
					X.put(taxon, new ClusterTable());
				}
				R = leafcode;
				X.get(taxon).setColumn(3, leafcode);
				leafcode++;
			}
			else
			{
				this.numInternalEdges++;
				int k = i - po.getWeight();
				PostorderObject leftmost = postorder.get(k); 
				int L = X.get(leftmost.getRootedPhylogeny().getTaxon()).getColumn(3);
				int loc;

				if (i+1 == postorder.size())
				{
					loc = R;
				}
				else
				{
					PostorderObject next = postorder.get(i+1);
  				if (next.getWeight() == 0)
  				{
  					loc = R;
  				}
  				else
  				{
  					loc = L;
  				}
				}
				
				Taxon taxon = taxonIndexer.i2o(loc);
				X.get(taxon).setColumn(1, L);
				X.get(taxon).setColumn(2, R);
				
				double we = po.getBranchLength();
				X.get(taxon).setCost(we);
				initialCost += we;
			}

		}

	}

	public Pair<Integer, Double> computePartitionMetric(RootedPhylogeny other)
	{
		// compute both the unweighted and weighted partition metric -- return both
		int numIdentical = 0;
		int internalEdges = 0;
		double dist = initialCost;

		// 1. get postorder sequence traversal of other
		List<PostorderObject> postorder = other.getPostorderTraversal();
		
		// 2. for each, internal node, get L and R -- use stack
		Stack<Cluster> clusters = new Stack<Cluster>();
		
		// 3. check for membership by checking the row L and R of X
		for (int i = 0; i < postorder.size(); i++)
		{
			PostorderObject node = postorder.get(i);
			RootedPhylogeny tree = node.getRootedPhylogeny();
			Taxon taxon = tree.getTaxon();

			if (tree.isLeaf())
			{
				int j = this.taxonIndexer.o2i(taxon);
				Cluster cluster = new Cluster(j, j, 1, 1);
				clusters.push(cluster);
			}
			else
			{
				internalEdges++;
				
				Cluster cluster0 = new Cluster(Integer.MAX_VALUE, 0, 0, 1);
				int w = node.getWeight();
				while (w != 0)
				{
					Cluster cluster = clusters.pop();
					cluster0.L = Math.min(cluster0.L, cluster.L);
					cluster0.R = Math.max(cluster0.R, cluster.R);
					cluster0.N = cluster0.N + cluster.N;
					cluster0.W = cluster0.W + cluster.W;
					w = w - cluster.W;

					if (w < 0)
						throw new RuntimeException();
				}

				clusters.push(cluster0);

				// check if this new cluster exists in the table X
				// 4. TODO: bring in edge weight (branch length) somehow
				Pair<Boolean, Double> member = this.checkClusterMembership(cluster0);
				double we = node.getBranchLength();
				if (member.getFirst().booleanValue())
				{
					numIdentical++;
					dist += Math.abs(member.getSecond().doubleValue() - we);
					dist -= member.getSecond().doubleValue(); // subtract the cost of the original tree because it is being counted twice
				}
				else
				{
					dist += we;
				}

			}
			
		}
				
		int unweighted = this.numInternalEdges + internalEdges - 2*numIdentical; 
		return new Pair<Integer, Double>(unweighted, dist);
	}
	
	// returns membership test (boolean ) as the first item and the weight of the internal edge as the second item
	private Pair<Boolean, Double> checkClusterMembership(Cluster c)
	{
		if (c.N == (c.R - c.L + 1))
		{
  		Taxon leftmostTaxon = this.taxonIndexer.i2o(c.L);
  		boolean leftMembership = this.X.containsKey(leftmostTaxon);
  		if (leftMembership)
  		{
  			ClusterTable clusterTable = this.X.get(leftmostTaxon);
  			if (clusterTable.getColumn(1) == c.L && clusterTable.getColumn(2) == c.R) 
  				return new Pair<Boolean, Double>(true, clusterTable.getCost());
  		}

  		Taxon rightmostTaxon = this.taxonIndexer.i2o(c.R);
  		boolean rightMembership = this.X.containsKey(rightmostTaxon);
  		if (rightMembership)
  		{
  			ClusterTable clusterTable = this.X.get(rightmostTaxon);
  			if (clusterTable.getColumn(1) == c.L && clusterTable.getColumn(2) == c.R) 
  				return new Pair<Boolean, Double>(true, clusterTable.getCost());
  		}
  		
		}

		return new Pair<Boolean, Double>(false, 0.0);
	}


	
	public static class Cluster
	{
		public int L;
		public int R;
		public int N;
		public int W;
		
		public Cluster(int L, int R, int N, int W)
		{
			this.L = L;
			this.R = R;
			this.N = N;
			this.W = W;
		}
	}
	
	public static class ClusterTable
	{
		private int [] X = new int[4];
		private double cost = 0.0;
		public void setCost(double cost)
		{
			this.cost = cost;
		}
		public double getCost()
		{
			return this.cost;
		}
		
		public void setColumn(int col, int value)
		{
			X[col] = value;
		}
		
		public int getColumn(int col)
		{
			return X[col];
		}
	}
	
	public static void main(String [] args)
	{
		List<Taxon> leaves = new ArrayList<Taxon>();
		leaves.add(new Taxon("T1"));
		leaves.add(new Taxon("T2"));
		leaves.add(new Taxon("T3"));
		leaves.add(new Taxon("T4"));
		//leaves.add(new Taxon("T5"));
		//leaves.add(new Taxon("T6"));
		RootedPhylogeny original = Coalescent.sampleFromCoalescent(new Random(1), leaves);
		//RootedPhylogeny other = GenerateData.sampleRootedPhylogeny(new Random(2), leaves, rate);
		RootedPhylogeny other = Coalescent.sampleFromCoalescent(new Random(3), leaves);
		System.out.println(original.getTreeString());
		System.out.println(other.getTreeString());
		
		PartitionMetric pm = new PartitionMetric(original);
		
		// compute the distance (un-weighted for now)
		Pair<Integer, Double> distances = pm.computePartitionMetric(other);
		System.out.println("unweighted dist=" + distances.getFirst().intValue() + ", weighted dist=" + distances.getSecond().doubleValue());
	}
	
}
