package phylo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.Stack;

import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

import phylo.models.Coalescent;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

// rooted, directed tree (direction towards the leaves)
// the nbrs of a node are children (parent is not a nbr as the edge is directed downwards)
// the difference in the height yields the branch length
public class RootedPhylogeny
{
  private Taxon taxon = null;
  // store the left and right tree, as well as the branch length
  private Pair<RootedPhylogeny, Double> t1 = null;
  private Pair<RootedPhylogeny, Double> t2 = null;

	private double height;
	private String treeString = null;
	private String dataString = null;

	private double logLik = 0.0;

	public RootedPhylogeny(Taxon taxon)
	{
		this.taxon = taxon;
		this.height = 0.0;
	}

	// two subtrees and the branch lengths from this node to the subtrees
	public RootedPhylogeny(Taxon taxon, RootedPhylogeny t1, RootedPhylogeny t2, double b1, double b2, double height)
	{
		this(taxon);
		this.t1 = Pair.create(t1, b1);
		this.t2 = Pair.create(t2, b2);
		this.height = height;
	}
	
	public boolean isLeaf()
	{
		return (t1 == null && t2 == null);
	}

	public List<RootedPhylogeny> getChildren()
	{
		List<RootedPhylogeny> children = new ArrayList<RootedPhylogeny>();
		children.add(t1.getFirst());
		children.add(t2.getFirst());
		return children;
	}
	
  public double getHeight()
  {
  	return this.height;
  }
  
  public double getLeftBranchLength()
  {
  	return this.t1.getSecond();
  }
  
  public double getRightBranchLength()
  {
  	return this.t2.getSecond();
  }

  public RootedPhylogeny getLeftChild()
  {
  	return this.t1.getFirst();
  }

  public RootedPhylogeny getRightChild()
  {
  	return this.t2.getFirst();
  }

  public Taxon getTaxon()
  {
  	return this.taxon;
  }
  
  public List<RootedPhylogeny> getPreorderTraversal()
  {
  	List<RootedPhylogeny> traversal = new ArrayList<>();
  	traversalHelper(this, traversal);
  	return traversal;
  }
  
  public Counter<UnorderedPair<Taxon, Taxon>> getPairwiseDistances(Indexer<Taxon> indexer)
  {
  	Counter<UnorderedPair<Taxon, Taxon>> distance = new Counter<UnorderedPair<Taxon,Taxon>>();
  	
  	// initialize the distances to 0
  	for (int i = 0; i < indexer.size(); i++)
  	{
  		Taxon taxon = indexer.i2o(i);
  		for (int j = i + 1; j < indexer.size(); j++)
  		{
  			Taxon other = indexer.i2o(j);
  			distance.setCount(UnorderedPair.of(taxon, other), 0.0);
  		}
  	}
  	
  	// one recursion approach
  	pairwiseDistanceRecursion(this, distance, 0.0, indexer);
  	
  	// this seems to be too slow... implement a faster algorithm for extracting pairwise distances
  	//pairwiseDistanceHelper(this, pairwise);
  	return distance;
  }
    
  private Set<Taxon> pairwiseDistanceRecursion(RootedPhylogeny tree, Counter<UnorderedPair<Taxon, Taxon>> distance, double bl, Indexer<Taxon> indexer)
  {
  	if (tree.isLeaf())
  	{
			Taxon taxon = tree.getTaxon();
  		for (int i = 0; i < indexer.size(); i++)
  		{
  			Taxon other = indexer.i2o(i);
  			
  			if (taxon.getName().equals(other.getName()))
  				continue;
  			
  			// surely this taxon will not be in the set as all leaf taxon are visited exactly once
  			UnorderedPair<Taxon, Taxon> unorderedPair = UnorderedPair.of(taxon, other);
  			double count = distance.getCount(unorderedPair) + bl;
  			distance.setCount(unorderedPair, count);
  		}
    	Set<Taxon> set = new HashSet<Taxon>();
  		set.add(taxon);
  		return set;
  	}
  	else
  	{
  		Set<Taxon> set = pairwiseDistanceRecursion(tree.getLeftChild(), distance, tree.getLeftBranchLength(), indexer);
  		Set<Taxon> set2 = pairwiseDistanceRecursion(tree.getRightChild(), distance, tree.getRightBranchLength(), indexer);
  		set.addAll(set2);
  		for (Taxon taxon : set)
  		{
    		for (int i = 0; i < indexer.size(); i++)
    		{
    			Taxon other = indexer.i2o(i);
    			
    			if (!set.contains(other))
    			{
      			UnorderedPair<Taxon, Taxon> unorderedPair = UnorderedPair.of(taxon, other);
      			double count = distance.getCount(unorderedPair) + bl;
      			distance.setCount(unorderedPair, count);
    			}
    		}
  		}
  		return set;
  	}
  }
  
  private void pairwiseDistanceHelper(RootedPhylogeny tree, Counter<UnorderedPair<Taxon, Taxon>> pairwise)
  {
  	if (tree.isLeaf())
  		return;
  	
  	RootedPhylogeny left = tree.getLeftChild();
  	List<RootedPhylogeny> leftLeaves = left.getPreorderTraversal();
  	
  	RootedPhylogeny right = tree.getRightChild();
  	List<RootedPhylogeny> rightLeaves = right.getPreorderTraversal();
  	
  	for (RootedPhylogeny l : leftLeaves)
  	{
  		if (!l.isLeaf())
  			continue;
  		
  		for (RootedPhylogeny r: rightLeaves)
  		{
  			if (!r.isLeaf())
  				continue;
  			
  			pairwise.setCount(UnorderedPair.of(l.getTaxon(), r.getTaxon()), tree.getHeight());
  		}
  	}
  	
  	// recurse
  	pairwiseDistanceHelper(left, pairwise);
  	pairwiseDistanceHelper(right, pairwise);
  }
  
  private void traversalHelper(RootedPhylogeny tree, List<RootedPhylogeny> traversal)
  {
  	traversal.add(tree);
  	treeString += "(" + tree.toString();

  	if (tree.isLeaf())
  	{
  		treeString += ")";
  		return;
  	}
  	
  	traversalHelper(tree.t1.getFirst(), traversal);
  	traversalHelper(tree.t2.getFirst(), traversal);
  	treeString += ")";
  }
  
  public void setLogLikelihood(double logLik)
  {
  	this.logLik = logLik;
  }
  
  public double logLikelihood()
  {
  	return logLik;
  }

  @Override
  public String toString()
  {
  	return taxon.getName() + ":" + Precision.round(getHeight(), 6);
  }
  
  public String getDataString()
  {
  	if (dataString == null)
  	{
    	dataString = "";
      List<RootedPhylogeny> preorder = getPreorderTraversal();
      for (RootedPhylogeny phylo : preorder)
      {
      	if (!phylo.isLeaf())
      	{
        	RootedPhylogeny left = phylo.getLeftChild();
        	RootedPhylogeny right = phylo.getRightChild();
        	dataString += "(" + phylo.toString() + ":" + phylo.getTaxon().getSequence() + ") -> ";
        	dataString += "(" + left.toString() + ":" + left.getTaxon().getSequence() + ")\n";
        	dataString += "(" + phylo.toString() + ":" + phylo.getTaxon().getSequence() + ") -> ";
        	dataString += "(" + right.toString() + ":" + right.getTaxon().getSequence() + ")\n";
      	}
      }
  	}
  	return dataString;
  }
  
  public String getNewickFormat()
  {
  	return (getNewickFormatHelper(0.0) + ";");
  }
  
  private String getNewickFormatHelper(double bl)
  {
  	if (this.isLeaf())
  		return this.taxon.toString() + ":" + bl;
  	else
  	{
  		StringBuilder sb = new StringBuilder();
  		sb.append("(");
  		sb.append(this.getLeftChild().getNewickFormatHelper(this.getLeftBranchLength()));
  		sb.append(",");
  		sb.append(this.getRightChild().getNewickFormatHelper(this.getRightBranchLength()));
  		sb.append(")");
  		sb.append(":");
  		sb.append(bl);
  		return sb.toString();
  	}
  }
  
  public String getTreeString()
  {
  	if (treeString == null)
  	{
  		treeString = "";
  		this.getPreorderTraversal();
  	}
  	return treeString;
  }
  
  public Hashtable<Taxon, RootedPhylogeny> getCluster()
  {
  	Hashtable<Taxon, RootedPhylogeny> cluster = new Hashtable<Taxon, RootedPhylogeny>();
  	
  	getClusterHelper(this, cluster);
  	// traverse the tree and store the leaves
  	return cluster;
  }

  private void getClusterHelper(RootedPhylogeny phylogeny, Hashtable<Taxon, RootedPhylogeny> cluster)
  {
  	if (phylogeny.isLeaf())
  	{
  		cluster.put(phylogeny.getTaxon(), phylogeny);
  		return;
  	}

  	getClusterHelper(phylogeny.getLeftChild(), cluster);
  	getClusterHelper(phylogeny.getRightChild(), cluster);
  }

  public List<PostorderObject> getPostorderTraversal()
  {
  	List<PostorderObject> postorder = new ArrayList<PostorderObject>();
  	postorderTraversalHelper(this, postorder, 0.0);

  	return postorder;
  }

  private int postorderTraversalHelper(RootedPhylogeny phylo, List<PostorderObject> postorder, double bl)
  {
  	int leftCount = 0, rightCount = 0;
  	if (!phylo.isLeaf())
  	{
    	leftCount = postorderTraversalHelper(phylo.getLeftChild(), postorder, phylo.getLeftBranchLength()) + 1;
    	rightCount = postorderTraversalHelper(phylo.getRightChild(), postorder, phylo.getRightBranchLength()) + 1;
  	}

		postorder.add(new PostorderObject(phylo, bl, leftCount + rightCount));
		return (leftCount + rightCount);
  }
  
  public static class PostorderObject
  {
  	public RootedPhylogeny phylo;
  	public double branchLength;
  	public int weight;
  	
  	public PostorderObject(RootedPhylogeny phylo, double bl, int w)
  	{
  		this.phylo = phylo;
  		this.branchLength = bl;
  		this.weight = w;
  	}
  	
  	public RootedPhylogeny getRootedPhylogeny()
  	{
  		return phylo;
  	}
  	
  	public int getWeight()
  	{
  		return weight;
  	}
  	
  	public double getBranchLength()
  	{
  		return this.branchLength;
  	}
  }
  
  public static void main(String [] args)
  {
  	int T = 20;
  	List<Taxon> leaves = new ArrayList<Taxon>();
  	Indexer<Taxon> indexer = new Indexer<Taxon>();
  	for (int t = 0; t < T; t++)
  	{
    	Taxon taxon = new Taxon("T" + t);
    	leaves.add(taxon);
    	indexer.addToIndex(taxon);
  	}
  	
  	RootedPhylogeny phylo = Coalescent.sampleFromCoalescent(new Random(1), leaves);
  	System.out.println(phylo.getTreeString());

  	/*
  	List<PostorderObject> postorder = phylo.getPostorderTraversal();
  	for (PostorderObject tree : postorder)
  	{
  		System.out.println(tree.phylo.toString() + "\t" + tree.branchLength + "\t" + tree.weight);
  	}
  	*/
  	
  	System.out.println(phylo.getNewickFormat());
  	
  	long start = System.currentTimeMillis();
  	System.out.println(phylo.getPairwiseDistances(indexer));
  	long end = System.currentTimeMillis();
  	double time = (end - start)/1000.0;
  	System.out.println("pairwise distance computation time: " + time + " sec");
  	
  }

}
