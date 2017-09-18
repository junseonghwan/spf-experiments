package spf;


import bayonet.math.NumericalUtils;


/**
 * Keeps the sufficient statistics of a potentially large collection of 
 * weighted particles.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class CompactPopulation
{
	private boolean storeWeights = false;
	private double [] logWeights;
  private int nParticles = 0;
  
	public CompactPopulation() { 
		this.storeWeights = false;
	}
	
	public CompactPopulation(boolean storeWeights, int maxNumberOfVirtualParticles) { 
		this.storeWeights = storeWeights;
		if (storeWeights) {
			this.logWeights = new double[maxNumberOfVirtualParticles];
		}
	}
  /**
   * Sum of the unnormalized weights (log scale)
   */
  private double logSum = Double.NEGATIVE_INFINITY;
  
  /**
   * Sum of the squares-of-unnormalized-weights (log scale)
   */
  private double logSumOfSquares = Double.NEGATIVE_INFINITY;
  
  public int getNumberOfParticles() 
  { 
    return nParticles; 
  }
  
  public synchronized void insertLogWeight(double logWeight)
  {
    if (storeWeights) {
  	  logWeights[nParticles] = logWeight;
    }

    nParticles++;
    logSum          = NumericalUtils.logAdd(logSum,              logWeight);
    logSumOfSquares = NumericalUtils.logAdd(logSumOfSquares, 2 * logWeight);
    //System.out.println("logw: " + logWeight + ", " + logSum);
  }

  public synchronized void insertLogWeight(double logWeight, int particleIdx)
  {
    nParticles++;
    logSum          = NumericalUtils.logAdd(logSum,              logWeight);
    logSumOfSquares = NumericalUtils.logAdd(logSumOfSquares, 2 * logWeight);
    //System.out.println("logw: " + logWeight + ", " + logSum);
    
    if (storeWeights) {
    	logWeights[particleIdx] = logWeight;
    }
  }
  
  /**
   * @return The effective sampling size (ESS).
   */
  public double ess()
  {
    return Math.exp(2 * logSum - logSumOfSquares);
  }
  
  /**
   * @return Standard importance sampling estimator of normalization (log scale)
   */
  public double logZEstimate()
  {
    return logSum - Math.log(nParticles);
  }

  public double getLogSum()
  {
    return logSum;
  }

  public double getLogSumOfSquares()
  {
    return logSumOfSquares;
  }
  
  /*
  public double [] getLogWeights() { 
	  return logWeights; 
  }
  */
  public double getLogWeight(int idx) {
	  return logWeights[idx];
  }
}