package simplesmc.pmcmc;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;

import blang.annotations.FactorComponent;
import blang.factors.Factor;
import simplesmc.AbstractSMCAlgorithm;


/**
 * A likelihood (factor in a factor graph) approximated using
 * PMCMC (particle MCMC)
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <P>
 */
public class PMCMCFactor<P> implements Factor
{
  /**
   * The parameters.
   * 
   * They are of type WithSignature because all we 
   * care about them is to monitor if they have changed 
   * (we want to make sure that if a move is rejected, 
   * we do not recompute Z again).
   */
  @FactorComponent
  public final WithSignature params;
  
  private final AbstractSMCAlgorithm<P> smcAlgorithm;

  /**
   * Used to remember the likelihood values associated with 
   * the recent calls.
   */
  private final Cache<Long, Double> cache = CacheBuilder.newBuilder().maximumSize(10).build();
  
  @Override
  public double logDensity()
  {
    final long currentSignature = params.signature();
    
    if (cache.getIfPresent(currentSignature) != null)
      return cache.getIfPresent(currentSignature);
    
    smcAlgorithm.sample();
    final double result = smcAlgorithm.logNormEstimate();
    cache.put(currentSignature, result);
    
    return result;
  }

  public PMCMCFactor(WithSignature params, AbstractSMCAlgorithm<P> smcAlgorithm)
  {
    this.params = params;
    this.smcAlgorithm = smcAlgorithm;
  }
}
