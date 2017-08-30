package spf;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.math.NumericalUtils;
import bayonet.smc.ParticlePopulation;


/**
 * Performs one cycle of importance sampling + resampling in a streaming fashion.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public class StreamingPropagator<S>
{
  public final ProposalWithRestart<S> proposal;
  public final SPFOptions options;
  private Consumer<S> processor = null;
  
  public StreamingPropagator(ProposalWithRestart<S> proposal, SPFOptions options)
  {
    this.proposal = proposal;
    this.options = options;
  }
  
  public void addProcessor(Consumer<S> newProcessor)
  {
    if (newProcessor == null)
      throw new RuntimeException();
    if (this.processor == null)
      this.processor = newProcessor;
    else
      this.processor = this.processor.andThen(newProcessor);
  }

  public StreamingPropagator(ProposalWithRestart<S> proposal)
  {
    this(proposal, new SPFOptions());
  }
  
  public static class PropagationResult<S>
  {
    /**
     * The virtual (implicit) particles.
     */
    public final CompactPopulation population;
    
    /**
     * The concrete samples obtained by resampling the virtual particles.
     */
    public final List<S> samples;
    
    private PropagationResult(CompactPopulation population, List<S> samples)
    {
      this.population = population;
      this.samples = samples;
    }

    public ParticlePopulation<S> getParticlePopulation()
    {
    	// resampling is carried out at every iteration for SPF
    	return ParticlePopulation.buildEquallyWeighted(samples, 0.0);
    }

  }
  
  /**
   * Performs one cycle of importance sampling + resampling.
   * @return
   */
  public PropagationResult<S> execute(int currentSmcIteration)
  {
    CompactPopulation population = new CompactPopulation(options.storeParticleWeights);
    long start = System.currentTimeMillis();
    propose(
        population, 
        options.targetedRelativeESS, 
        options.numberOfConcreteParticles, 
        options.maxNumberOfVirtualParticles,
        currentSmcIteration);
    long end = System.currentTimeMillis();
    if (options.verbose)
      System.out.println(
            "nVirtual=" + population.getNumberOfParticles() + ", "
          + "nConcrete=" + options.numberOfConcreteParticles + ", "
          + "ess=" + population.ess() + ", "
          + "relative_ess=" + (population.ess()/options.numberOfConcreteParticles) + ", "
          + "proposal_time=" + (end - start)/1000.0 + " seconds.");
    start = System.currentTimeMillis();
    double [] sortedCumulativeProbabilitiesForFinalResampling = 
        options.resamplingScheme.getSortedCumulativeProbabilities(
            options.resamplingRandom, 
            options.numberOfConcreteParticles);
    List<S> samples = null;
    if (options.storeParticleWeights) {
    	samples = resampleWithWeights(
    			population, 
    			sortedCumulativeProbabilitiesForFinalResampling, 
    			currentSmcIteration);
    } else {
	    samples = resample(
	        population, 
	        sortedCumulativeProbabilitiesForFinalResampling,
	        currentSmcIteration);
    }
    end = System.currentTimeMillis();
    System.out.println("resampling_time=" + (end - start)/1000.0 + " seconds.");
    return new PropagationResult<>(population, samples);
  }
  
  private List<S> resampleWithWeights(
	      CompactPopulation population,
	      double [] sortedCumulativeProbabilities,
	      int currentSmcIteration)
	  {
	    StreamingProposalWithWeights<S> proposal = (StreamingProposalWithWeights<S>)this.proposal;
	    if (proposal.numberOfCalls() != 0)
	    {
	      proposal = (StreamingProposalWithWeights<S>)proposal.restart();

	      if (proposal.numberOfCalls() != 0)
	        throw new RuntimeException("restart() incorrectly implemented");
	    }

	    final double logSum = population.getLogSum();
	    final int nParticles = population.getNumberOfParticles();
	    final int popAfterCollapse = sortedCumulativeProbabilities.length;
	    final List<S> result = new ArrayList<>(popAfterCollapse);
	    CompactPopulation sanityCheck = new CompactPopulation(false);

	    double normalizedPartialSum = 0.0;
	    S candidate = null;
	    int numUnique = 0;
	    int j = 0;
	    for (int i = 0; i < popAfterCollapse; i++)
	    {
	      double nextCumulativeProbability = sortedCumulativeProbabilities[i];
	      // sum normalized weights until we get to the next resampled cumulative probability
	      while (normalizedPartialSum < nextCumulativeProbability)
	      {
	    	  int before = proposal.numberOfCalls();
	    	  double logw = population.getLogWeights().get(j);
	    	  final double normalizedWeight = Math.exp(logw - logSum);
	    	  normalizedPartialSum += normalizedWeight;
	    	  if (normalizedPartialSum < nextCumulativeProbability) {
	    		  proposal.advanceStream(currentSmcIteration);
		    	  sanityCheck.insertLogWeight(logw);
	    	  } else {
	    		  Pair<Double, S> nextLogWeightSamplePair = proposal.nextLogWeightSamplePair(currentSmcIteration);
	    		  candidate = nextLogWeightSamplePair.getRight();
	    		  logw = nextLogWeightSamplePair.getLeft();
		    	  sanityCheck.insertLogWeight(nextLogWeightSamplePair.getLeft());
		    	  if (!NumericalUtils.isClose(population.getLogWeights().get(j), nextLogWeightSamplePair.getLeft(), 1e-6)) {
		    		  System.out.println("The log weights do not match: " + population.getLogWeights().get(j) + ", " + nextLogWeightSamplePair.getLeft() + ", " + j);
		    		  throw new RuntimeException();
		    	  }
	    	  }
	    	  if (proposal.numberOfCalls() != before + 1)
	    		  throw new RuntimeException("The method numberOfCalls() was incorrectly implemented in the proposal");
	    	  j++;
	    	  
	    	  //System.out.println(currentSmcIteration + ", " + j + ", " + logw);
	      }
	      if (result.size() > 0 && candidate != result.get(result.size()-1))
	    	  numUnique += 1;
	      // we have found one particle that survived the collapse
	      result.add(candidate);
	    }

	    // replay the last few calls of the proposal sequence to make sure things were indeed behaving deterministically
	    for (int i = proposal.numberOfCalls(); i < nParticles; i++) {
	    	proposal.advanceStream(currentSmcIteration);
	    	final double normalizedWeight = Math.exp(population.getLogWeights().get(i) - logSum);
	    	normalizedPartialSum += normalizedWeight;
	    	sanityCheck.insertLogWeight(population.getLogWeights().get(i));
	    }

	    System.out.println("# unique particles: " + numUnique);

	    if (sanityCheck.getLogSum() != logSum || sanityCheck.getLogSumOfSquares() != population.getLogSumOfSquares()) 
	      throw new RuntimeException("The provided proposal does not behave deterministically: " + sanityCheck.getLogSum() + " vs " + logSum);
	    
	    return result;
	  }

  
  /**
   * Perform resampling by replaying randomness to instantiate
   * concrete version of the particles that survive the resampling step.
   * 
   * @param proposal
   * @param sortedCumulativeProbabilities See ResamplingScheme
   * @return The list of resampled, equi-weighted particles
   */
  private List<S> resample(
      CompactPopulation population,
      double [] sortedCumulativeProbabilities,
      int currentSmcIteration)
  {
    ProposalWithRestart<S> proposal = this.proposal;
    if (proposal.numberOfCalls() != 0)
    {
      proposal = proposal.restart();

      if (proposal.numberOfCalls() != 0)
        throw new RuntimeException("restart() incorrectly implemented");
    }
    
    final double logSum = population.getLogSum();
    final int nParticles = population.getNumberOfParticles();
    final int popAfterCollapse = sortedCumulativeProbabilities.length;
    final List<S> result = new ArrayList<>(popAfterCollapse);
    CompactPopulation sanityCheck = new CompactPopulation(false);
    
    double normalizedPartialSum = 0.0;
    S candidate = null;
    int numUnique = 0;
    for (int i = 0; i < popAfterCollapse; i++)
    {
      double nextCumulativeProbability = sortedCumulativeProbabilities[i];
      // sum normalized weights until we get to the next resampled cumulative probability
      while (normalizedPartialSum < nextCumulativeProbability) 
      {
        int before = proposal.numberOfCalls();
        Pair<Double, S> nextLogWeightSamplePair = proposal.nextLogWeightSamplePair(currentSmcIteration);
        if (proposal.numberOfCalls() != before + 1)
          throw new RuntimeException("The method numberOfCalls() was incorrectly implemented in the proposal");
        candidate = nextLogWeightSamplePair.getRight();
        final double normalizedWeight = Math.exp(nextLogWeightSamplePair.getLeft() - logSum);
        normalizedPartialSum += normalizedWeight;
        sanityCheck.insertLogWeight(nextLogWeightSamplePair.getLeft());
      }
      if (result.size() > 0 && candidate != result.get(result.size()-1))
        	numUnique += 1;
        // we have found one particle that survived the collapse
        result.add(candidate);
    }
    
    System.out.println("num_unique=" + numUnique);
    
    // replay the last few calls of the proposal sequence to make sure things were indeed behaving deterministically
    for (int i = proposal.numberOfCalls(); i < nParticles; i++)
      sanityCheck.insertLogWeight(proposal.nextLogWeight(currentSmcIteration));
    if (sanityCheck.getLogSum() != logSum || sanityCheck.getLogSumOfSquares() != population.getLogSumOfSquares()) 
      throw new RuntimeException("The provided proposal does not behave deterministically: " + sanityCheck.getLogSum() + " vs " + logSum);
    
    return result;
  }
  
  /**
   * Grow this population by using a proposal distribution.
   * 
   * The number of particles proposed is determined as follows:
   * First, the proposal will be called at least minNumberOfParticles times.
   * Then, after minNumberOfParticles have been proposed, the growth will continue 
   * until the first of these two condition is met:
   * - maxNumberOfParticles is exceeded
   * - the relative ESS exceeds the targetedRelativeESS (here relative ESS is defined
   *   as ESS divided by minNumberOfParticles, since minNumberOfParticle will 
   *   correspond to the number of 'concrete' particles)
   * 
   * @param proposal
   * @param targetedRelativeESS
   * @param minNumberOfParticles
   * @param maxNumberOfParticles
   */
  private void  propose( 
    CompactPopulation population,
    double targetedRelativeESS,
    int minNumberOfParticles,
    int maxNumberOfParticles,
    int currentSmcIteration)
  {
    while 
        ( 
          population.getNumberOfParticles() < minNumberOfParticles || 
          (
              population.getNumberOfParticles() < maxNumberOfParticles && 
              population.ess() / minNumberOfParticles < targetedRelativeESS
          )
        )
      population.insertLogWeight(proposal.nextLogWeight(currentSmcIteration));
  }
  
}
