package simplesmc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.SplittableRandom;
import org.apache.commons.lang3.tuple.Pair;

import bayonet.smc.ParticlePopulation;
import briefj.BriefParallel;
import phylo.FelsensteinPruningSystBiol2012;
import phylo.PhyloOptions;


/**
 * An SMC algorithm using multi-threading for proposing and suitable
 * for abstract 'SMC samplers' problems as well as more classical ones.
 * 
 * Also performs adaptive re-sampling by monitoring ESS.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <P> The type (class) of the individual particles
 */
public class SMCAlgorithm<P> extends AbstractSMCAlgorithm<P>
{
  public final SMCProblemSpecification<P> proposal;
  private final SMCOptions options;
  
  private double logZEstimate = 0.0;
  
  /**
   * This is used to ensure that the result is deterministic even in a 
   * multi-threading context: each particle index has its own unique random 
   * stream
   */
  private final Random[] randoms;
  
  private List<Double> effectiveSampleSize;
  private List<Double> logZs;
  
  private ParticleProcessor<P> processor = null;
  
  /**
   * Compute the SMC algorithm
   * 
   * @return The particle population at the last step
   */
  public ParticlePopulation<P> sample()
  {
    //System.out.println("params: " + PhyloOptions.calc.getModel().getModelParameters().asCommaSeparatedLine());

    int nSMCIterations = proposal.nIterations();
    
    timeInSeconds = new ArrayList<>(nSMCIterations); 
    effectiveSampleSize = new ArrayList<>(nSMCIterations);
    logZs = new ArrayList<>(nSMCIterations);
    long start = System.currentTimeMillis();
    ParticlePopulation<P> currentPopulation = propose(null, 0);
    effectiveSampleSize.add(currentPopulation.getESS());
    logZs.add(currentPopulation.logNormEstimate());
    if (currentPopulation.getRelativeESS() < options.essThreshold)
  	  currentPopulation = currentPopulation.resample(options.random, options.resamplingScheme);
    long end = System.currentTimeMillis();
    double samplingTime = (end-start)/1000.0;
    timeInSeconds.add(samplingTime);
    if (options.verbose) {
        System.out.println("Iter 1");
    	System.out.println("Sampling time: " + samplingTime + ", ESS: " + effectiveSampleSize.get(0));
    }
    if (processor != null)
  	  processor.process(0, currentPopulation);
    
    for (int currentIteration = 1; currentIteration < nSMCIterations; currentIteration++)
    {
    	start = System.currentTimeMillis();
    	currentPopulation = propose(currentPopulation, currentIteration);
    	effectiveSampleSize.add(currentPopulation.getESS());
    	logZs.add(currentPopulation.logNormEstimate());
    	if (currentPopulation.getRelativeESS() < options.essThreshold && currentIteration < nSMCIterations - 1)
    		currentPopulation = currentPopulation.resample(options.random, options.resamplingScheme);
    	end = System.currentTimeMillis();
    	samplingTime = (end-start)/1000.0;
    	timeInSeconds.add(samplingTime);
    	if (options.verbose) {
        	System.out.println("Iter " + (currentIteration + 1));
    		System.out.println("Sampling time: " + samplingTime + ", ESS: " + effectiveSampleSize.get(currentIteration));
    	}
    	if (processor != null)
    		processor.process(currentIteration, currentPopulation);
    }

    logZEstimate = currentPopulation.logNormEstimate();
    return currentPopulation;
  }
  
  public double logNormEstimate() { return logZEstimate; }
  public List<Double> getLogNorms() { return logZs; }
  
  /**
   * Calls the proposal options.nParticles times, form the new weights, and return the new population.
   * 
   * If the provided currentPopulation is null, use the initial distribution, otherwise, use the 
   * transition. Both are specified by the proposal object.
   * 
   * @param currentPopulation The population of particles before the proposal
   * @param currentIteration The iteration of the particles used as starting points for the proposal step
   * @return
   */
  private ParticlePopulation<P> propose(final ParticlePopulation<P> currentPopulation, final int currentIteration)
  {
    final boolean isInitial = currentPopulation == null;
    
    final double [] logWeights = new double[options.nParticles];
    @SuppressWarnings("unchecked")
    final P [] particles = (P[]) new Object[options.nParticles];
    
    BriefParallel.process(options.nParticles, options.nThreads, particleIndex ->
    {
      Pair<Double, P> proposed = isInitial ?
        proposal.proposeInitial(randoms[particleIndex]) :
        proposal.proposeNext(currentIteration, randoms[particleIndex], currentPopulation.particles.get(particleIndex));
      logWeights[particleIndex] = 
        proposed.getLeft().doubleValue() + 
        (isInitial ? 0.0 : Math.log(currentPopulation.getNormalizedWeight(particleIndex)));
      particles[particleIndex] = (proposed.getRight());
    });
    
    return ParticlePopulation.buildDestructivelyFromLogWeights(
        logWeights, 
        Arrays.asList(particles),
        isInitial ? 0.0 : currentPopulation.logScaling);
  }

  public SMCAlgorithm(SMCProblemSpecification<P> proposal, SMCOptions options)
  {
	  this(proposal, options, null);
  }

  public SMCAlgorithm(SMCProblemSpecification<P> proposal, SMCOptions options, ParticleProcessor<P> processor)
  {
    this.proposal = proposal;
    this.options = options;
    this.randoms = new Random[options.nParticles];
    SplittableRandom splitRandom = new SplittableRandom(options.random.nextLong());
    for (int i = 0; i < options.nParticles; i++)
    	this.randoms[i] = new Random(splitRandom.split().nextLong());
    
    this.processor = processor;
  }

  public List<Double> effectiveSampleSize() { return effectiveSampleSize; }
}
