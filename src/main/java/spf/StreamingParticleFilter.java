package spf;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.smc.ParticlePopulation;
import simplesmc.AbstractSMCAlgorithm;
import simplesmc.ParticleProcessor;
import simplesmc.SMCProblemSpecification;
import spf.StreamingPropagator.PropagationResult;

/**
 * Re-implementation of StreamingBootstrapFilter
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 * @param <P>
 */
public class StreamingParticleFilter<P> extends AbstractSMCAlgorithm<P>
{
  public final SMCProblemSpecification<P> problemSpec;
  private final SPFOptions options;
  private double logZ = 0.0;
  private List<Double> nImplicitParticles;
  private List<Double> relESS;
  private ParticleProcessor<P> processor;

  public StreamingParticleFilter(SMCProblemSpecification<P> problemSpec, SPFOptions options)
  {
    this(problemSpec, options, null);
  }
  
  public StreamingParticleFilter(SMCProblemSpecification<P> problemSpec, SPFOptions options, ParticleProcessor<P> processor)
  {
    this.problemSpec = problemSpec;
    this.options = options;
    this.processor = processor;
  }

  public ParticlePopulation<P> sample()
  {
    int nSMCIterations = problemSpec.nIterations();

    // instantiate new arraylist each time
    nImplicitParticles = new ArrayList<>(nSMCIterations);
    relESS = new ArrayList<>();
    timeInSeconds = new ArrayList<>(nSMCIterations); 
    long start = 0, end = 0;

    // initial distribution
    start = System.currentTimeMillis();
    StreamingBootstrapProposal proposal = getInitialDistributionProposal();
    StreamingPropagator<P> propagator = new StreamingPropagator<P>(proposal, options);
    PropagationResult<P> propResults = propagator.execute(0);
    logZ = propResults.population.logZEstimate();
    end = System.currentTimeMillis();
    nImplicitParticles.add((double)propResults.population.getNumberOfParticles());
    relESS.add(propResults.population.ess()/options.numberOfConcreteParticles);
    timeInSeconds.add((end - start)/1000.);
    if (processor != null)
  	  processor.process(0, propResults.getParticlePopulation());

    // recursion
    for (int i = 1; i < nSMCIterations; i++)
    {
      start = System.currentTimeMillis();
      proposal = new StreamingBootstrapProposal(options.mainRandom.nextLong(), propResults.getParticlePopulation());
      propagator = new StreamingPropagator<>(proposal, options);
      propResults = propagator.execute(i);
      logZ += propResults.population.logZEstimate();
      end = System.currentTimeMillis();
      nImplicitParticles.add((double)propResults.population.getNumberOfParticles());
      relESS.add(propResults.population.ess()/options.numberOfConcreteParticles);
      timeInSeconds.add((end - start)/1000.);
      if (processor != null)
    	  processor.process(i, propResults.getParticlePopulation());
    }

    return propResults.getParticlePopulation();
  }
  
  public double logNormEstimate() { return logZ; }
  public List<Double> nImplicitParticles() { return nImplicitParticles; }
  public List<Double> relESS() { return relESS; }
  
  private StreamingBootstrapProposal getInitialDistributionProposal()
  {
    return new StreamingBootstrapProposal(options.mainRandom.nextLong(), null);
  }

  /**
   * Adapt the more abstract machinery of the lazy proposal/propagator to the
   *  simpler bootstrap filter case.
   *  
   * This type is used both for the initialization (in which case oldLatents is null), and
   * for the recursion steps. One instance is created at each SMC generation.
   *  
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  private class StreamingBootstrapProposal implements ProposalWithRestart<P>
  {
    private final long seed;
    private final PermutationStream permutationStream;
    private final Random random;
    
    ParticlePopulation<P> previousPopulation;
    private int nCalls = 0;

    private StreamingBootstrapProposal(long seed, ParticlePopulation<P> previousPopulation)
    {
      this.seed = seed;
      this.random = new Random(seed * 171);
      this.previousPopulation = previousPopulation;
      this.permutationStream = previousPopulation == null ? null : new PermutationStream(previousPopulation.nParticles(), random);
    }

    @Override
    public Pair<Double, P> nextLogWeightSamplePair(int currentSmcIteration)
    {
      // terminology: old means the SMC generation before current (null if we are doing initial)
      //              cur means the current SMC generation
      Pair<Double, P> curLatent = isInitial() 
          ? problemSpec.proposeInitial(random)
          : problemSpec.proposeNext(currentSmcIteration, random, sampleOldLatent());
      nCalls++;
      return curLatent;
    }
    
    private boolean isInitial() 
    { 
      return previousPopulation == null; 
    }

    private P sampleOldLatent()
    {
      return previousPopulation.particles.get(permutationStream.popIndex());
    }

    @Override
    public int numberOfCalls()
    {
      return nCalls;
    }

    @Override
    public ProposalWithRestart<P> restart()
    {
      return new StreamingBootstrapProposal(seed, previousPopulation);
    }
  }

}
