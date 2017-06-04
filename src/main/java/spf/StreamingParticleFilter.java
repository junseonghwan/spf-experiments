package spf;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.smc.ParticlePopulation;
import simplesmc.AbstractSMCAlgorithm;
import simplesmc.SMCProblemSpecification;
import spf.StreamingPropagator.PropagationResult;

public class StreamingParticleFilter<P> extends AbstractSMCAlgorithm<P>
{
  public final SMCProblemSpecification<P> problemSpec;
  private final SPFOptions options;
  private final Random mainRandom = new Random(1);
  private double logZ = 0.0;
  private List<Integer> nImplicitParticles;

  public StreamingParticleFilter(SMCProblemSpecification<P> problemSpec, SPFOptions options)
  {
    this.problemSpec = problemSpec;
    this.options = options;
  }

  public ParticlePopulation<P> sample()
  {
    int nSMCIterations = problemSpec.nIterations();
    nImplicitParticles = new ArrayList<>(nSMCIterations);

    // initial distribution
    StreamingBootstrapProposal proposal = getInitialDistributionProposal();
    StreamingPropagator<P> propagator = new StreamingPropagator<P>(proposal, options);
    PropagationResult<P> propResults = propagator.execute(0);
    logZ = propResults.population.logZEstimate();
    nImplicitParticles.add(propResults.population.getNumberOfParticles());

    // recursion
    for (int i = 1; i < nSMCIterations; i++)
    {
      proposal = new StreamingBootstrapProposal(mainRandom.nextLong(), propResults.getParticlePopulation());
      propagator = new StreamingPropagator<>(proposal, options);
      propResults = propagator.execute(i);
      logZ += propResults.population.logZEstimate();
      nImplicitParticles.add(propResults.population.getNumberOfParticles());
    }

    return propResults.getParticlePopulation();
  }

  public double logNormEstimate() { return logZ; }
  public List<Integer> nImplicitParticles() { return nImplicitParticles; }
  
  private StreamingBootstrapProposal getInitialDistributionProposal()
  {
    return new StreamingBootstrapProposal(mainRandom.nextLong(), null);
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
