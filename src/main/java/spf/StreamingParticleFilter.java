package spf;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.smc.ParticlePopulation;
import phylo.FelsensteinPruningSystBiol2012;
import phylo.PhyloOptions;
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
  private List<Double> ess;
  private List<Double> logZs;
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
    System.out.println("params: " + ((FelsensteinPruningSystBiol2012)PhyloOptions.calc).getModel().getModelParameters().asCommaSeparatedLine());

    // instantiate new arraylist each time
    nImplicitParticles = new ArrayList<>(nSMCIterations);
    relESS = new ArrayList<>();
    ess = new ArrayList<>();
    logZs = new ArrayList<>();
    timeInSeconds = new ArrayList<>(nSMCIterations);
    long start = 0, end = 0;

    if (options.verbose)
    	System.out.println("Iter=" + 1);
    // initial distribution
    start = System.currentTimeMillis();
    ProposalWithRestart<P> proposal = getInitialDistributionProposal();
    StreamingPropagator<P> propagator = new StreamingPropagator<P>(proposal, options);
    PropagationResult<P> propResults = propagator.execute(0);
    end = System.currentTimeMillis();
    logZ = propResults.population.logZEstimate();
    logZs.add(logZ);
    nImplicitParticles.add((double)propResults.population.getNumberOfParticles());
    relESS.add(propResults.population.ess()/options.numberOfConcreteParticles);
    ess.add(propResults.population.ess());
    timeInSeconds.add((end - start)/1000.);
    if (processor != null)
  	  processor.process(0, propResults.getParticlePopulation());

    // recursion
    for (int i = 1; i < nSMCIterations; i++)
    {

      start = System.currentTimeMillis();
      if (options.storeParticleWeights)
    	  proposal = new StreamingProposalWithWeights<>(options.mainRandom.nextLong(), problemSpec, propResults.getParticlePopulation(), options.maxNumberOfVirtualParticles);
      else
    	  proposal = new StreamingBootstrapProposal(options.mainRandom.nextLong(), propResults.getParticlePopulation());
      propagator = new StreamingPropagator<>(proposal, options);
      propResults = propagator.execute(i);
      end = System.currentTimeMillis();
      double logZr = propResults.population.logZEstimate();
      logZs.add(logZr);
      logZ += logZr;
      nImplicitParticles.add((double)propResults.population.getNumberOfParticles());
      relESS.add(propResults.population.ess()/options.numberOfConcreteParticles);
      ess.add(propResults.population.ess());
      timeInSeconds.add((end - start)/1000.);
      if (processor != null)
    	  processor.process(i, propResults.getParticlePopulation());
    }

    return propResults.getParticlePopulation();
  }
  
  public List<Double> getLogZs() { return logZs; }
  public double logNormEstimate() { return logZ; }
  public List<Double> nImplicitParticles() { return nImplicitParticles; }
  public List<Double> relESS() { return relESS; }
  public List<Double> ess() { return ess; }
  public List<Double> getExecutionTimes() { return timeInSeconds; }
  
  private ProposalWithRestart<P> getInitialDistributionProposal()
  {
	  if (options.storeParticleWeights)
		  return new StreamingProposalWithWeights<>(options.mainRandom.nextLong(), problemSpec, null, options.maxNumberOfVirtualParticles);
	  else
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
    public Pair<Double, P> nextLogWeightSamplePair(int currentSmcIteration, int particleIdx)
    {
      // terminology: old means the SMC generation before current (null if we are doing initial)
      //              cur means the current SMC generation
      Pair<Double, P> curLatent = isInitial() 
          ? problemSpec.proposeInitial(random)
          : problemSpec.proposeNext(currentSmcIteration, random, sampleOldLatent());
      nCalls++;
      return curLatent;
    }
    
    @Override
    public double nextLogWeight(int currentSmcIteration, int particleIdx)
    {
        double curWeight = isInitial() 
                ? problemSpec.proposeInitialStream(random)
                : problemSpec.proposeNextStream(currentSmcIteration, random, sampleOldLatent());
            nCalls++;
        return curWeight;    	
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
