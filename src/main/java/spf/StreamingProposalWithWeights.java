package spf;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import simplesmc.SMCProblemSpecification;
import bayonet.smc.ParticlePopulation;

/**
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 * @param <S>
 */
public class StreamingProposalWithWeights<S> implements ProposalWithRestart<S> 
{
    private final long seed;
    private final PermutationStream permutationStream;
    private final Random random;
    
    ParticlePopulation<S> previousPopulation;
    SMCProblemSpecification<S> problemSpec;
    private int nCalls = 0;

    public StreamingProposalWithWeights(long seed, SMCProblemSpecification<S> problemSpec, ParticlePopulation<S> previousPopulation)
    {
      this.seed = seed;
      this.random = new Random(seed * 171);
      this.problemSpec = problemSpec;
      this.previousPopulation = previousPopulation;
      this.permutationStream = previousPopulation == null ? null : new PermutationStream(previousPopulation.nParticles(), new Random(random.nextLong()));
    }

	@Override
	public Pair<Double, S> nextLogWeightSamplePair(int currentSmcIteration) {
		Random rand = new Random(getNextSeed());
		Pair<Double, S> curLatent = currentSmcIteration == 0 
	          ? problemSpec.proposeInitial(rand)
	          : problemSpec.proposeNext(currentSmcIteration, rand, sampleOldLatent());
	      return curLatent;
	}

    @Override
    public double nextLogWeight(int currentSmcIteration)
    {
		Random rand = new Random(getNextSeed());
        double curWeight = currentSmcIteration == 0 
                ? problemSpec.proposeInitialStream(rand)
                : problemSpec.proposeNextStream(currentSmcIteration, rand, sampleOldLatent());
        return curWeight;    	
    }

	@Override
	public int numberOfCalls() {
		return nCalls;
	}

	@Override
	public ProposalWithRestart<S> restart() {
		return new StreamingProposalWithWeights<S>(seed, problemSpec, previousPopulation);
	}
	
	private long getNextSeed()
	{
		nCalls++;
		return random.nextLong();
	}
	
	public void advanceStream(int currentSmcIteration)
	{
		if (currentSmcIteration > 0)
			sampleOldLatent();
		getNextSeed();
	}
	
    private S sampleOldLatent()
    {
      return previousPopulation.particles.get(permutationStream.popIndex());
    }


}
