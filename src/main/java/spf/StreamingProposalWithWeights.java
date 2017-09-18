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
	private long [] randomLongs = null;

    public StreamingProposalWithWeights(long seed, SMCProblemSpecification<S> problemSpec, ParticlePopulation<S> previousPopulation, int maxNumberOfVirtualParticles)
    {
      this.seed = seed;
      this.random = new Random(seed * 171);
      this.problemSpec = problemSpec;
      this.previousPopulation = previousPopulation;
      this.permutationStream = previousPopulation == null ? null : new PermutationStream(previousPopulation.nParticles(), new Random(random.nextLong()));
      if (randomLongs == null) {
          this.randomLongs = new long[maxNumberOfVirtualParticles];
          generateLongSeeds();
      }
    }
    
    public StreamingProposalWithWeights(long seed, SMCProblemSpecification<S> problemSpec, ParticlePopulation<S> previousPopulation, long [] randomSeeds)
    {
      this.seed = seed;
      this.random = new Random(seed * 171);
      this.problemSpec = problemSpec;
      this.previousPopulation = previousPopulation;
      this.permutationStream = previousPopulation == null ? null : new PermutationStream(previousPopulation.nParticles(), new Random(random.nextLong()));
      this.randomLongs = randomSeeds;
    }

    private void generateLongSeeds()
    {
    	Random rand = new Random(this.random.nextLong());
    	for (int i = 0; i < randomLongs.length; i++) {
    		randomLongs[i] = rand.nextLong();
    	}
    }

	@Override
	public Pair<Double, S> nextLogWeightSamplePair(int currentSmcIteration, int particleIdx) {
		Random rand = new Random(getSeed(particleIdx));
		Pair<Double, S> curLatent = currentSmcIteration == 0 
	          ? problemSpec.proposeInitial(rand)
	          : problemSpec.proposeNext(currentSmcIteration, rand, sampleOldLatent(rand, particleIdx));
	      return curLatent;
	}

    @Override
    public double nextLogWeight(int currentSmcIteration, int particleIdx)
    {
    	long seed = getSeed(particleIdx);
		Random rand = new Random(seed);
        double curWeight = currentSmcIteration == 0 
                ? problemSpec.proposeInitialStream(rand)
                : problemSpec.proposeNextStream(currentSmcIteration, rand, sampleOldLatent(rand, particleIdx));
        return curWeight;
    }

	@Override
	public int numberOfCalls() {
		return nCalls;
	}

	@Override
	public ProposalWithRestart<S> restart() {
		return new StreamingProposalWithWeights<S>(seed, problemSpec, previousPopulation, this.randomLongs);
	}

	private synchronized long getSeed(int particleIdx)
	{
		nCalls++;
		return randomLongs[particleIdx];
	}

	public void advanceStream(int currentSmcIteration, int particleIdx)
	{
		long seed = getSeed(particleIdx);
		Random rand = new Random(seed);
		if (currentSmcIteration > 0)
			sampleOldLatent(rand, particleIdx);
	}
	
    private S sampleOldLatent(Random rand, int particleIdx)
    {
      return previousPopulation.particles.get(permutationStream.popIndex(rand.nextInt(permutationStream.size())));    	
    }


}
