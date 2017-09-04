package experiments.phylo;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import phylo.PartialCoalescentState;
import phylo.Taxon;
import simplesmc.SMCProblemSpecification;

public class PriorPostProblemSpecification implements SMCProblemSpecification<PartialCoalescentState> {

	private PartialCoalescentState initialState;
	
	public PriorPostProblemSpecification(List<Taxon> taxa)
	{
		this.initialState = PartialCoalescentState.getInitial(taxa);
	}

	@Override
	public Pair<Double, PartialCoalescentState> proposeNext(
	    int currentSmcIteration, Random random,
	    PartialCoalescentState currentParticle) {
		
		return currentParticle.priorPost(random, false);
	}
	
	@Override
  public double proposeNextStream(int currentSmcIteration, Random random, PartialCoalescentState currentParticle)
  {
		return currentParticle.priorPost(random, true).getLeft();
  }

	@Override
	public Pair<Double, PartialCoalescentState> proposeInitial(Random random) {
		return proposeNext(0, random, initialState);
	}

	@Override
	public double  proposeInitialStream(Random random)
  {
		return proposeNextStream(0, random, initialState);
  }


	@Override
	public int nIterations() {
		return initialState.numTrees() - 1;
	}

}
