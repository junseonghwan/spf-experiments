package experiments.phylo;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import phylo.PartialCoalescentState;
import phylo.Taxon;
import simplesmc.SMCProblemSpecification;

public class PriorPriorProblemSpecification implements SMCProblemSpecification<PartialCoalescentState> {

	private PartialCoalescentState initialState;
	
	public PriorPriorProblemSpecification(List<Taxon> taxa)
	{
		this.initialState = PartialCoalescentState.getInitial(taxa);
	}

	@Override
	public Pair<Double, PartialCoalescentState> proposeNext(int currentSmcIteration, Random random, PartialCoalescentState currentParticle) {
		PartialCoalescentState newState = currentParticle.coalesce(random);
		double logLik = newState.logLikelihood();
		return Pair.of(logLik, newState);
	}

	@Override
	public Pair<Double, PartialCoalescentState> proposeInitial(Random random) {
		return proposeNext(0, random, initialState);
	}

	@Override
	public int nIterations() {
		return initialState.numTrees() - 1;
	}

}
