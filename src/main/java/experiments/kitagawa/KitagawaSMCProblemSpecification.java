package experiments.kitagawa;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.KitagawaModel;
import bayonet.distributions.Normal;
import simplesmc.SMCProblemSpecification;

public class KitagawaSMCProblemSpecification implements SMCProblemSpecification<Double>
{
	private final KitagawaModel model;
	private final double [] y;
	
	public KitagawaSMCProblemSpecification(KitagawaModel model, double [] y) 
	{
		this.model = model;
		this.y = y;
	}

	@Override
	public Pair<Double, Double> proposeNext(int currentSmcIteration, Random random, Double currentParticle) {
		double mu = currentParticle/2.0 + 25*currentParticle / (1 + Math.pow(currentParticle, 2.0)) + 8*Math.cos(1.2*currentSmcIteration);
		double xstar = Normal.generate(random, mu, model.getVarV());
		// evaluate the weight: w(x_r) = p(y_r|x_r)
		double logw = Normal.logDensity(y[currentSmcIteration], Math.pow(xstar,  2.0)/20.0, model.getVarW());
		return Pair.of(logw, xstar);
	}

	@Override
	public Pair<Double, Double> proposeInitial(Random random) {
		double xstar = Normal.generate(random, 0.0, 5.0);
		// evaluate the weight: w(x_1) = p(y_1|x_1)
		double logw = Normal.logDensity(y[0], Math.pow(xstar,  2.0)/20.0, model.getVarW());
		return Pair.of(logw, xstar);
	}
	
	@Override
	public int nIterations() {
		return y.length;
	}

}
