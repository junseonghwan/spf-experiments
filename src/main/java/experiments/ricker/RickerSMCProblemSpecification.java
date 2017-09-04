package experiments.ricker;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import dynamic.models.RickerModel;
import bayonet.distributions.Normal;
import bayonet.distributions.Poisson;
import bayonet.distributions.Uniform;
import simplesmc.SMCProblemSpecification;

public class RickerSMCProblemSpecification implements SMCProblemSpecification<Double> 
{
	private int T;
	private List<Integer> obs;
	private RickerModel model;

	public RickerSMCProblemSpecification(int T, RickerModel model, List<Integer> obs)
	{
		this.T = T;
		this.model = model;
		this.obs = obs;
	}

	@Override
	public Pair<Double, Double> proposeNext(int currentSmcIteration, Random random, Double currentParticle) 
	{
		// propose the error term from Normal
		double error = Normal.generate(random, 0.0, model.getParamValue("var"));
		double newParticle = model.getParamValue("r") * currentParticle * Math.exp(-currentParticle + error);

		// evaluate the likelihood
		double logWeight = Poisson.logDensity(obs.get(currentSmcIteration), model.getParamValue("phi")*newParticle);
		return Pair.of(logWeight, newParticle);
	}

	@Override
	public Pair<Double, Double> proposeInitial(Random random) 
	{
		// set uniform distribution on [0, 10) as the prior
		double newParticle = Uniform.generate(random, 0.0, 10.0);
		double logWeight = Poisson.logDensity(obs.get(0), model.getParamValue("phi")*newParticle);
		return Pair.of(logWeight, newParticle);
	}

	@Override
	public int nIterations() {
		return T + 1;
	}
	
}
