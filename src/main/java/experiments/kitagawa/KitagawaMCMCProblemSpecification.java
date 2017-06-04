package experiments.kitagawa;

import java.util.Random;

import bayonet.distributions.InverseGamma;
import bayonet.distributions.Normal;
import pmcmc.MCMCProblemSpecification;

public class KitagawaMCMCProblemSpecification implements MCMCProblemSpecification<KitagawaParams> 
{
	private final double shape;
	private final double rate;
	private final double sd_v;
	private final double sd_w;
	
	public KitagawaMCMCProblemSpecification(double shape, double rate, double sd_v, double sd_w) {
		this.shape = shape;
		this.rate = rate;
		this.sd_v = sd_v;
		this.sd_w = sd_w;
	}

	@Override
	public KitagawaParams propose(Random random, KitagawaParams curr) {
		// propose new values from Normal distribution with standard deviation provided in Section 3.1 of PMCMC
		KitagawaParams pstar = new KitagawaParams();
		pstar.var_v.setValue(Normal.generate(random, curr.var_v.getValue(), sd_v));
		pstar.var_w.setValue(Normal.generate(random, curr.var_w.getValue(), sd_w));
		return pstar;
	}

	@Override
	public double logProposalDensity(KitagawaParams curr, KitagawaParams prev) {
		return Normal.logDensity(curr.var_v.getValue(), prev.var_v.getValue(), sd_v) + Normal.logDensity(curr.var_w.getValue(), prev.var_w.getValue(), sd_w);
	}

	@Override
	public double logPriorDensity(KitagawaParams params) {
		return InverseGamma.logDensity(params.var_v.getValue(), rate, shape) + InverseGamma.logDensity(params.var_w.getValue(), rate, shape);
	}

	@Override
	public KitagawaParams initialize(Random random) {
		KitagawaParams params = new KitagawaParams();
		params.var_v.setValue(InverseGamma.generate(random, rate, shape));
		params.var_w.setValue(InverseGamma.generate(random, rate, shape));
		return params;
	}

}
