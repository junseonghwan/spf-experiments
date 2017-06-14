package experiments.kitagawa;

import java.util.Random;

import distributions.InverseGamma;
import bayonet.distributions.Normal;
import pmcmc.MCMCProblemSpecification;

public class KitagawaMCMCProblemSpecification implements MCMCProblemSpecification<KitagawaParams> 
{
	private final double shape;
	private final double rate;
	private final double sd_v2;
	private final double sd_w2;
	
	public KitagawaMCMCProblemSpecification(double shape, double rate, double sd_v, double sd_w) {
		this.shape = shape;
		this.rate = rate;
		this.sd_v2 = Math.pow(sd_v, 2.0);
		this.sd_w2 = Math.pow(sd_w, 2.0);
	}

	@Override
	public KitagawaParams propose(Random random, KitagawaParams curr) {
		// propose new values from Normal distribution with standard deviation provided in Section 3.1 of PMCMC
		KitagawaParams pstar = new KitagawaParams();
		double new_var_v = Math.pow(Normal.generate(random, Math.sqrt(curr.var_v.getValue()), sd_v2), 2);
		double new_var_w = Math.pow(Normal.generate(random, Math.sqrt(curr.var_w.getValue()), sd_w2), 2);
		pstar.var_v.setValue(new_var_v);
		pstar.var_w.setValue(new_var_w);
		return pstar;
	}

	@Override
	public double logProposalDensity(KitagawaParams curr, KitagawaParams prev) {
		return Normal.logDensity(curr.var_v.getValue(), prev.var_v.getValue(), sd_v2) + Normal.logDensity(curr.var_w.getValue(), prev.var_w.getValue(), sd_w2);
	}
	
	@Override
	public double logPriorDensity(KitagawaParams params) {
		return InverseGamma.logDensity(params.var_v.getValue(), rate, shape) + InverseGamma.logDensity(params.var_w.getValue(), rate, shape);
	}

	@Override
	public KitagawaParams initialize(Random random) {
		KitagawaParams params = new KitagawaParams();
		/*
		params.var_v.setValue(InverseGamma.generate(random, rate, shape));
		params.var_w.setValue(InverseGamma.generate(random, rate, shape));
		*/
		params.var_v.setValue(10.0);
		params.var_w.setValue(10.0);
		return params;
	}

}
