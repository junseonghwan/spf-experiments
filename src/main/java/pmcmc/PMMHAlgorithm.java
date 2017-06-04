package pmcmc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import briefj.run.Results;
import simplesmc.AbstractSMCAlgorithm;
import util.OutputHelper;

/**
 * Framework for inference using PMMH. 
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class PMMHAlgorithm<P extends ModelParameters, S> 
{
	private final AbstractSMCAlgorithm<S> smcAlgorithm;
	private final MCMCProblemSpecification<P> mcmcProblemSpecification;
	private final PMCMCOptions options;
	private final P params;
	private final List<PMCMCProcessor<P>> processors;
	private int nAccepts = 0;
	private SummaryStatistics marginalLikelihoodStat = new SummaryStatistics();
	private LogZProcessor<P> logZProcessor;

	public PMMHAlgorithm(
			P params,
			AbstractSMCAlgorithm<S> smcAlgorithm, 
			MCMCProblemSpecification<P> mcmcProblemSpecification,
			PMCMCOptions options)
	{
		this(params, smcAlgorithm, mcmcProblemSpecification, options, null);
	}

	public PMMHAlgorithm(
			P params,
			AbstractSMCAlgorithm<S> smcAlgorithm, 
			MCMCProblemSpecification<P> mcmcProblemSpecification,
			PMCMCOptions options,
			List<PMCMCProcessor<P>> processors)
	{
		this.smcAlgorithm = smcAlgorithm;
		this.mcmcProblemSpecification = mcmcProblemSpecification;
		this.options = options;
		this.params = params;
		this.processors = processors;
		
		if (options.collectMarginalLikelihoods) {
			logZProcessor = new LogZProcessor<>();
		}
	}

	public void sample()
	{
		// run SMC algorithm for the given parameter p to obtain 
		smcAlgorithm.sample();
		double logZCurr = smcAlgorithm.logNormEstimate();
		double logPriorCurr = mcmcProblemSpecification.logPriorDensity(params);
		
		for (int iter = 0; iter < options.nIter; iter++)
		{
			// propose new values for the parameters
			P pstar = mcmcProblemSpecification.propose(options.random, params);
			// update the params with the newly proposed pstar
			params.update(pstar);
			// run SMC algorithm to get an estimate of the marginal likelihood
			smcAlgorithm.sample();
			double logZStar = smcAlgorithm.logNormEstimate();
			double logPriorStar = mcmcProblemSpecification.logPriorDensity(pstar);

			// compute the acceptance ratio
			double a = logZStar + logPriorStar + mcmcProblemSpecification.logProposalDensity(params, pstar);
			a -= (logZCurr + logPriorCurr + mcmcProblemSpecification.logProposalDensity(pstar, params));
			a = Math.exp(a);
			double u = options.random.nextDouble();
			if (u < a) {
				// accept the proposal
				logZCurr = logZStar;
				logPriorCurr = logPriorStar;
				nAccepts++;
			} else {
				// revert the params
				params.revert();
			}
			
			if (iter >= options.burnIn && iter % options.thinningPeriod == 0) {
				System.out.println(iter + ", " + params + ", " + logZCurr);
				System.out.println(iter + ", " + pstar + ", " + logZStar);
				for (PMCMCProcessor<P> processor : processors)
					processor.process(pstar); 
			}

			// collect the marginal likelihood stat for each proposed value of pstar
			marginalLikelihoodStat.addValue(logZStar);
			if (options.collectMarginalLikelihoods)
				logZProcessor.process(pstar, logZStar);

		}

		// output the parameters
		File results = Results.getResultFolder();
		for (PMCMCProcessor<P> processor : processors)
			processor.output(new File(results, processor.getClass().getName() + ".csv"));

		// output statistics
		List<String> contents = new ArrayList<>();
		contents.add(marginalLikelihoodStat.toString());
		contents.add("nAccepts: " + nAccepts);
		contents.add("nIter: " + options.nIter);
		OutputHelper.writeLines(new File(results, "summary.txt"), contents);
		
		// output the marginal likelihoods
		if (options.collectMarginalLikelihoods)
			logZProcessor.output(new File(results, "logZ.csv"));
	}
	
	public int nAccepts() { return nAccepts; }

}
