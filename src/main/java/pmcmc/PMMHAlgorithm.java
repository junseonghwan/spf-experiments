package pmcmc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import briefj.run.Results;
import simplesmc.AbstractSMCAlgorithm;
import simplesmc.SMCAlgorithm;
import spf.StreamingParticleFilter;
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
	private final MCMCProposal<P> mcmcProposal;
	private final ProbabilityDistribution<P> prior;
	private final PMCMCOptions options;
	private final Model<P> model;
	private final List<PMCMCProcessor<P>> processors;
	private int nAccepts = 0;
	private SummaryStatistics marginalLikelihoodStat = new SummaryStatistics();
	private LogZProcessor<P> logZProcessor;
	private List<SummaryStatistics> smcStatistics = null;
	private List<SummaryStatistics> smcTimingStatistics = null;
	private List<P> paramsDuringBurnIn;
	private boolean adaptiveMCMC = false;

	public PMMHAlgorithm(
			Model<P> model,
			AbstractSMCAlgorithm<S> smcAlgorithm, 
			MCMCProposal<P> mcmcProblemSpecification,
			ProbabilityDistribution<P> prior,
			PMCMCOptions options,
			boolean adaptiveMCMC)
	{
		this(model, smcAlgorithm, mcmcProblemSpecification, prior, options, null, null, adaptiveMCMC);
	}

	public PMMHAlgorithm(
			Model<P> model,
			AbstractSMCAlgorithm<S> smcAlgorithm, 
			MCMCProposal<P> mcmcProblemSpecification,
			ProbabilityDistribution<P> prior,
			PMCMCOptions options,
			List<PMCMCProcessor<P>> processors,
			LogZProcessor<P> logZProcessor,
			boolean adaptiveMCMC)
	{
		this.model = model;
		this.smcAlgorithm = smcAlgorithm;
		this.mcmcProposal = mcmcProblemSpecification;
		this.prior = prior;
		this.options = options;
		this.logZProcessor = logZProcessor;
		this.adaptiveMCMC = adaptiveMCMC;
		if (adaptiveMCMC)
			this.paramsDuringBurnIn = new ArrayList<>();			
		
		if (processors == null)
			this.processors = new ArrayList<>();
		else
			this.processors = processors;
		
		// add the default output processor
		if (smcAlgorithm instanceof SMCAlgorithm)
			this.processors.add(new PMCMCDefaultOutputProcessor<>("smc"));
		else if (smcAlgorithm instanceof StreamingParticleFilter)
			this.processors.add(new PMCMCDefaultOutputProcessor<>("spf"));
	}

	public void sample()
	{
		// run SMC algorithm for the initial parameter p to obtain 
		smcAlgorithm.sample();
		double logZCurr = smcAlgorithm.logNormEstimate();
		double logPriorCurr = prior.logDensity(model.getModelParameters());

		for (int iter = 0; iter < options.nIter; iter++)
		{
			// propose new values for the parameters
			P pcurr = model.getModelParameters();
			P pstar = mcmcProposal.propose(options.random, model.getModelParameters());
			// quick check to see if pstar is in the support set:
			if (Double.isInfinite(prior.logDensity(pstar)))
					continue;
			// compute the acceptance ratio
			double q = mcmcProposal.logProposalDensity(pcurr, pstar);
			double qstar = mcmcProposal.logProposalDensity(pstar, pcurr);
			if (Double.isInfinite(qstar))
				continue;

			// update the params with the newly proposed pstar and run SMC algorithm to get an estimate of the marginal likelihood
			model.updateModelParameters(pstar);
			smcAlgorithm.sample();
			updateSMCStatistics(smcAlgorithm);

			double logZStar = smcAlgorithm.logNormEstimate();
			double logPriorStar = prior.logDensity(pstar);

			double a = logZStar + logPriorStar + q;
			a -= (logZCurr + logPriorCurr + qstar);
			double acceptanceProb = Math.min(1.0, Math.exp(a));
			double u = options.random.nextDouble();
			System.out.println("PMCMC Iter " + iter + ": ");
			System.out.println("p=" + pcurr.asCommaSeparatedLine() + ", " +
					"p*=" + pstar.asCommaSeparatedLine() + ", " + 
					"logZ=" + logZCurr + ", " + 
					"logZ*=" + logZStar + ", " + 
					"q=" + q + ", " +
					"q*=" + qstar + ", " +
					"logPrior=" + logPriorCurr + ", " + 
					"logPrior*=" + logPriorStar + ", " +
					"a=" + a + ", " + ", " + 
					"P(accept)=" + acceptanceProb + ", " +
					"u= " + u);
			if (u < acceptanceProb) {
				// accept the proposal
				logZCurr = logZStar;
				logPriorCurr = logPriorStar;
				nAccepts++;
			} else {
				// revert the params
				model.revert();
			}
			
			if (iter < options.burnIn && adaptiveMCMC) {
				paramsDuringBurnIn.add(model.getModelParameters());
			} else if (iter == options.burnIn && adaptiveMCMC) {
				mcmcProposal.adapt(paramsDuringBurnIn);
			}
			if (processors != null && iter >= options.burnIn && iter % options.thinningPeriod == 0) {
				for (PMCMCProcessor<P> processor : processors)
					processor.process(pstar); 
			}

			// collect the marginal likelihood stat for each proposed value of pstar
			marginalLikelihoodStat.addValue(logZStar);
			if (logZProcessor != null)
				logZProcessor.process(iter, pstar, logZStar);
		}

		// output the parameters 
		File results = Results.getResultFolder();
		for (PMCMCProcessor<P> processor : processors)
			processor.output(new File(results, processor.outputPrefix() + "." + processor.getClass().getName() + ".csv"));

		// output statistics
		List<String> contents = new ArrayList<>();
		//contents.add(marginalLikelihoodStat.toString());
		contents.add("nAccepts: " + nAccepts);
		contents.add("nIter: " + options.nIter);
		contents.add("burnIn: " + options.burnIn);
		OutputHelper.writeLines(new File(results, "summary.txt"), contents);

		// output the marginal likelihoods
		if (logZProcessor != null)
			logZProcessor.output(new File(results, logZProcessor.getOutputPrefix() + ".logZ.csv"));

		// output smc statistics for each iteration: timing results,
		List<String> smcOutputLines = new ArrayList<>();
		smcOutputLines.add("Iter, Avg, Var, TimeAvg, TimeVar");
		for (int i = 0; i < smcStatistics.size(); i++)
		{
			SummaryStatistics stat = smcStatistics.get(i);
			SummaryStatistics timingStat = smcTimingStatistics.get(i);

			// form the output line
			smcOutputLines.add(i + ", " + stat.getMean() + ", " + stat.getVariance() + ", " + timingStat.getMean() + ", " + timingStat.getVariance());
		}
		OutputHelper.writeLines(new File(results, "smcStat.csv"), smcOutputLines);
	}

	public int nAccepts() { return nAccepts; }

	private void updateSMCStatistics(AbstractSMCAlgorithm<S> abstractSMC)
	{
		List<Double> stat;

		if (smcAlgorithm instanceof SMCAlgorithm)
		{
			SMCAlgorithm<S> smc = (SMCAlgorithm<S>)abstractSMC;
			stat = smc.effectiveSampleSize();
		}
		else if (smcAlgorithm instanceof StreamingParticleFilter)
		{
			StreamingParticleFilter<S> spf = (StreamingParticleFilter<S>)abstractSMC;
			stat = spf.nImplicitParticles();			
		}
		else
			throw new RuntimeException("Only the standard SMC and SPF algorithms are supported.");
		
		// process stat
		List<Double> timingInSeconds = abstractSMC.timeInSeconds();
		if (stat.size() != timingInSeconds.size())
			throw new RuntimeException("SMC statistics collection and timing results collection are incompatible! Bug in the program.");
		
		if (smcStatistics == null || smcTimingStatistics == null)
		{
			// initialize the stat collection storage
			smcStatistics = new ArrayList<>(stat.size());
			smcTimingStatistics = new ArrayList<>(stat.size());
			for (int i = 0; i < stat.size(); i++)
			{
				smcStatistics.add(new SummaryStatistics());
				smcTimingStatistics.add(new SummaryStatistics());
			}
		}
		
		for (int i = 0; i < stat.size(); i++)
		{
			smcStatistics.get(i).addValue(stat.get(i));
			smcTimingStatistics.get(i).addValue(timingInSeconds.get(i));
		}
	}
}
