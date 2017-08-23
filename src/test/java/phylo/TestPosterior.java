package phylo;

import org.junit.Test;

public class TestPosterior 
{
	@Test
	/**
	 * Follow the validation framework of Cook, Gelman, and Rubin (2006).
	 * 1. Generate the tree and then the data.
	 * 2. Run SMC to draw samples from the posterior distribution given the data.
	 * 3. Compute the quantile of the true height amongst the heights from trees sampled using SMC.
	 * 4. Repeat the above for Nrep times.
	 * 5. Compute the chisq test statistic and then, compute the p-value. Extreme p-value indicates error in the posterior sampling code. 
	 */
	public void validatePosteriorCode()
	{
	}


}
