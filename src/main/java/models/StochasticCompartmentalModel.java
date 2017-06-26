package models;

import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

/**
 * Implements the compartmental model for susceptible (S), infected (I), recovered (R) model with k recovered compartments
 * Reference: https://kingaa.github.io/sbied/stochsim/stochsim.html
 * 
 * dN_{SI}/dt = rateSI(t) * S(t)
 * dN_{IR}/dt = rateIR * I(t)
 * 
 * @author Seong-Hwan Jun (s2jun.uw@gmail.com)
 *
 */
public class StochasticCompartmentalModel 
{
	public static Pair<double [], double []> simulate(Random random,
			int S, int I, int R,
			double birthRate, double [] deathRates, 
			double contactRate, double rateIR)
	{
		return null;
	}

}
