package distributions;

import java.util.Random;

import bayonet.distributions.Gamma;
import bayonet.distributions.UnivariateRealDistribution;
import bayonet.math.SpecialFunctions;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;

public class InverseGamma<P extends InverseGamma.Parameters> implements GenerativeFactor, UnivariateRealDistribution
{
	/**
	  * The variable on which this density is defined on.
	  */
	@FactorArgument(makeStochastic=true)
	public final RealVariable realization;

	/**
	 * The parameter of this exponential density.
	 */
	@FactorComponent public final P parameters;

	public static interface Parameters
	{
		public double getRate();
		public double getShape();
	}
	
	public static class ShapeRateParameterization implements Parameters
	{
	    @FactorArgument 
	    public final RealVariable rate;
	    @FactorArgument 
	    public final RealVariable shape;

		
		public ShapeRateParameterization(RealVariable shape, RealVariable rate) {
	    	if (shape.getValue() <= 0.0 || rate.getValue() <= 0.0)
	    		throw new RuntimeException("Shape and rate parameters must be positive!");

	    	this.shape = shape;
	    	this.rate = rate;
		}
		
		@Override
		public double getRate() {
			return rate.getValue();
		}

		@Override
		public double getShape() {
			return shape.getValue();
		}
		
	}

	public InverseGamma(RealVariable realization, P parameters)
	{
 		this.realization = realization;
 		this.parameters = parameters;
 	}

	@Override
	public double logDensity() {
		return logDensity(realization.getValue(), parameters.getRate(), parameters.getShape());		
	}
	
	public static double logDensity(double x, double rate, double shape) {
		if (x <= 0.0)
			return Double.NEGATIVE_INFINITY;
		
		return (shape * Math.log(rate) + (-shape - 1) * Math.log(x) - rate/x - SpecialFunctions.lnGamma(shape));
	}

	@Override
	public RealVariable getRealization() {
		return realization;
	}

	@Override
	public void generate(Random random) {
		realization.setValue(generate(random, parameters.getRate(), parameters.getShape()));
	}

	public static double generate(Random random, double rate, double shape) {
		double x = Gamma.generate(random, rate, shape);
		return 1/x;
	}
	
	public static InverseGamma<ShapeRateParameterization> on(RealVariable realization, RealVariable shape, RealVariable rate)
	{
		return new InverseGamma<ShapeRateParameterization>(realization, new ShapeRateParameterization(shape, rate));
	}

}
