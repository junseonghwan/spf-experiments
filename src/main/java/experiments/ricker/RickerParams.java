package experiments.ricker;

import blang.variables.RealVariable;
import pmcmc.ModelParameters;

public class RickerParams implements ModelParameters
{
	public final RealVariable phi = new RealVariable(0.0);
	public final RealVariable r = new RealVariable(0.0);
	public final RealVariable var = new RealVariable(0.0);
	public RickerParams old = null;
	
	public RickerParams(double phi, double r, double var) { 
		this.phi.setValue(phi);
		this.r.setValue(r);
		this.var.setValue(var);
	}
	
	public void saveCurrent()
	{
		old = new RickerParams(this.phi.getValue(), this.r.getValue(), this.var.getValue());
	}

	@Override
	public void update(ModelParameters p) 
	{
		saveCurrent();
		RickerParams newParams = (RickerParams)p;
		this.phi.setValue(newParams.phi.getValue());
		this.r.setValue(newParams.r.getValue());
		this.var.setValue(newParams.var.getValue());
	}

	@Override
	public void revert() 
	{
		if (old == null)
			throw new RuntimeException("Nothing to revert!");

		this.phi.setValue(old.phi.getValue());
		this.r.setValue(old.r.getValue());
		this.var.setValue(old.var.getValue());
	}

	@Override
	public String asCommaSeparatedLine() 
	{
		return toString();
	}
	
	@Override 
	public String toString()
	{
		return phi.getValue() + ", " + r.getValue() + ", " + var.getValue();		
	}

}
