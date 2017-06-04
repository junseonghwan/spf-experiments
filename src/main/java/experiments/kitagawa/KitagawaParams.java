package experiments.kitagawa;

import pmcmc.ModelParameters;
import blang.variables.RealVariable;
import simplesmc.pmcmc.WithSignature;

public class KitagawaParams implements WithSignature, ModelParameters
{
	public final RealVariable var_v = new RealVariable(1.0);
	public final RealVariable var_w = new RealVariable(1.0);
	public KitagawaParams old;

	private KitagawaParams copy() {
		KitagawaParams copy = new KitagawaParams();
		copy.var_v.setValue(this.var_v.getValue());
		copy.var_w.setValue(this.var_w.getValue());
		return copy;
	}

	@Override
	public long signature() {
		long b = Double.hashCode(var_v.getValue());
		long c = Double.hashCode(var_w.getValue()) << 32;
		return b + c;
	}

	@Override
	public void update(ModelParameters p) {
		this.old = copy();
		KitagawaParams params = (KitagawaParams)p;
		this.var_v.setValue(params.var_v.getValue());
		this.var_w.setValue(params.var_w.getValue());
	}

	@Override
	public void revert() {
		this.var_v.setValue(old.var_v.getValue());
		this.var_w.setValue(old.var_w.getValue());
		this.old = null;
	}
	
	@Override
	public String toString()
	{
		return var_v.getValue() + ", " + var_w.getValue();
	}

	@Override
	public String asCommaSeparatedLine() {
		return var_v.getValue() + ", " + var_w.getValue();
	}
}
