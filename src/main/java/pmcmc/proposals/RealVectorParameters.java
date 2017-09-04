package pmcmc.proposals;

import pmcmc.ModelParameters;

public class RealVectorParameters implements ModelParameters {

	private int dim;
	private double [] vec;
	public RealVectorParameters(int dim) {
		this.dim = dim;
	}

	public RealVectorParameters(double [] vec) {
		this.vec = vec;
		this.dim = vec.length;
	}

	@Override
	public String asCommaSeparatedLine() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < dim; i++) {
			sb.append(vec[i]);
			if (i < dim - 1)
				sb.append(", ");
		}
		return sb.toString();
	}
	
	public double [] getVector() {
		return vec;
	}
	
	public int getDim() { return dim; }

}
