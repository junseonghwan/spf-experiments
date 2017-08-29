package phylo;

import org.jblas.DoubleMatrix;

public interface EvolutionaryModel 
{
	public DoubleMatrix getRateMatrix();
	public double [] getStationaryDistribution();
}
