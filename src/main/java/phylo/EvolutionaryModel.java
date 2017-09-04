package phylo;

import org.jblas.DoubleMatrix;

import pmcmc.Model;
import pmcmc.ModelParameters;

public interface EvolutionaryModel<P extends ModelParameters> extends Model<P>
{
	public DoubleMatrix getRateMatrix();
	public double [] getStationaryDistribution();
}
