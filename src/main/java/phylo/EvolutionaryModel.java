package phylo;

public interface EvolutionaryModel 
{
	public double [][] getRateMatrix();
	public double [] getStationaryDistribution();
}
