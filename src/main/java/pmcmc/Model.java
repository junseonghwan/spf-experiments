package pmcmc;

public interface Model<P extends ModelParameters> 
{
	public P getModelParameters();
	public void updateModelParameters(P p);
	public void revert();
}
