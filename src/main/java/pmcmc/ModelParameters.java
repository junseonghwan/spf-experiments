package pmcmc;

public interface ModelParameters
{
	public void update(ModelParameters p);
	public void revert();
	public String asCommaSeparatedLine();
}
