package phylo;

import java.util.Random;

import org.apache.commons.math3.random.RandomGenerator;

public class CustomRandomGenerator implements RandomGenerator
{
	private Random random;
	
	public CustomRandomGenerator(long seed)
	{
		random = new Random(seed);
	}

	@Override
  public boolean nextBoolean() {
    return random.nextBoolean();
  }

	@Override
  public void nextBytes(byte[] bytes) {
		random.nextBytes(bytes);
  }

	@Override
  public double nextDouble() {
    return random.nextDouble();
  }

	@Override
  public float nextFloat() {
    return random.nextFloat();
  }

	@Override
  public double nextGaussian() {
    return random.nextGaussian();
  }

	@Override
  public int nextInt() {
    return random.nextInt();
  }

	@Override
  public int nextInt(int n) {
    return random.nextInt(n);
  }

	@Override
  public long nextLong() {
    return random.nextLong();
  }

	@Override
  public void setSeed(int seed) {
    random.setSeed(seed);
  }

	@Override
  public void setSeed(int[] arg0) {
  }

	@Override
  public void setSeed(long seed) {
		random.setSeed(seed);
  }
	
}
