package models;

import java.util.Random;

import bayonet.distributions.Normal;

public class OrnsteinUhlenbeckProcess 
{
	
	/**
	 * Simulate from dY_t = (\theta1 - \theta2 Y_t) dt + \\theta3 d W_t using Euler-Maruyama
	 * 
	 * @param random
	 * @param tbegin
	 * @param tend
	 * @param numPoints
	 * @param y0
	 * @param theta1
	 * @param theta2
	 * @param theta3
	 */
	public static double [] simulate(Random random, double tbegin, double tend, int numPoints, double y0, double theta1, double theta2, double theta3)
	{
		double dt = (tend - tbegin)/numPoints;
		double [] y = new double[numPoints];
		y[0] = y0;
		
		for (int i = 1; i < numPoints; i++)
		{
			double z = Normal.generate(random, 0.0, dt);
			y[i] = y[i-1] + (theta1 + theta2 * y[i-1]) + theta3 * z;
		}
		return y;
	}
	
	public static void main(String [] args)
	{
		Random random = new Random(1);
		double [] y = OrnsteinUhlenbeckProcess.simulate(random, 0, 10, 10, 0, 0.2, 0.3, 0.4);
		for (double val : y)
			System.out.println(val);
	}

}
