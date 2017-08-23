package phylo;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class PhyloUtils 
{
	public static double [][] getTransitionMatrix(EvolutionaryModel model, double t)
	{
		DoubleMatrix Q = new DoubleMatrix(model.getRateMatrix());
		DoubleMatrix Qt = Q.mul(t);
		// exponentiate the rate matrix
		DoubleMatrix P = MatrixFunctions.expm(Qt);
		
		// check P is proper transition matrix
		if (!LikelihoodCalculatorExpFam.check(P))
			throw new RuntimeException("Not a propoer transition matrix");
		
		return P.toArray2();
	}

}
