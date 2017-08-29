package phylo;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

public class PhyloUtils 
{
	public static double [][] getTransitionMatrix(EvolutionaryModel model, double t)
	{
		DoubleMatrix P = MatrixFunctions.expm(model.getRateMatrix().mul(t));		
		// check P is proper transition matrix
		if (!LikelihoodCalculatorExpFam.check(P))
			throw new RuntimeException("Not a propoer transition matrix");
		
		return P.toArray2();
	}

}
