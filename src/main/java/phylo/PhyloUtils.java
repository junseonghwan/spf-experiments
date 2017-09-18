package phylo;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import briefj.BriefIO;
import briefj.Indexer;
import pmcmc.proposals.RealVectorParameters;

public class PhyloUtils 
{
	public static double [][] getTransitionMatrix(EvolutionaryModel<RealVectorParameters> model, double t)
	{
		DoubleMatrix P = MatrixFunctions.expm(model.getRateMatrix().mul(t));		
		// check P is proper transition matrix
		if (!LikelihoodCalculatorExpFam.check(P))
			throw new RuntimeException("Not a propoer transition matrix");
		
		return P.toArray2();
	}

	public static List<Taxon> readData(String file)
	{
		List<Taxon> taxa = new ArrayList<>();
		Indexer<Taxon> taxonIndexer = new Indexer<>();
		for (String line : BriefIO.readLines(file))
		{
			String [] row = line.split("\\s+");
			if (row.length == 0)
				continue;
			Taxon taxon = new Taxon(row[0].trim());
			if (!taxonIndexer.containsObject(taxon)) {
				taxonIndexer.addToIndex(taxon);
				taxa.add(taxon);
			}
			for (int i = 1; i < row.length; i++)
				taxa.get(taxonIndexer.o2i(taxon)).appendSequence(row[i].trim());
		}
		
		return taxa;
	}
}
