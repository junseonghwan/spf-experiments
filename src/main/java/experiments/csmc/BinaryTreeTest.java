package experiments.csmc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.MultivariateSummaryStatistics;

import bayonet.smc.ParticlePopulation;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import phylo.RootedPhylogeny;
import phylo.Taxon;
import simplesmc.SMCAlgorithm;
import simplesmc.SMCOptions;
import simplesmc.SMCProblemSpecification;
import util.OutputHelper;

public class BinaryTreeTest implements Runnable
{
	@Option public static final int numLeaves = 4;
	@OptionSet(name="smc") public static final SMCOptions options = new SMCOptions();

	public static Indexer<Taxon> indexer = new Indexer<>();
	public static List<Taxon> leaves = new ArrayList<>();
	public static UnorderedPair<Taxon, Taxon> t12;
	public static UnorderedPair<Taxon, Taxon> t34;
	public static UnorderedPair<Taxon, Taxon> t13;
	public static UnorderedPair<Taxon, Taxon> t14;

	static {
		for (int i = 0; i < numLeaves; i++) {
			leaves.add(new Taxon("T" + i));
			indexer.addToIndex(leaves.get(i));
		}

		t12 = UnorderedPair.of(leaves.get(0), leaves.get(1));
		t34 = UnorderedPair.of(leaves.get(2), leaves.get(3));
		t13 = UnorderedPair.of(leaves.get(0), leaves.get(2));
		t14 = UnorderedPair.of(leaves.get(0), leaves.get(3));
	}

	@Override
	public void run() 
	{
		PartialBinaryForestState initial = PartialBinaryForestState.getInitial(leaves);

		SMCProblemSpecification<PartialBinaryForestState> proposalWrong = new SMCProblemSpecification<BinaryTreeTest.PartialBinaryForestState>() {

			@Override
			public Pair<Double, PartialBinaryForestState> proposeNext(int currentSmcIteration, Random random,
					PartialBinaryForestState currentParticle) {
				PartialBinaryForestState newState = PartialBinaryForestState.copy(currentParticle);
				// select a pair at random
				RootedPhylogeny t1 = newState.trees.remove(random.nextInt(newState.trees.size()));
				RootedPhylogeny t2 = newState.trees.remove(random.nextInt(newState.trees.size()));
				double height = Math.max(t1.getHeight(), t2.getHeight()) + 1;
				double b1 = height - t1.getHeight();
				double b2 = height - t2.getHeight();
				RootedPhylogeny tm = new RootedPhylogeny(new Taxon("t" + currentSmcIteration), t1, t2, b1, b2, height, true);
				newState.trees.add(tm);
				return Pair.of(0.0, newState);
			}

			@Override
			public Pair<Double, PartialBinaryForestState> proposeInitial(Random random) {
				return proposeNext(0, random, initial);
			}

			@Override
			public int nIterations() {
				return numLeaves - 1;
			}
		};

		SMCProblemSpecification<PartialBinaryForestState> proposalCorrect = new SMCProblemSpecification<BinaryTreeTest.PartialBinaryForestState>() {

			@Override
			public Pair<Double, PartialBinaryForestState> proposeNext(int currentSmcIteration, Random random,
					PartialBinaryForestState currentParticle) {
				PartialBinaryForestState newState = PartialBinaryForestState.copy(currentParticle);
				// select a pair at random
				RootedPhylogeny t1 = newState.trees.remove(random.nextInt(newState.trees.size()));
				RootedPhylogeny t2 = newState.trees.remove(random.nextInt(newState.trees.size()));
				double height = Math.max(t1.getHeight(), t2.getHeight()) + 1;
				double b1 = height - t1.getHeight();
				double b2 = height - t2.getHeight();
				RootedPhylogeny tm = new RootedPhylogeny(new Taxon("t" + currentSmcIteration), t1, t2, b1, b2, height, true);
				newState.trees.add(tm);
				// count the number of non-trivial trees
				int numNonTrivialTrees = 0;
				for (RootedPhylogeny t : newState.trees)
				{
					if (!t.isLeaf())
						numNonTrivialTrees += 1;
				}
				double logw = -Math.log(numNonTrivialTrees);
				return Pair.of(logw, newState);
			}

			@Override
			public Pair<Double, PartialBinaryForestState> proposeInitial(Random random) {
				return proposeNext(0, random, initial);
			}

			@Override
			public int nIterations() {
				return numLeaves - 1;
			}
		};

		int reps = 5;
		int [] nParticles = new int[]{100, 200, 400, 800, 1600, 3200, 6400};
		double [][] estimatesIncorrect = new double[nParticles.length][5];
		double [][] estimatesCorrect = new double[nParticles.length][5];
		File resultDir = Results.getResultFolder();

		for (int i = 0; i < nParticles.length; i++)
		{
			options.nParticles = nParticles[i];
			options.essThreshold = 1.0;
			MultivariateSummaryStatistics statI = new MultivariateSummaryStatistics(2, true);
			MultivariateSummaryStatistics statC = new MultivariateSummaryStatistics(2, true);
			for (int rep  = 0; rep < reps; rep++) 
			{
				SMCAlgorithm<PartialBinaryForestState> smcCorrect = new SMCAlgorithm<>(proposalCorrect, options);
				ParticlePopulation<PartialBinaryForestState> samplesCorrect = smcCorrect.sample();

				SMCAlgorithm<PartialBinaryForestState> smcIncorrect = new SMCAlgorithm<>(proposalWrong, options);
				ParticlePopulation<PartialBinaryForestState> samplesIncorrect = smcIncorrect.sample();

				statI.addValue(count(samplesIncorrect.particles));
				statC.addValue(count(samplesCorrect.particles));
			}
			estimatesIncorrect[i][0] = nParticles[i];
			estimatesIncorrect[i][1] = statI.getMean()[0];
			estimatesIncorrect[i][2] = statI.getStandardDeviation()[0];
			estimatesIncorrect[i][3] = statI.getMean()[1];
			estimatesIncorrect[i][4] = statI.getStandardDeviation()[1];
			
			estimatesCorrect[i][0] = nParticles[i];
			estimatesCorrect[i][1] = statC.getMean()[0];
			estimatesCorrect[i][2] = statC.getStandardDeviation()[0];
			estimatesCorrect[i][3] = statC.getMean()[1];
			estimatesCorrect[i][4] = statC.getStandardDeviation()[1];
		}

		OutputHelper.writeDoubleArray2DAsCSV(new File(resultDir, "estimates-correct.csv"), 
				new String[]{"nParticles", "mu_balanced", "sd_balanced", "mu_unbalanced", "sd_unbalanced"}, estimatesCorrect);
		OutputHelper.writeDoubleArray2DAsCSV(new File(resultDir, "estimates-incorrect.csv"),
				new String[]{"nParticles", "mu_balanced", "sd_balanced", "mu_unbalanced", "sd_unbalanced"}, estimatesIncorrect);
	}
	
	private static double [] count(List<PartialBinaryForestState> states)
	{
		double [] estimates = new double[2];
		int N = states.size();
		for (PartialBinaryForestState state : states)
		{
			Counter<UnorderedPair<Taxon, Taxon>> counter = state.trees.get(0).getPairwiseDistances(indexer);
			// check for the following states 
			// 1. where the distance between (T1, T2) = 2 and (T3, T4) = 2
			// 2. where the distance between (T1, T2) = 2 and (T1, T3) = 4
			if (counter.getCount(t12) == 2 && counter.getCount(t34) == 2)
				estimates[0] += 1;
			if (counter.getCount(t12) == 2 && counter.getCount(t13) == 4 && counter.getCount(t14) == 6)
				estimates[1] += 1;
		}

		estimates[0] /= N;
		estimates[1] /= N;
		System.out.println(estimates[0] + ", " + estimates[1]);
		return estimates;
	}

	public static void main(String [] args)
	{
		Mains.instrumentedRun(args, new BinaryTreeTest());
	}
	
	public static class PartialBinaryForestState
	{
		List<RootedPhylogeny> trees;
		public static PartialBinaryForestState getInitial(List<Taxon> leaves)
		{
			PartialBinaryForestState state = new PartialBinaryForestState();
			state.trees = new ArrayList<>();
			for (Taxon taxon : leaves)
			{
				state.trees.add(new RootedPhylogeny(taxon));				
			}
			return state;
		}
		
		public static PartialBinaryForestState copy(PartialBinaryForestState src)
		{
			PartialBinaryForestState newState = new PartialBinaryForestState();
			newState.trees = new ArrayList<>(src.trees);
			return newState;
		}
	}
	

}
