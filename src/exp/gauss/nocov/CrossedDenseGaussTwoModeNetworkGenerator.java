package exp.gauss.nocov;

import umontreal.iro.lecuyer.rng.LFSR113;
import model.gauss.crossing.diffVar.NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkSampler;

public class CrossedDenseGaussTwoModeNetworkGenerator {

	public static void main(String[] args) {

		NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkSampler sampler = new NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkSampler();

		String modelParameters = args[0];
		String delim = args[1];
		if (Integer.parseInt(delim) == 0)
			delim = " ";
		else
			delim = "\t";
		String confirmedParameters = args[2];
		int numOfRows = Integer.parseInt(args[3]);
		int numOfColumns = Integer.parseInt(args[4]);
		String outputZFilename = args[5];
		String outputWFilename = args[6];
		String outputYFilename = args[7];
		int seed0 = Integer.parseInt(args[8].trim());
		int[] seeds = { 2, 8, 16, 128 };
		if (seed0 > 0)
			seeds[0] = -seed0;
		else
			seeds[0] = seed0;

		// Set model parameters
		System.out.println("Read model parameters");
		sampler.setParameters(modelParameters, delim);

		// Write model parameters to confirm the simulation setting
		System.out.println("Write model parameters");
		sampler.printModelParameters2File(confirmedParameters, delim);

		// Create a random generator
		LFSR113 randomGenerator = new LFSR113();
		randomGenerator.setSeed(seeds);

		// Generate a random network and save it to files
		System.out.println("Start generating the random network!");
		sampler.generateSampleWITHOUTCovariates(randomGenerator, numOfRows,
				numOfColumns, outputZFilename, outputWFilename,
				outputYFilename, delim);
		System.out.println("Finish generating the random network!");

	}
}
