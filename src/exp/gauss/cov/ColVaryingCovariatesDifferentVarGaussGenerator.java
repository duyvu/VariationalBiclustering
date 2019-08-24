package exp.gauss.cov;

import model.gauss.crossing.diffVar.ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkSampler;
import umontreal.iro.lecuyer.rng.LFSR113;

public class ColVaryingCovariatesDifferentVarGaussGenerator {
	public static void main(String[] args) throws Exception {

		ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkSampler sampler = new ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkSampler();

		String modelParameters = args[0];
		String delim = args[1];
		if (Integer.parseInt(delim) == 0)
			delim = " ";
		else
			delim = "\t";
		String confirmedModelFilename = args[2];
		String outputZFilename = args[3];
		String outputWFilename = args[4];
		String outputYFilename = args[5];
		int seed0 = Integer.parseInt(args[6].trim()) * 100;
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
		sampler.printModelParameters2File(confirmedModelFilename, delim);

		// Create a random generator
		LFSR113 randomGenerator = new LFSR113();
		randomGenerator.setSeed(seeds);

		// Generate a random network and save it to files
		System.out.println("Start generating the random network!");
		sampler.generateSampleWITHFixedCovariates(randomGenerator,
				outputZFilename, outputWFilename, outputYFilename, delim);
		System.out.println("Finish generating the random network!");

	}
}
