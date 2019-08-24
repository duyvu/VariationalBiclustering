package exp.binary.nocov;

import java.io.File;
import java.io.FileFilter;
import java.io.PrintStream;

import em.TwoModeEMEngine;

import model.binary.additive.BinaryCovariateAdditiveBinaryTwoModeNetworkSampler;
import umontreal.iro.lecuyer.rng.LFSR113;

public class AdditiveSparseBinaryTwoModeNetworkSimulator {

	protected static BinaryCovariateAdditiveBinaryTwoModeNetworkSampler sampler = new BinaryCovariateAdditiveBinaryTwoModeNetworkSampler();

	protected static String delim;
	protected static int seed0;

	protected static String outputDir;
	protected static String subDirString;
	protected static String outputZFilename;
	protected static String outputWFilename;
	protected static String outputRowCovFilename;
	protected static String outputColCovFilename;
	protected static String outputEdgeCovFilename;
	protected static String outputYFilename;

	protected static TwoModeEMEngine emEngine = new TwoModeEMEngine();

	protected static String configurationEMFilename;
	protected static int numRowGroups;
	protected static int numColumnGroups;

	public static void main(String[] args) {

		generateData(args);
		runClustering();

	}

	protected static void generateData(String[] args) {

		// Model parameters to generate data
		String modelParameters = args[0];
		// How to parse the model parameter file
		delim = args[1];
		if (Integer.parseInt(delim) == 0)
			delim = " ";
		else
			delim = "\t";

		// Matrix size
		int numOfRows = Integer.parseInt(args[2]);
		int numOfColumns = Integer.parseInt(args[3]);

		// Sample ID which is the seed number for the random generator
		seed0 = 2 * Integer.parseInt(args[4].trim());
		int[] seeds = { 2, 8, 16, 128 };
		if (seed0 > 0)
			seeds[0] = -seed0;
		else
			seeds[0] = seed0;

		// Output directory
		outputDir = args[5].trim();

		// Configuration EM file
		configurationEMFilename = args[6].trim();

		subDirString = outputDir + "/" + seed0;
		boolean subDirSuccess = (new File(subDirString)).mkdirs();
		if (!subDirSuccess) {
			System.out.println("Can not create sub directory " + subDirString);
		}

		// Output files
		String confirmedParameters = subDirString + "/confirmedParameters.dat";
		outputZFilename = subDirString + "/Z.txt";
		outputWFilename = subDirString + "/W.txt";
		outputRowCovFilename = subDirString + "/RowCov.txt";
		outputColCovFilename = subDirString + "/ColCov.txt";
		outputEdgeCovFilename = subDirString + "/EdgeCov.txt";
		outputYFilename = subDirString + "/Y.txt";

		// Set model parameters
		System.out.println("Read model parameters");
		sampler.setParameters(modelParameters, delim);

		// Store the cluster configuration
		numRowGroups = sampler.getNumOfRowGroups();
		numColumnGroups = sampler.getNumOfColumnGroups();

		// Create a random generator
		LFSR113 randomGenerator = new LFSR113();
		randomGenerator.setSeed(seeds);

		// Write model parameters to confirm the simulation setting
		System.out.println("Write model parameters");
		sampler.printModelParameters2File(confirmedParameters, delim);

		// Generate a random network and save it to files
		System.out.println("Start generating the random network!");

		sampler.generateSampleWITHRandomCovariates(randomGenerator, numOfRows,
				numOfColumns, outputZFilename, outputWFilename,
				outputRowCovFilename, outputColCovFilename,
				outputEdgeCovFilename, outputYFilename, delim);
		System.out.println("Finish generating the random network!");
	}

	protected static void runClustering() {

		emEngine.configure(configurationEMFilename, "\t");

		emEngine.setDataFiles(TwoModeEMEngine.SPARSE_DATA_FORMAT, delim,
				outputYFilename, outputRowCovFilename, outputColCovFilename,
				outputEdgeCovFilename);

		emEngine.setNumOfRowGroups(numRowGroups);

		emEngine.setNumOfColumnGroups(numColumnGroups);

		emEngine.setRandomSeed(seed0, 0);

		emEngine.finalizeModel();

		String formatString = "%.20f";
		PrintStream outputStream = System.out;
		emEngine.runEM(outputStream, formatString);

		double[] misRates = emEngine.computeMisRates(outputZFilename,
				outputWFilename);
		System.out.println("misRates:");
		for (int k = 0; k < misRates.length; k++)
			System.out.println(misRates[k]);

		// Write every statistics to files
		try {
			PrintStream MisRatesOUT = new PrintStream(new File(outputDir
					+ "/MisRates-" + seed0 + ".txt"));
			for (int k = 0; k < misRates.length; k++)
				MisRatesOUT.print(misRates[k] + "\t");
			MisRatesOUT.close();

		} catch (Exception e) {
			System.out
					.println("Can not open files to write the generated random network");
			System.exit(-1);
		}

		// Clean the generate data files
		// List the files using our FileFilter
		File dir = new File(subDirString);
		File[] files = dir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File file) {
				if (file.getName().toLowerCase().endsWith(".txt")) {
					return true;
				}
				return false;
			}
		});
		for (File f : files) {
			System.out.println("Deleting file: " + f.getName());
			f.delete();
		}

	}
}
