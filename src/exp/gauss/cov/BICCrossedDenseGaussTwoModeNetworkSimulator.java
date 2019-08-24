package exp.gauss.cov;

import java.io.File;
import java.io.FileFilter;
import java.io.PrintStream;

import em.TwoModeEMEngine;

public class BICCrossedDenseGaussTwoModeNetworkSimulator extends
		CrossedDenseGaussTwoModeNetworkSimulator {

	public static void main(String[] args) {

		generateData(args);

		runClustering();

	}

	protected static void runClustering() {

		// Write every statistics to files
		try {
			PrintStream BICsOUT = new PrintStream(new File(outputDir + "/BICs-"
					+ seed0 + ".txt"));
			PrintStream bestBICOUT = new PrintStream(new File(outputDir
					+ "/bestBIC-" + seed0 + ".txt"));

			double bestRowGroups = -1;
			double bestColGroups = -1;
			double bestBIC = Double.MAX_VALUE;

			for (int runningRowGroups = 1; runningRowGroups <= 2 * numRowGroups; runningRowGroups++) {
				for (int runningColGroups = 1; runningColGroups <= 2 * numColumnGroups; runningColGroups++) {

					emEngine = new TwoModeEMEngine();

					emEngine.configure(configurationEMFilename, "\t");

					emEngine.setDataFiles(TwoModeEMEngine.MATRIX_DATA_FORMAT,
							delim, outputYFilename, outputRowCovFilename,
							outputColCovFilename, outputEdgeCovFilename);

					emEngine.setNumOfRowGroups(runningRowGroups);

					emEngine.setNumOfColumnGroups(runningColGroups);

					emEngine.setRandomSeed(seed0, 0);

					emEngine.finalizeModel();

					String formatString = "%.20f";
					PrintStream outputStream = System.out;
					emEngine.runEM(outputStream, formatString);

					double BIC = emEngine.computeBIC();
					System.out.println("BIC[\t" + runningRowGroups + ",\t"
							+ runningColGroups + "]:\t" + BIC);
					BICsOUT.println(runningRowGroups + "," + runningColGroups
							+ "," + BIC);

					if (BIC < bestBIC) {
						bestBIC = BIC;
						bestRowGroups = runningRowGroups;
						bestColGroups = runningColGroups;
					}

				}
			}
			bestBICOUT.print(bestRowGroups + "," + bestColGroups + ","
					+ bestBIC);

			// Close output files
			BICsOUT.close();
			bestBICOUT.close();

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

		} catch (Exception e) {
			System.out
					.println("Can not open files to write the simulation results");
			System.exit(-1);
		}

	}

}
