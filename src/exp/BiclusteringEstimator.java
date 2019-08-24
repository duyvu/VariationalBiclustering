package exp;

import java.io.PrintStream;

import em.TwoModeEMEngine;

public class BiclusteringEstimator {

	public static void main(String[] args) {

		String configurationFilename = args[0];
		String delim = args[1];
		if (Integer.parseInt(delim) == 0)
			delim = " ";
		else
			delim = "\t";
		int numOfRowGroups = Integer.parseInt(args[2]);
		int numOfColumnGroups = Integer.parseInt(args[3]);
		int seed0 = Integer.parseInt(args[4].trim());
		if (seed0 > 0)
			seed0 = -seed0;

		TwoModeEMEngine emEngine = new TwoModeEMEngine();

		emEngine.configure(configurationFilename, delim);

		emEngine.setNumOfRowGroups(numOfRowGroups);

		emEngine.setNumOfColumnGroups(numOfColumnGroups);

		emEngine.setRandomSeed(seed0, 0);

		emEngine.finalizeModel();

		String formatString = "%.20f";
		PrintStream outputStream = System.out;
		emEngine.runEM(outputStream, formatString);
	}

}
