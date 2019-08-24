package model.gauss.crossing.sameVar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.StringTokenizer;

import exception.UnsupportedMethodException;

import graph.DenseGraph;

import model.TwoModeNetworkSampler;

import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.rng.RandomStream;

import utils.ArrayMethods;
import utils.PrintingMethods;
import utils.SamplingMethods;

public class NoCovariateCrossingDenseGaussTwoModeNetworkSampler extends
		NoCovariateCrossingDenseGaussTwoModeNetworkModel implements
		TwoModeNetworkSampler {

	@Override
	/* Read model parameters from a file */
	public void setParameters(String modelFileName, String delim) {
		BufferedReader modelReader;
		try {
			modelReader = new BufferedReader(new FileReader(new File(
					modelFileName)));

			StringTokenizer tokenizer = null;

			numOfRowGroups = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("numOfRowGroups = " + numOfRowGroups);
			// Read row mixing probabilities alpha
			alpha = new double[numOfRowGroups];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int k = 0; k < numOfRowGroups; k++)
				if (tokenizer.hasMoreTokens())
					alpha[k] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of row parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			System.out.println("alpha: ");
			ArrayMethods.printDoubleArray(alpha, " ", System.out);

			numOfColumnGroups = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("numOfColumnGroups = " + numOfColumnGroups);
			// Read column mixing probabilities beta
			beta = new double[numOfColumnGroups];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int l = 0; l < numOfColumnGroups; l++)
				if (tokenizer.hasMoreTokens())
					beta[l] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of column parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			System.out.println("beta: ");
			ArrayMethods.printDoubleArray(beta, " ", System.out);

			// Read theta parameters
			theta = new double[numOfRowGroups][numOfColumnGroups];
			for (int k = 0; k < numOfRowGroups; k++) {
				try {
					tokenizer = new StringTokenizer(modelReader.readLine(),
							delim);
				} catch (Exception e) {
					System.out
							.println("The number of row parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
				for (int l = 0; l < numOfColumnGroups; l++)
					if (tokenizer.hasMoreTokens())
						theta[k][l] = Double.parseDouble(tokenizer.nextToken());
					else {
						System.out
								.println("The number of column parameters is not matched to what was declared in file "
										+ modelFileName);
						System.exit(-1);
					}
			}
			System.out.println("theta: ");
			ArrayMethods.printDoubleMatrix(theta, " ", System.out);

			// Read sigma2 parameter
			try {
				tokenizer = new StringTokenizer(modelReader.readLine(), delim);
				if (tokenizer.hasMoreTokens())
					sigma2 = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out.println("Error when reading sigma2 in file "
							+ modelFileName);
					System.exit(-1);
				}
			} catch (Exception e) {
				System.out.println("Error when reading sigma2 in file "
						+ modelFileName);
				System.exit(-1);
			}
			System.out.println("sigma2: " + sigma2);

			modelReader.close();
		} catch (Exception e) {
			System.out.println("Can not read the model file " + modelFileName);
			System.exit(-1);
		}
	}

	@Override
	/*
	 * Print model parameters to a file which is a way to cross-check results
	 * later
	 */
	public void printModelParameters2File(String confirmedModelFilename,
			String delim) {

		try {
			PrintStream outputStream = new PrintStream(new File(
					confirmedModelFilename));

			// Write the number of row groups
			outputStream.println(numOfRowGroups);

			// Write alpha parameters
			PrintingMethods.printDoubleVector(alpha, outputStream, delim);

			// Write the number of column groups
			outputStream.println(numOfColumnGroups);

			// Write beta parameters
			PrintingMethods.printDoubleVector(beta, outputStream, delim);

			// Write theta parameters
			PrintingMethods.printDoubleMatrix(theta, outputStream, delim);

			// Write sigma2 parameter
			outputStream.println(sigma2);

			outputStream.close();

		} catch (Exception e) {
			System.out.println("Could not open the file "
					+ confirmedModelFilename);
			System.exit(-1);
		}

	}

	@Override
	/*
	 * Generate a sample and save it to files where there is no covariate on
	 * rows, columns, or edges
	 */
	public void generateSampleWITHOUTCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputYFilename, String delim) {

		// Generate row membership variables
		int[] rowMembership = new int[numOfRows];
		int[] rowStartingIndices = new int[numOfRowGroups];
		int[] rowGroupCounter = new int[numOfRowGroups];
		generateMembership(numOfRows, numOfRowGroups, alpha, sampler,
				rowMembership, rowStartingIndices, rowGroupCounter);

		// Generate column membership variables
		int[] columnMembership = new int[numOfColumns];
		int[] columnStartingIndices = new int[numOfColumnGroups];
		int[] columnGroupCounter = new int[numOfColumnGroups];
		generateMembership(numOfColumns, numOfColumnGroups, beta, sampler,
				columnMembership, columnStartingIndices, columnGroupCounter);

		// Generate edges in an efficient manner
		DenseGraph randomGraph = new DenseGraph(numOfRows, numOfColumns);
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				generateBlock(randomGraph, k, l, rowStartingIndices[k],
						rowGroupCounter[k], columnStartingIndices[l],
						columnGroupCounter[l], sampler);

		// Write every thing to the stream
		try {
			// Row membership variables
			PrintStream outputZ = new PrintStream(new File(outputZFilename));
			outputZ.println(rowMembership.length);
			for (int i = 0; i < rowMembership.length; i++)
				outputZ.println(rowMembership[i]);
			outputZ.close();
			// Column membership variables
			PrintStream outputW = new PrintStream(new File(outputWFilename));
			outputW.println(columnMembership.length);
			for (int j = 0; j < columnMembership.length; j++)
				outputW.println(columnMembership[j]);
			outputW.close();
			// Edge variables
			PrintStream outputY = new PrintStream(new File(outputYFilename));
			outputY.println(randomGraph.getNumOfRows());
			outputY.println(randomGraph.getNumOfColumns());
			for (int i = 0; i < numOfRows; i++) {
				for (int j = 0; j < numOfColumns; j++)
					outputY.print(randomGraph.getElement(i, j) + delim);
				outputY.println();
			}
			outputY.close();
		} catch (Exception e) {
			System.out
					.println("Can not open files to write the generated random network");
			System.exit(-1);
		}

	}

	protected void generateMembership(int numOfNodes, int numOfGroups,
			double[] p, RandomStream sampler, int[] membership,
			int[] startingIndices, int[] groupCounter) {
		ArrayMethods
				.copyIntegerArray(SamplingMethods.drawMultipleMultinomial(
						sampler, numOfNodes, p), groupCounter);
		int startIndex = 0;
		for (int k = 0; k < numOfGroups; k++) {
			startingIndices[k] = startIndex;
			int endIndex = startIndex + groupCounter[k];
			for (int i = startIndex; i < endIndex; i++)
				membership[i] = k;
			startIndex = endIndex;
		}
	}

	protected void generateBlock(DenseGraph randomGraph, int k, int l,
			int kStartIndex, int kN, int lStartIndex, int lN,
			RandomStream sampler) {

		System.out.println("Generating the block " + k + " x " + l);
		System.out.println("kStartIndex = " + kStartIndex);
		System.out.println("kN = " + kN);
		System.out.println("lStartIndex = " + lStartIndex);
		System.out.println("lN = " + lN);

		try {
			NormalGen generator = new NormalGen(sampler, theta[k][l],
					Math.sqrt(sigma2));

			for (int i = kStartIndex; i < (kStartIndex + kN); i++)
				for (int j = lStartIndex; j < (lStartIndex + lN); j++)
					randomGraph.setElement(i, j, generator.nextDouble());
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}

	public void generateSampleWITHRandomCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputRowCovFilename,
			String outputColCovFilename, String outputEdgeCovFilename,
			String outputYFilename, String delim) {
		System.out.println(this.getClass().getName()
				+ " does not support covariate structures!");
		System.exit(-1);
	}

	public void generateSampleWITHFixedCovariates(RandomStream sampler,
			String outputZFilename, String outputWFilename,
			String outputYFilename, String delim)
			throws UnsupportedMethodException {
		throw new UnsupportedMethodException();
	}
}
