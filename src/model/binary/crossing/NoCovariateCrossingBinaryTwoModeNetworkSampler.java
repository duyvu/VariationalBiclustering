package model.binary.crossing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.StringTokenizer;

import exception.UnsupportedMethodException;

import model.TwoModeNetworkSampler;

import umontreal.iro.lecuyer.randvar.BinomialGen;
import umontreal.iro.lecuyer.rng.RandomStream;

import utils.ArrayMethods;
import utils.PrintingMethods;
import utils.SamplingMethods;

import graph.MatrixElement;
import graph.SparseGraph;

public class NoCovariateCrossingBinaryTwoModeNetworkSampler extends
		NoCovariateCrossingBinaryFrequentistTwoModeNetworkModel implements
		TwoModeNetworkSampler {

	/* Read model parameters from a file */
	public void setParameters(String modelFileName, String delim) {
		BufferedReader modelReader;
		try {
			modelReader = new BufferedReader(new FileReader(new File(
					modelFileName)));

			StringTokenizer tokenizer = null;

			numOfRowGroups = Integer.parseInt(modelReader.readLine().trim());
			alpha = new double[numOfRowGroups];
			// Read row mixing probabilities beta
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

			numOfColumnGroups = Integer.parseInt(modelReader.readLine().trim());
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

			// Read theta parameters
			theta = new double[numOfRowGroups][numOfColumnGroups];
			for (int k = 0; k < numOfRowGroups; k++) {
				tokenizer = new StringTokenizer(modelReader.readLine(), delim);
				for (int l = 0; l < numOfColumnGroups; l++)
					if (tokenizer.hasMoreTokens())
						theta[k][l] = Double.parseDouble(tokenizer.nextToken());
					else {
						System.out
								.println("The number of row parameters is not matched to what was declared in file "
										+ modelFileName);
						System.exit(-1);
					}
			}

		} catch (Exception e) {
			System.out.println("Can not read the model file " + modelFileName);
			System.exit(-1);
		}
	}

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

			// Write pi1 matrix
			double[][] pi1 = new double[numOfRowGroups][numOfColumnGroups];
			computeP1(pi1);
			PrintingMethods.printDoubleMatrix(pi1, outputStream, delim);

			outputStream.close();

		} catch (Exception e) {
			System.out.println("Could not open the file "
					+ confirmedModelFilename);
			System.exit(-1);
		}

	}

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

		// System.out.println("rowGroupCounter: ");
		// for (int k = 0; k < rowGroupCounter.length; k++)
		// System.out.print(rowGroupCounter[k] + "\t");
		// System.out.println();

		// Generate column membership variables
		int[] columnMembership = new int[numOfColumns];
		int[] columnStartingIndices = new int[numOfColumnGroups];
		int[] columnGroupCounter = new int[numOfColumnGroups];
		generateMembership(numOfColumns, numOfColumnGroups, beta, sampler,
				columnMembership, columnStartingIndices, columnGroupCounter);

		// System.out.println("columnGroupCounter: ");
		// for (int l = 0; l < columnGroupCounter.length; l++)
		// System.out.print(columnGroupCounter[l] + "\t");
		// System.out.println();

		// Compute the probability matrix pi(Y = 1)
		double[][] pi1 = new double[numOfRowGroups][numOfColumnGroups];
		computeP1(pi1);

		// Generate edges in an efficient manner
		SparseGraph randomGraph = new SparseGraph(numOfRows, numOfColumns);
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				generateBlock(randomGraph, k, l, rowStartingIndices[k],
						rowGroupCounter[k], columnStartingIndices[l],
						columnGroupCounter[l], pi1, sampler);

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
			outputY.println(randomGraph.getNumberOfNonzeroEdges());
			Iterator<MatrixElement> it = randomGraph.elementIterator();
			if (it != null)
				for (; it.hasNext();) {
					MatrixElement element = it.next();
					int i = element.getRow();
					int j = element.getCol();
					outputY.println(i + delim + j);
				}
			outputY.close();
		} catch (Exception e) {
			System.out
					.println("Can not open files to write the generated random network");
			System.exit(-1);
		}

	}

	protected void computeP1(double[][] pi1) {
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				double norm = Math.exp(theta[k][l]);
				pi1[k][l] = norm / (1 + norm);
			}
	}

	protected void generateMembership(int numOfNodes, int numOfGroups,
			double[] p, RandomStream sampler, int[] membership,
			int[] startingIndices, int[] groupCounter) {
		ArrayMethods
				.copyIntegerArray(SamplingMethods.drawMultipleMultinomial(
						sampler, numOfNodes, p), groupCounter);

		// System.out.println("groupCounter: ");
		// for (int k = 0; k < groupCounter.length; k++)
		// System.out.print(groupCounter[k] + "\t");

		int startIndex = 0;
		for (int k = 0; k < numOfGroups; k++) {
			startingIndices[k] = startIndex;
			int endIndex = startIndex + groupCounter[k];
			for (int i = startIndex; i < endIndex; i++)
				membership[i] = k;
			startIndex = endIndex;
		}
	}

	protected void generateBlock(SparseGraph randomGraph, int k, int l,
			int kStartIndex, int kN, int lStartIndex, int lN, double[][] pi1,
			RandomStream sampler) {

		System.out.println("Generating the block " + k + " x " + l);
		System.out.println("kStartIndex = " + kStartIndex);
		System.out.println("kN = " + kN);
		System.out.println("lStartIndex = " + lStartIndex);
		System.out.println("lN = " + lN);

		// Draw the number of non-zero edges
		int potentialEdges = kN * lN;
		int nonzeroEdges = BinomialGen.nextInt(sampler, potentialEdges,
				pi1[k][l]);

		System.out.println("potentialEdges = " + potentialEdges);
		System.out.println("nonzeroEdges = " + nonzeroEdges);

		// Select each non-zero edge randomly using the accept-reject method
		for (int edgeIndex = 0; edgeIndex < nonzeroEdges; edgeIndex++) {
			try {
				// randomly select non-zero edge y_ij
				int i, j;
				do {
					i = kStartIndex
							+ (int) Math.floor(sampler.nextDouble() * kN);
					j = lStartIndex
							+ (int) Math.floor(sampler.nextDouble() * lN);
					// System.out.println("Propose " + i + "\t" + j);
				} while (randomGraph.getElement(i, j) != 0);
				// System.out.println("Accept the edge " + i + " -> " + j);
				randomGraph.setElement(i, j, 1.0);
			} catch (Exception e) {
				System.out.println("This exception could not happen since"
						+ randomGraph.getClass().getName()
						+ " supports setElement()");
				e.printStackTrace();
			}
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
