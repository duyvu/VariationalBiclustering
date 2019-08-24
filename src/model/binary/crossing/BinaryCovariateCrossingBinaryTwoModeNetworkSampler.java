package model.binary.crossing;

import exception.UnsupportedMethodException;
import graph.MatrixElement;
import graph.SparseGraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.StringTokenizer;

import umontreal.iro.lecuyer.rng.RandomStream;
import utils.ArrayMethods;
import utils.PrintingMethods;
import utils.SamplingMethods;
import model.TwoModeNetworkSampler;

public class BinaryCovariateCrossingBinaryTwoModeNetworkSampler extends
		CovariateCrossingBinaryTwoModeNetworkModel implements
		TwoModeNetworkSampler {

	protected double[] rowBinaryRate;
	protected double[] columnBinaryRate;
	protected double[] edgeBinaryRate;

	@Override
	/* Read model parameters from a file */
	public void setParameters(String modelFileName, String delim) {
		BufferedReader modelReader;
		try {
			modelReader = new BufferedReader(new FileReader(new File(
					modelFileName)));

			StringTokenizer tokenizer = null;

			numOfRowGroups = Integer.parseInt(modelReader.readLine().trim());
			// Read row mixing probabilities beta
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

			numOfRowCovariates = Integer
					.parseInt(modelReader.readLine().trim());
			rowBeta = new double[numOfRowCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfRowCovariates; p++)
				if (tokenizer.hasMoreTokens())
					rowBeta[p] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of row beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			rowBinaryRate = new double[numOfRowCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfRowCovariates; p++)
				if (tokenizer.hasMoreTokens())
					rowBinaryRate[p] = Double
							.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of row beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}

			numOfColCovariates = Integer
					.parseInt(modelReader.readLine().trim());
			columnBeta = new double[numOfColCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfColCovariates; p++)
				if (tokenizer.hasMoreTokens())
					columnBeta[p] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of column beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			columnBinaryRate = new double[numOfColCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfColCovariates; p++)
				if (tokenizer.hasMoreTokens())
					columnBinaryRate[p] = Double.parseDouble(tokenizer
							.nextToken());
				else {
					System.out
							.println("The number of column beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}

			numOfEdgeCovariates = Integer.parseInt(modelReader.readLine()
					.trim());
			edgeBeta = new double[numOfEdgeCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfEdgeCovariates; p++)
				if (tokenizer.hasMoreTokens())
					edgeBeta[p] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of edge beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			edgeBinaryRate = new double[numOfEdgeCovariates];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int p = 0; p < numOfEdgeCovariates; p++)
				if (tokenizer.hasMoreTokens())
					edgeBinaryRate[p] = Double.parseDouble(tokenizer
							.nextToken());
				else {
					System.out
							.println("The number of edge beta is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}

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

			// Write row coefficients
			outputStream.println(numOfRowCovariates);
			PrintingMethods.printDoubleVector(rowBeta, outputStream, delim);
			PrintingMethods.printDoubleVector(rowBinaryRate, outputStream,
					delim);

			// Write column coefficients
			outputStream.println(numOfColCovariates);
			PrintingMethods.printDoubleVector(columnBeta, outputStream, delim);
			PrintingMethods.printDoubleVector(columnBinaryRate, outputStream,
					delim);

			// Write edge coefficients
			outputStream.println(numOfEdgeCovariates);
			PrintingMethods.printDoubleVector(edgeBeta, outputStream, delim);
			PrintingMethods.printDoubleVector(edgeBinaryRate, outputStream,
					delim);

			outputStream.close();

		} catch (Exception e) {
			System.out.println("Could not open the file "
					+ confirmedModelFilename);
			System.exit(-1);
		}

	}

	public void generateSampleWITHOUTCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputYFilename, String delim) {
		System.out.println(this.getClass().getName()
				+ " does not support covariate structures!");
		System.exit(-1);
	}

	@Override
	/*
	 * Generate a sample and save it to files where there is no covariate on
	 * rows, columns, or edges
	 */
	public void generateSampleWITHRandomCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputRowCovFilename,
			String outputColCovFilename, String outputEdgeCovFilename,
			String outputYFilename, String delim) {

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

		// Generate row covariates: rowCovariates
		rowCovariates = new double[numOfRows][numOfRowCovariates];
		drawSideCovariates(sampler, rowBinaryRate, rowCovariates);

		// Generate column covariates: columnCovariates
		columnCovariates = new double[numOfColumns][numOfColCovariates];
		drawSideCovariates(sampler, columnBinaryRate, columnCovariates);

		// Generate edge covariates: edgeCovariates;
		edgeCovariates = new double[numOfRows][numOfColumns][numOfEdgeCovariates];
		drawEdgeCovariates(sampler, edgeBinaryRate, edgeCovariates);

		// Generate edges
		SparseGraph randomGraph = new SparseGraph(numOfRows, numOfColumns);
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
			// Row covariates
			PrintStream outputRowCov = new PrintStream(new File(
					outputRowCovFilename));
			outputRowCov.println(numOfRowCovariates);
			ArrayMethods.printDoubleMatrix(rowCovariates, delim, outputRowCov);
			outputRowCov.close();
			// Column covariates
			PrintStream outputColCov = new PrintStream(new File(
					outputColCovFilename));
			outputColCov.println(numOfColCovariates);
			ArrayMethods.printDoubleMatrix(columnCovariates, delim,
					outputColCov);
			outputColCov.close();
			// Edge covariates
			PrintStream outputEdgeCov = new PrintStream(new File(
					outputEdgeCovFilename));
			outputEdgeCov.println(numOfEdgeCovariates);
			ArrayMethods.printDouble3DMatrix(edgeCovariates, delim,
					outputEdgeCov);
			outputEdgeCov.close();
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

	protected void drawSideCovariates(RandomStream sampler,
			double[] binaryRate, double[][] covariates) {
		for (int i = 0; i < covariates.length; i++)
			for (int p = 0; p < binaryRate.length; p++)
				covariates[i][p] = SamplingMethods.drawBernoulli(sampler,
						binaryRate[p]);
	}

	protected void drawEdgeCovariates(RandomStream sampler,
			double[] binaryRate, double[][][] covariates) {
		for (int i = 0; i < covariates.length; i++)
			for (int j = 0; j < covariates[0].length; j++)
				for (int p = 0; p < binaryRate.length; p++)
					covariates[i][j][p] = SamplingMethods.drawBernoulli(
							sampler, binaryRate[p]);
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

	protected void generateBlock(SparseGraph randomGraph, int k, int l,
			int kStartIndex, int kN, int lStartIndex, int lN,
			RandomStream sampler) {

		System.out.println("Generating the block " + k + " x " + l);
		System.out.println("kStartIndex = " + kStartIndex);
		System.out.println("kN = " + kN);
		System.out.println("lStartIndex = " + lStartIndex);
		System.out.println("lN = " + lN);

		try {
			for (int i = kStartIndex; i < (kStartIndex + kN); i++)
				for (int j = lStartIndex; j < (lStartIndex + lN); j++) {
					/* Compute the probability of the edge ij */
					double linCom = computeLinearCoefficientCombination(i, j,
							k, l);
					double expLinCom = Math.exp(linCom);
					double pi = expLinCom / (1.0 + expLinCom);
					/* Draw the edge value */
					if (SamplingMethods.drawBernoulli(sampler, pi) == 1)
						randomGraph.setElement(i, j, 1.0);
				}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}

	public void generateSampleWITHFixedCovariates(RandomStream sampler,
			String outputZFilename, String outputWFilename,
			String outputYFilename, String delim)
			throws UnsupportedMethodException {
		throw new UnsupportedMethodException();
	}

}
