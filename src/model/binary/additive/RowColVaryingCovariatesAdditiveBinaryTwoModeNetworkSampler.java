package model.binary.additive;

import exception.UnsupportedMethodException;
import graph.MatrixElement;
import graph.SparseAccess;
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

public class RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkSampler extends
		RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkModel implements
		TwoModeNetworkSampler {

	protected String rowCovsFile;
	protected String rowCoefsFile;

	protected String columnCovsFile;
	protected String columnCoefsFile;

	protected int numOfRows;
	protected int numOfColumns;

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

			// Write the number of rows
			outputStream.println(numOfRows);

			// Write the number of row groups
			outputStream.println(numOfRowGroups);

			// Write alpha parameters
			PrintingMethods.printDoubleVector(alpha, outputStream, delim);

			// Write theta parameters
			PrintingMethods.printDoubleVector(theta, outputStream, delim);

			// Write the number of columns
			outputStream.println(numOfColumns);

			// Write the number of column groups
			outputStream.println(numOfColumnGroups);

			// Write beta parameters
			PrintingMethods.printDoubleVector(beta, outputStream, delim);

			// Write gamma parameters
			PrintingMethods.printDoubleVector(gamma, outputStream, delim);

			// Write the path to row covariates and coefficients
			outputStream.println(rowCovsFile);
			outputStream.println(rowCoefsFile);

			// Write the path to column covariates and coefficients
			outputStream.println(columnCovsFile);
			outputStream.println(columnCoefsFile);

			outputStream.close();

		} catch (Exception e) {
			System.out.println("Could not open the file "
					+ confirmedModelFilename);
			System.exit(-1);
		}

	}

	public void setParameters(String modelFileName, String delim) {
		BufferedReader modelReader;
		try {
			modelReader = new BufferedReader(new FileReader(new File(
					modelFileName)));

			StringTokenizer tokenizer = null;

			System.out.println("Reading rows");

			numOfRows = Integer.parseInt(modelReader.readLine().trim());
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
			// Read theta parameters
			theta = new double[numOfRowGroups];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int k = 0; k < numOfRowGroups; k++)
				if (tokenizer.hasMoreTokens())
					theta[k] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of row parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}

			System.out.println("Reading columns");

			numOfColumns = Integer.parseInt(modelReader.readLine().trim());
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
			// Read gamma parameters
			gamma = new double[numOfColumnGroups];
			tokenizer = new StringTokenizer(modelReader.readLine(), delim);
			for (int l = 0; l < numOfColumnGroups; l++)
				if (tokenizer.hasMoreTokens())
					gamma[l] = Double.parseDouble(tokenizer.nextToken());
				else {
					System.out
							.println("The number of column parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}

			System.out.println("Reading rows");

			/* Read the row covariate file */
			rowCovsFile = modelReader.readLine().trim();
			rowCoefsFile = modelReader.readLine().trim();
			if (rowCovsFile != null && rowCoefsFile != null) {

				System.out.println("Reading the row covariate file "
						+ rowCovsFile);

				BufferedReader rowCovReader = new BufferedReader(
						new FileReader(new File(rowCovsFile)));

				// Read the number of row covariates
				numOfRowCovariates = Integer.parseInt(rowCovReader.readLine()
						.trim());

				if (numOfRowCovariates > 0) {
					System.out.println("Number of Row Covariates is "
							+ numOfRowCovariates);
					rowCovariates = new double[numOfRows][numOfRowCovariates];
					// Read rows one by one
					ArrayMethods.readMatrix(rowCovReader, delim, rowCovariates);
				}

				rowCovReader.close();

				System.out.println("Reading row coefficients");

				if (numOfRowCovariates > 0) {
					rowBeta = new double[numOfColumns][numOfRowCovariates];
					BufferedReader rowCoefsReader = new BufferedReader(
							new FileReader(new File(rowCoefsFile)));
					ArrayMethods.readMatrix(rowCoefsReader, delim, rowBeta);
					rowCoefsReader.close();
				}

			}

			System.out.println("Reading columns");

			/* Read the column covariate file */
			columnCovsFile = modelReader.readLine().trim();
			columnCoefsFile = modelReader.readLine().trim();
			if (columnCovsFile != null && columnCoefsFile != null) {

				System.out.println("Reading the column covariate file "
						+ columnCovsFile);

				BufferedReader colCovReader = new BufferedReader(
						new FileReader(new File(columnCovsFile)));

				// Read the number of row covariates
				numOfColCovariates = Integer.parseInt(colCovReader.readLine()
						.trim());

				if (numOfColCovariates > 0) {
					System.out.println("Number of Column Covariates is "
							+ numOfColCovariates);
					columnCovariates = new double[numOfColumns][numOfColCovariates];
					// Read rows one by one
					ArrayMethods.readMatrix(colCovReader, delim,
							columnCovariates);
				}

				colCovReader.close();

				System.out.println("Reading column coefficients");

				if (numOfColCovariates > 0) {
					columnBeta = new double[numOfRows][numOfColCovariates];
					BufferedReader colCoefsReader = new BufferedReader(
							new FileReader(new File(columnCoefsFile)));
					ArrayMethods.readMatrix(colCoefsReader, delim, columnBeta);
					colCoefsReader.close();
				}

			}

		} catch (Exception e) {
			System.out.println("Can not read the model file " + modelFileName);
			System.exit(-1);
		}
	}

	public void generateSampleWITHOUTCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputYFilename, String delim) {
		System.out.println("This method is not supported by "
				+ this.getClass().getName());
	}

	public void generateSampleWITHRandomCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputRowCovFilename,
			String outputColCovFilename, String outputEdgeCovFilename,
			String outputYFilename, String delim) {
		System.out.println("This method is not supported by "
				+ this.getClass().getName());
	}

	protected boolean checkDegrees(SparseGraph randomGraph) {

		double[] userDegrees = new double[numOfRows];
		ArrayMethods.resetDoubleArray(userDegrees);

		Iterator<MatrixElement> it = ((SparseAccess) randomGraph)
				.elementIterator();
		if (it != null)
			for (; it.hasNext();) {
				MatrixElement element = it.next();
				int i = element.getRow();
				userDegrees[i]++;
			}

		for (int i = 0; i < numOfRows; i++)
			if (userDegrees[i] < 20)
				return false;

		return true;
	}

	public void generateSampleWITHFixedCovariates(RandomStream sampler,
			String outputZFilename, String outputWFilename,
			String outputYFilename, String delim)
			throws UnsupportedMethodException {

		int[] rowMembership = new int[numOfRows];
		int[] columnMembership = new int[numOfColumns];
		SparseGraph randomGraph;

		int trialCount = 0;
		do {

			trialCount++;
			System.out.println("Running the trial " + trialCount);

			// Generate row membership variables
			generateMembership(numOfRows, alpha, sampler, rowMembership);

			// Generate column membership variables
			generateMembership(numOfColumns, beta, sampler, columnMembership);

			// Generate edges
			randomGraph = new SparseGraph(numOfRows, numOfColumns);
			for (int i = 0; i < numOfRows; i++)
				for (int j = 0; j < numOfColumns; j++)
					generateEdge(randomGraph, i, j, rowMembership[i],
							columnMembership[j], sampler);

		} while (!checkDegrees(randomGraph));

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

	protected void generateMembership(int numOfNodes, double[] p,
			RandomStream sampler, int[] membership) {
		for (int i = 0; i < numOfNodes; i++)
			membership[i] = SamplingMethods.drawSingleMultinomial(sampler, p);
	}

	protected void generateEdge(SparseGraph randomGraph, int i, int j, int k,
			int l, RandomStream sampler) {

		// System.out.println("Generating the edge " + i + " - " + j);

		/* Compute the probability of the edge ij */
		double linCom = computeLinearCoefficientCombination(i, j, k, l);
		double expLinCom = Math.exp(linCom);
		double pi = expLinCom / (1.0 + expLinCom);
		/* Draw the edge value */
		if (SamplingMethods.drawBernoulli(sampler, pi) == 1)
			try {
				randomGraph.setElement(i, j, 1.0);
			} catch (UnsupportedMethodException e) {
				e.printStackTrace();
			}

	}

}
