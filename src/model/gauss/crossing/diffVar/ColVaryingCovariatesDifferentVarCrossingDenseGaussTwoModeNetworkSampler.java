package model.gauss.crossing.diffVar;

import exception.UnsupportedMethodException;
import graph.DenseGraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.StringTokenizer;

import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.rng.RandomStream;
import utils.ArrayMethods;
import utils.PrintingMethods;
import utils.SamplingMethods;
import model.TwoModeNetworkSampler;

public class ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkSampler
		extends
		ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel
		implements TwoModeNetworkSampler {

	protected String rowCovsFile;
	protected String rowCoefsFile;

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

			// Write the number of columns
			outputStream.println(numOfColumns);

			// Write the number of column groups
			outputStream.println(numOfColumnGroups);

			// Write beta parameters
			PrintingMethods.printDoubleVector(beta, outputStream, delim);

			// Write theta parameters
			PrintingMethods.printDoubleMatrix(theta, outputStream, delim);

			// Write sigma2 parameter
			PrintingMethods.printDoubleMatrix(sigma2, outputStream, delim);

			// Write the path to row covariates and coefficients
			outputStream.println(rowCovsFile);
			outputStream.println(rowCoefsFile);

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

			System.out.println("Reading row parameters");

			numOfRows = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("Number of Rows is " + numOfRows);
			numOfRowGroups = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("Number of Row Groups is " + numOfRowGroups);
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

			System.out.println("Reading column parameters");

			numOfColumns = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("Number of Columns is " + numOfColumns);
			numOfColumnGroups = Integer.parseInt(modelReader.readLine().trim());
			System.out.println("Number of Column Groups is "
					+ numOfColumnGroups);

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

			System.out.println("Reading theta parameters");

			// Read theta parameters
			theta = new double[numOfRowGroups][numOfColumnGroups];
			for (int k = 0; k < numOfRowGroups; k++) {
				String line = modelReader.readLine();
				if (line != null) {
					tokenizer = new StringTokenizer(line, delim);
					for (int l = 0; l < numOfColumnGroups; l++) {
						if (tokenizer.hasMoreTokens())
							theta[k][l] = Double.parseDouble(tokenizer
									.nextToken());
						else {
							System.out
									.println("The number of column parameters is not matched to what was declared in file "
											+ modelFileName);
							System.exit(-1);
						}
					}
				} else {
					System.out
							.println("The number of row parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			}

			System.out.println("Reading sigma2 parameters");

			// Read theta parameters
			sigma2 = new double[numOfRowGroups][numOfColumnGroups];
			for (int k = 0; k < numOfRowGroups; k++) {
				String line = modelReader.readLine();
				if (line != null) {
					tokenizer = new StringTokenizer(line, delim);
					for (int l = 0; l < numOfColumnGroups; l++) {
						if (tokenizer.hasMoreTokens())
							sigma2[k][l] = Double.parseDouble(tokenizer
									.nextToken());
						else {
							System.out
									.println("The number of column parameters is not matched to what was declared in file "
											+ modelFileName);
							System.exit(-1);
						}
					}
				} else {
					System.out
							.println("The number of row parameters is not matched to what was declared in file "
									+ modelFileName);
					System.exit(-1);
				}
			}

			System.out.println("Reading row covariates and coefficients");

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

	public void generateSampleWITHFixedCovariates(RandomStream sampler,
			String outputZFilename, String outputWFilename,
			String outputYFilename, String delim)
			throws UnsupportedMethodException {

		int[] rowMembership = new int[numOfRows];
		int[] columnMembership = new int[numOfColumns];
		DenseGraph randomGraph = new DenseGraph(numOfRows, numOfColumns);

		// Generate row membership variables
		generateMembership(numOfRows, alpha, sampler, rowMembership);

		// Generate column membership variables
		generateMembership(numOfColumns, beta, sampler, columnMembership);

		// Generate edges
		for (int i = 0; i < numOfRows; i++)
			for (int j = 0; j < numOfColumns; j++)
				generateEdge(randomGraph, i, j, rowMembership[i],
						columnMembership[j], sampler);

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

	protected void generateMembership(int numOfNodes, double[] p,
			RandomStream sampler, int[] membership) {
		for (int i = 0; i < numOfNodes; i++)
			membership[i] = SamplingMethods.drawSingleMultinomial(sampler, p);
	}

	protected void generateEdge(DenseGraph randomGraph, int i, int j, int k,
			int l, RandomStream sampler) {

		// System.out.println("Generating the edge " + i + " - " + j);

		/* Draw the edge value */
		double mu = computeLinearCoefficientCombination(i, j, k, l);
		try {
			randomGraph.setElement(i, j,
					NormalGen.nextDouble(sampler, mu, Math.sqrt(sigma2[k][l])));
		} catch (UnsupportedMethodException e) {
			e.printStackTrace();
		}

	}

}
