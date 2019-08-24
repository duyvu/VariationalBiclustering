package model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.StringTokenizer;

import em.EMEngine;
import exception.UnsupportedMethodException;

import utils.ArrayMethods;

import graph.Graph;
import graph.MatrixElement;

public abstract class TwoModeNetworkModel extends NetworkModel {

	public static final double logDoubleMax = Math.log(Double.MAX_VALUE) - 1;

	protected int numOfRowGroups;
	protected double[][] z;
	protected double[][] prevZ;
	protected double[] alpha;
	protected double[] prevAlpha;

	protected int numOfColumnGroups;
	protected double[][] w;
	protected double[][] prevW;
	protected double[] beta;
	protected double[] prevBeta;

	protected Graph graph;

	protected int numOfRowCovariates;
	protected double[][] rowCovariates;

	protected int numOfColCovariates;
	protected double[][] columnCovariates;

	protected int numOfEdgeCovariates;
	protected double[][][] edgeCovariates;

	public abstract void setHyperParameters(
			HashMap<String, String> configurationMap);

	/* Set network data and all other cached data like degree sequences */
	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {

		System.out.println("Network data file is " + networkDataFile);

		try {
			BufferedReader networkReader = new BufferedReader(new FileReader(
					new File(networkDataFile)));

			// Read the number of rows
			int numOfRows = Integer.parseInt(networkReader.readLine().trim());
			System.out.println("Number of Rows is " + numOfRows);

			// Read the number of columns
			int numOfColumns = Integer
					.parseInt(networkReader.readLine().trim());
			System.out.println("Number of Columns is " + numOfColumns);

			// If the sparse format is used
			if (dataFormat == EMEngine.SPARSE_DATA_FORMAT) {

				LinkedList<MatrixElement> edges = new LinkedList<MatrixElement>();

				// Read the number of non-zero edges
				int numOfEdges = Integer.parseInt(networkReader.readLine()
						.trim());

				// Read edges one by one
				readEdgeList(networkReader, dataDelim, numOfEdges, edges);

				// Create the data graph which depends on the subclasses
				createGraph(numOfRows, numOfColumns, edges);
			}
			// If the matrix format is used
			else if (dataFormat == EMEngine.MATRIX_DATA_FORMAT) {

				double[][] edges = new double[numOfRows][numOfColumns];

				// Read rows one by one
				readMatrix(networkReader, dataDelim, edges);

				// Create the data graph which depends on the subclasses
				createGraph(numOfRows, numOfColumns, edges);
			} else {
				System.out
						.println("Data format is not supported! Data format must be 0 for matrix and 1 for sparse data.");
				System.exit(-1);
			}

			networkReader.close();

			/* Read the row covariate file */
			if (rowCovFile != null) {

				System.out.println("Reading the row covariate file "
						+ rowCovFile);

				BufferedReader rowCovReader = new BufferedReader(
						new FileReader(new File(rowCovFile)));

				// Read the number of row covariates
				numOfRowCovariates = Integer.parseInt(rowCovReader.readLine()
						.trim());

				if (numOfRowCovariates > 0) {
					System.out.println("Number of Row Covariates is "
							+ numOfRowCovariates);
					rowCovariates = new double[numOfRows][numOfRowCovariates];
					// Read rows one by one
					ArrayMethods.readMatrix(rowCovReader, dataDelim,
							rowCovariates);
				}

				rowCovReader.close();
			}

			/* Read the column covariate file */
			if (colCovFile != null) {

				System.out.println("Reading the column covariate file "
						+ colCovFile);

				BufferedReader colCovReader = new BufferedReader(
						new FileReader(new File(colCovFile)));

				// Read the number of row covariates
				numOfColCovariates = Integer.parseInt(colCovReader.readLine()
						.trim());

				if (numOfColCovariates > 0) {
					System.out.println("Number of Column Covariates is "
							+ numOfColCovariates);
					columnCovariates = new double[numOfColumns][numOfColCovariates];
					// Read rows one by one
					ArrayMethods.readMatrix(colCovReader, dataDelim,
							columnCovariates);
				}

				colCovReader.close();
			}

			/* Read the edge covariate file */
			if (edgeCovFile != null) {

				System.out.println("Reading the edge covariate file "
						+ edgeCovFile);

				BufferedReader edgeCovReader = new BufferedReader(
						new FileReader(new File(edgeCovFile)));

				// Read the number of row covariates
				numOfEdgeCovariates = Integer.parseInt(edgeCovReader.readLine()
						.trim());

				if (numOfEdgeCovariates > 0) {
					System.out.println("Number of Edge Covariates is "
							+ numOfEdgeCovariates);
					edgeCovariates = new double[numOfRows][numOfColumns][numOfEdgeCovariates];
					// Read rows one by one
					ArrayMethods.read3DMatrix(edgeCovReader, dataDelim,
							edgeCovariates);
				}

				edgeCovReader.close();
			}

		} catch (Exception e) {
			System.out.println("I/O errors when reading the data file "
					+ networkDataFile);
			e.printStackTrace();
			System.exit(-1);
		}

	}

	/*
	 * Create the data graph which depends on the subclasses. If the subclasses
	 * are sparse models, sparse graphs are created
	 */
	protected abstract void createGraph(int numOfRows, int numOfColumns,
			LinkedList<MatrixElement> edges);

	/*
	 * Create the data graph which depends on the subclasses. If the subclasses
	 * are sparse models, sparse graphs are created
	 */
	protected abstract void createGraph(int numOfRows, int numOfColumns,
			double[][] edges);

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {

		z = new double[graph.getNumOfRows()][numOfRowGroups];
		prevZ = new double[graph.getNumOfRows()][numOfRowGroups];
		alpha = new double[numOfRowGroups];
		prevAlpha = new double[numOfRowGroups];

		w = new double[graph.getNumOfColumns()][numOfColumnGroups];
		prevW = new double[graph.getNumOfColumns()][numOfColumnGroups];
		beta = new double[numOfColumnGroups];
		prevBeta = new double[numOfColumnGroups];
	}

	/* Reset model parameters and every data structure for a new EM run */
	public void reset(boolean isRandom, boolean isMMinMStep) {
		if (isRandom) {
			System.out.println("RANDOMLY: Reset the model for a new EM run!");

			// Reset row membership variables
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					z[i][k] = randomGenerator.nextDouble();
			normalizeMembership(z);

			// Reset column membership variables
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					w[j][l] = randomGenerator.nextDouble();
			normalizeMembership(w);

		} else {
			System.out.println("EQUAL INIT: Reset the model for a new EM run!");

			// Reset row membership variables
			double rowNoise = .1 / numOfRowGroups;
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					z[i][k] = 1.0 / numOfRowGroups + rowNoise
							* (randomGenerator.nextDouble() - 0.5);
				}

			// Reset column membership variables
			double columnNoise = .1 / numOfColumnGroups;
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++) {
					w[j][l] = 1.0 / numOfColumnGroups + columnNoise
							* (randomGenerator.nextDouble() - 0.5);
				}

		}

		// ArrayMethods.printDoubleMatrix(z, " ", System.out);
		// ArrayMethods.printDoubleMatrix(w, " ", System.out);

		// Update hyper-parameters of membership variables
		updateLatentHyperParameters();

		// ArrayMethods.printDoubleArray(alpha, " ", System.out);
		// ArrayMethods.printDoubleArray(beta, " ", System.out);

		// Reset model parameters
		resetModelParamters();

		runMStep(0, isMMinMStep);

	}

	/* Reset model parameters */
	protected abstract void resetModelParamters();

	/* Compute the likelihood lower bound */
	public abstract double getLikelihoodLowerBound();

	protected abstract double computeSumOfLogProbs();

	public abstract void updateCachedParameters();

	/* Run a full EM iteration */
	public void runEMIteration(int emIteration, int maxESteps,
			double membershipPrecision, boolean isMMinEStep, boolean isMMinMStep) {
		System.out
				.println("runEMIteration(...) must be implemented for the model class "
						+ this.getClass().getName() + "!!!");
		System.exit(-1);
	}

	/* Update latent membership variables and their hyper-parameters */
	public void runEStep(int emIteration, int eStep, boolean isMM) {

		/* Refresh all cached parameters */
		updateCachedParameters();

		/*
		 * Update latent membership variables for rows and columns. Two options
		 * for E step: using MM updates or FP updates
		 */
		if (isMM) {
			runMM_EStep(emIteration, eStep);

		} else {
			runFP_EStep(emIteration, eStep);
		}

		/* Update hyper-parameters of latent groups */
		updateLatentHyperParameters();
	}

	protected abstract void runFP_EStep(int emIteraction, int eStep);

	protected abstract void runMM_EStep(int emIteraction, int eStep);

	protected void normalizeLogMembership2Membership(double[][] tau) {

		int numOfNodes = tau.length;
		int numOfGroups = tau[0].length;

		/* Sliding values to avoid overflow */
		for (int i = 0; i < numOfNodes; i++) {
			// First it holds the max value
			double slidingValue = tau[i][0];
			for (int k = 1; k < numOfGroups; k++)
				if (slidingValue < tau[i][k])
					slidingValue = tau[i][k];
			// Now it actually holds the sliding value
			slidingValue = logDoubleMax - slidingValue;
			// Slide values
			for (int k = 0; k < numOfGroups; k++)
				tau[i][k] += slidingValue;
		}

		/* Do the normalization */
		for (int i = 0; i < numOfNodes; i++) {

			double denominator = 0;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] = Math.exp(tau[i][k]);
				denominator += tau[i][k];
			}

			// Make sure all membership variables are greater than
			// membershipPrecision
			boolean again = false;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] /= denominator;
				if (tau[i][k] < minMembership) {
					tau[i][k] = minMembership;
					again = true;
				}
			}

			if (again) {
				denominator = 0;
				for (int k = 0; k < numOfGroups; k++) {
					denominator += tau[i][k];
				}
				for (int k = 0; k < numOfGroups; k++)
					tau[i][k] /= denominator;
			}

		}
	}

	protected void normalizeLogMembership2Membership(double[][] tau,
			int fromRow, int toRow) {

		int numOfGroups = tau[0].length;

		/* Sliding values to avoid overflow */
		for (int i = fromRow; i < toRow; i++) {
			// First it holds the max value
			double slidingValue = tau[i][0];
			for (int k = 1; k < numOfGroups; k++)
				if (slidingValue < tau[i][k])
					slidingValue = tau[i][k];
			// Now it actually holds the sliding value
			slidingValue = logDoubleMax - slidingValue;
			// Slide values
			for (int k = 0; k < numOfGroups; k++)
				tau[i][k] += slidingValue;
		}

		/* Do the normalization */
		for (int i = fromRow; i < toRow; i++) {

			double denominator = 0;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] = Math.exp(tau[i][k]);
				denominator += tau[i][k];
			}

			// Make sure all membership variables are greater than
			// membershipPrecision
			boolean again = false;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] /= denominator;
				if (tau[i][k] < minMembership) {
					tau[i][k] = minMembership;
					again = true;
				}
			}

			if (again) {
				denominator = 0;
				for (int k = 0; k < numOfGroups; k++) {
					denominator += tau[i][k];
				}
				for (int k = 0; k < numOfGroups; k++)
					tau[i][k] /= denominator;
			}

		}
	}

	protected void normalizeMembership(double[][] tau) {

		int numOfNodes = tau.length;
		int numOfGroups = tau[0].length;

		for (int i = 0; i < numOfNodes; i++) {

			double denominator = 0;
			for (int k = 0; k < numOfGroups; k++)
				denominator += tau[i][k];

			// Make sure all membership variables are greater than
			// membershipPrecision
			boolean again = false;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] /= denominator;
				if (tau[i][k] < minMembership) {
					tau[i][k] = minMembership;
					again = true;
				}
			}

			if (again) {
				denominator = 0;
				for (int k = 0; k < numOfGroups; k++) {
					denominator += tau[i][k];
				}
				for (int k = 0; k < numOfGroups; k++)
					tau[i][k] /= denominator;
			}

		}
	}

	protected void normalizeMembership(double[][] tau, int fromRow, int toRow) {

		int numOfGroups = tau[fromRow].length;

		for (int i = fromRow; i < toRow; i++) {

			double denominator = 0;
			for (int k = 0; k < numOfGroups; k++)
				denominator += tau[i][k];

			// Make sure all membership variables are greater than
			// membershipPrecision
			boolean again = false;
			for (int k = 0; k < numOfGroups; k++) {
				tau[i][k] /= denominator;
				if (tau[i][k] < minMembership) {
					tau[i][k] = minMembership;
					again = true;
				}
			}

			if (again) {
				denominator = 0;
				for (int k = 0; k < numOfGroups; k++) {
					denominator += tau[i][k];
				}
				for (int k = 0; k < numOfGroups; k++)
					tau[i][k] /= denominator;
			}

		}
	}

	public boolean areMembershipVariablesChangedSignificantly(
			double tauPrecision) {

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				if (Math.abs(z[i][k] - prevZ[i][k]) > tauPrecision)
					return true;

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs(w[j][l] - prevW[j][l]) > tauPrecision)
					return true;

		return false;
	}

	public double getLargestMembershipChange() {

		double largestMembershipChange = 0;

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				if (Math.abs(z[i][k] - prevZ[i][k]) > largestMembershipChange)
					largestMembershipChange = Math.abs(z[i][k] - prevZ[i][k]);

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs(w[j][l] - prevW[j][l]) > largestMembershipChange)
					largestMembershipChange = Math.abs(w[j][l] - prevW[j][l]);

		return largestMembershipChange;
	}

	public abstract void updateLatentHyperParameters();

	/* Update model parameters and their hyper-parameters */
	public void runMStep(int emIteraction, boolean isMM) {

		/* Refresh all cached parameters */
		updateCachedParameters();

		/*
		 * Update model parameters with two options: using MM updates or
		 * Gradient updates
		 */
		if (isMM)
			try {
				runMM_MStep(emIteraction);
			} catch (UnsupportedMethodException e) {
				System.out.println("MM updates are not supported in "
						+ this.getClass().getName());
				e.printStackTrace();
			}
		else
			try {
				runGradient_MStep(emIteraction);
			} catch (UnsupportedMethodException e) {
				System.out.println("Gradient updates are not supported in "
						+ this.getClass().getName());
				e.printStackTrace();
			}

		/* Update hyper-parameters of model parameters */
		updateModelHyperParameters();

	}

	/* Use MM updates to increase the likelihood lower bound in M step */
	protected abstract void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException;

	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected abstract void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException;

	/* Update hyper-parameters of model parameters */
	protected abstract void updateModelHyperParameters();

	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		for (int k = 0; k < numOfRowGroups; k++)
			if (Math.abs((alpha[k] - prevAlpha[k]) / alpha[k]) > relativeParameterPrecision)
				return true;

		for (int l = 0; l < numOfColumnGroups; l++)
			if (Math.abs((beta[l] - prevBeta[l]) / beta[l]) > relativeParameterPrecision)
				return true;

		System.out.print("alpha: ");
		ArrayMethods.printDoubleArray(alpha, " ", System.out);
		System.out.print("prevAlpha: ");
		ArrayMethods.printDoubleArray(prevAlpha, " ", System.out);
		System.out.print("beta: ");
		ArrayMethods.printDoubleArray(beta, " ", System.out);
		System.out.print("prevBeta: ");
		ArrayMethods.printDoubleArray(prevBeta, " ", System.out);

		System.out.println("relativeParameterPrecision: "
				+ relativeParameterPrecision);

		System.out.println("Alpha va Beta did not pass the test!");

		return false;
	}

	/* Print model parameters and the likelihood lower bound */
	public abstract void printModel(PrintStream outputStream,
			String formatString);

	/* Read the sparse data */
	protected void readEdgeList(BufferedReader dataReader, String dataDelim,
			int numOfEdges, LinkedList<MatrixElement> edges) {
		try {
			for (int e = 0; e < numOfEdges; e++) {
				StringTokenizer tokenizer = new StringTokenizer(dataReader
						.readLine().trim(), dataDelim);
				int row = Integer.parseInt(tokenizer.nextToken().trim());
				int col = Integer.parseInt(tokenizer.nextToken().trim());
				edges.add(new MatrixElement(row, col, 1.0));
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/* Read the matrix data */
	protected void readMatrix(BufferedReader dataReader, String dataDelim,
			double[][] edges) {
		try {
			for (int i = 0; i < edges.length; i++) {
				StringTokenizer tokenizer = new StringTokenizer(dataReader
						.readLine().trim(), dataDelim);
				for (int j = 0; j < edges[0].length; j++) {
					double value = Double.parseDouble(tokenizer.nextToken()
							.trim());
					if (value != 0)
						edges[i][j] = 1.0;
					else
						edges[i][j] = 0.0;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/* Getters and Setters */
	public int getNumOfRowGroups() {
		return numOfRowGroups;
	}

	public void setNumOfRowGroups(int numOfRowGroups) {
		this.numOfRowGroups = numOfRowGroups;

	}

	public int getNumOfColumnGroups() {
		return numOfColumnGroups;
	}

	public void setNumOfColumnGroups(int numOfColumnGroups) {
		this.numOfColumnGroups = numOfColumnGroups;
	}

	public double[] computeRowMSEs(double[] rowBeta) {
		System.out.println(this.getClass()
				+ " does not support computeRowMSEs(...)");
		System.exit(-1);
		return new double[0];
	}

	public double[] computeMisRates(String outputZFilename,
			String outputWFilename) {
		System.out.println(this.getClass()
				+ " does not support computeMisRates(...)");
		System.exit(-1);
		return new double[0];
	}

	public double computeBIC() {
		System.out.println(this.getClass()
				+ " does not support computeBIC(...)");
		System.exit(-1);
		return Double.NaN;
	}

}
