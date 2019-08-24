package model.binary.additive;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.TreeMap;

import org.apache.commons.math.analysis.DifferentiableMultivariateRealFunction;
import org.apache.commons.math.analysis.MultivariateRealFunction;
import org.apache.commons.math.analysis.MultivariateVectorialFunction;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.SimpleScalarValueChecker;
import org.apache.commons.math.optimization.general.NonLinearConjugateGradientOptimizer;

import utils.ArrayMethods;

import exception.UnsupportedMethodException;

import graph.ColumnElement;
import graph.MatrixElement;
import graph.RowElement;
import graph.SparseAccess;

public class NoCovariateAdditiveBinaryTwoModeNetworkModel extends
		AdditiveSparseBinaryTwoModeNetworkModel implements
		DifferentiableMultivariateRealFunction, MultivariateRealFunction {

	protected int[] rowDegrees;
	protected int[] columnDegrees;

	protected double[][] logPi0;
	protected double[][] logPi1;

	/* The lazy initialization is used here */
	public NoCovariateAdditiveBinaryTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
		computeDegrees();
	}

	/*
	 * This function should be called in the constructor so that degree
	 * statistics are only computed one time
	 */
	protected void computeDegrees() {
		rowDegrees = new int[graph.getNumOfRows()];
		ArrayMethods.resetIntegerArray(rowDegrees);
		columnDegrees = new int[graph.getNumOfColumns()];
		ArrayMethods.resetIntegerArray(columnDegrees);
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			for (; it.hasNext();) {
				MatrixElement element = it.next();
				rowDegrees[element.getRow()]++;
				columnDegrees[element.getCol()]++;
			}
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		logPi0 = new double[numOfRowGroups][numOfColumnGroups];
		logPi1 = new double[numOfRowGroups][numOfColumnGroups];
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		updateCachedParameters();
	}

	public double getLogProb(int y, int i, int j, int k, int l) {
		return (y == 1) ? logPi1[k][l] : logPi0[k][l];
	}

	protected double computeSumOfLogProbs() {
		return computeSumOfLogProbs(logPi0, logPi1);
	}

	protected double computeSumOfLogProbs(double[][] logPi0, double[][] logPi1) {
		double sumOfLogProbs = 0;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {

				double sumOfZW_Y1 = 0;
				Iterator<MatrixElement> it = ((SparseAccess) graph)
						.elementIterator();
				if (it != null)
					for (; it.hasNext();) {
						MatrixElement element = it.next();
						int i = element.getRow();
						int j = element.getCol();
						sumOfZW_Y1 += z[i][k] * w[j][l];
					}

				sumOfLogProbs += sumOfZW_Y1 * logPi1[k][l];
				sumOfLogProbs += ((graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]) - sumOfZW_Y1)
						* logPi0[k][l];
			}

		return sumOfLogProbs;
	}

	/* This function must use prevW when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		for (int l = 0; l < numOfColumnGroups; l++) {

			double sumOfW_Y1 = 0;
			Iterator<ColumnElement> it = ((SparseAccess) graph).rowIterator(i);
			if (it != null)
				for (; it.hasNext();) {
					ColumnElement element = it.next();
					int j = element.getCol();
					sumOfW_Y1 += prevW[j][l];
				}

			rowSumOfLogProbs += sumOfW_Y1 * logPi1[k][l];
			rowSumOfLogProbs += (graph.getNumOfColumns() * beta[l] - sumOfW_Y1)
					* logPi0[k][l];
		}

		return rowSumOfLogProbs;
	}

	/* This function must use prevZ when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		for (int k = 0; k < numOfRowGroups; k++) {

			double sumOfZ_Y1 = 0;
			Iterator<RowElement> it = ((SparseAccess) graph).columnIterator(j);
			if (it != null)
				for (; it.hasNext();) {
					RowElement element = it.next();
					int i = element.getRow();
					sumOfZ_Y1 += prevZ[i][k];
				}

			columnSumOfLogProbs += sumOfZ_Y1 * logPi1[k][l];
			columnSumOfLogProbs += (graph.getNumOfRows() * alpha[k] - sumOfZ_Y1)
					* logPi0[k][l];
		}

		return columnSumOfLogProbs;
	}

	public void updateCachedParameters() {
		computeLogPis(theta, gamma, logPi0, logPi1);
	}

	protected void computeLogPis(double[] theta, double[] gamma,
			double[][] logPi0, double[][] logPi1) {
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				double sum = theta[k] + gamma[l];
				double norm = -Math.log(1 + Math.exp(sum));
				logPi0[k][l] = norm;
				logPi1[k][l] = sum + norm;
			}
	}

	/* Use MM updates to increase the likelihood lower bound in M step */
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleArray(theta, prevTheta);
		ArrayMethods.copyDoubleArray(gamma, prevGamma);

		double[][] prevLambda = new double[numOfRowGroups][numOfColumnGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				prevLambda[k][l] = 1 + Math.exp(prevTheta[k] + prevGamma[l]);

		// numOfRowGroups >= 1 AND numOfColumnGroups >= 1

		// Estimate theta
		// To make the model identifiable, theta[0] = 0 while other parameters
		// are estimated.
		theta[0] = 0;
		// If numOfRowGroups = 1, we only need to estimate gamma
		if (numOfRowGroups > 1) {
			for (int k = 1; k < numOfRowGroups; k++) {
				// Compute the numerator
				double numerator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++)
					numerator += z[i][k] * rowDegrees[i];
				// Compute the denominator
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					double term = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						term += (Math.exp(prevGamma[l]) / prevLambda[k][l])
								* (graph.getNumOfColumns() * beta[l]);
					denominator += (z[i][k] / Math.exp(prevTheta[k])) * term;
				}
				theta[k] = .5 * Math.log(numerator / denominator);
			}
		}

		// Estimate gamma
		for (int l = 0; l < numOfColumnGroups; l++) {
			// Compute the numerator
			double numerator = 0;
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				numerator += w[j][l] * columnDegrees[j];
			// Compute the denominator
			double denominator = 0;
			for (int j = 0; j < graph.getNumOfColumns(); j++) {
				double term = 0;
				for (int k = 0; k < numOfRowGroups; k++) {
					term += (Math.exp(prevTheta[k]) / prevLambda[k][l])
							* (graph.getNumOfRows() * alpha[k]);
				}
				denominator += (w[j][l] / Math.exp(prevGamma[l])) * term;
			}
			gamma[l] = .5 * Math.log(numerator / denominator);
		}

	}

	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleArray(theta, prevTheta);
		ArrayMethods.copyDoubleArray(gamma, prevGamma);

		double[] params = new double[numOfRowGroups - 1 + numOfColumnGroups];
		packParameters(params, theta, gamma);

		NonLinearConjugateGradientOptimizer optimizer = new NonLinearConjugateGradientOptimizer(
				gradientFormula);
		optimizer.setMaxIterations(maxIterations);
		optimizer.setMaxEvaluations(maxEvaluations);
		optimizer.setConvergenceChecker(new SimpleScalarValueChecker(
				relativeThreshold, absoluteThreshold));

		RealPointValuePair optimum;
		try {
			optimum = optimizer.optimize(this, GoalType.MAXIMIZE, params);
			params = optimum.getPoint();
			unpackParameters(params, theta, gamma);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public MultivariateVectorialFunction gradient() {
		return new MultivariateVectorialFunction() {
			public double[] value(double[] beta) {
				try {
					return computeGradient(beta);
				} catch (Exception e) {
					e.printStackTrace();
					return null;
				}
			}
		};
	}

	public MultivariateRealFunction partialDerivative(final int k) {
		return new MultivariateRealFunction() {
			public double value(double[] beta) {
				try {
					return computeGradient(beta)[k];
				} catch (Exception e) {
					e.printStackTrace();
					return Double.NaN;
				}
			}
		};
	}

	public double value(double[] params) {

		double[] myTheta = new double[numOfRowGroups];
		double[] myGamma = new double[numOfColumnGroups];
		unpackParameters(params, myTheta, myGamma);

		double[][] myLogPi0 = new double[numOfRowGroups][numOfColumnGroups];
		double[][] myLogPi1 = new double[numOfRowGroups][numOfColumnGroups];
		computeLogPis(myTheta, myGamma, myLogPi0, myLogPi1);

		return computeSumOfLogProbs(myLogPi0, myLogPi1);
	}

	public double[] computeGradient(double[] params) {
		double[] grad = new double[numOfRowGroups - 1 + numOfColumnGroups];
		for (int r = 0; r < grad.length; r++)
			grad[r] = 0;

		double[] myTheta = new double[numOfRowGroups];
		double[] myGamma = new double[numOfColumnGroups];
		unpackParameters(params, myTheta, myGamma);

		/* The sum over non-zero edges of the grad */
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			for (; it.hasNext();) {
				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();
				/* For theta */
				for (int k = 1; k < numOfRowGroups; k++) {
					double temp = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						temp += w[j][l];
					grad[k - 1] += z[i][k] * temp;
				}
				/* For gamma */
				for (int l = 0; l < numOfColumnGroups; l++) {
					double temp = 0;
					for (int k = 0; k < numOfRowGroups; k++)
						temp += w[i][k];
					grad[numOfRowGroups - 1 + l] += w[j][l] * temp;
				}
			}

		double[][] norm = new double[numOfRowGroups][numOfColumnGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				norm[k][l] = Math.exp(theta[k] + gamma[l]);
				norm[k][l] = norm[k][l] / (1 + norm[k][l]);
			}

		/* The norm of the grad for theta */
		for (int k = 1; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				grad[k - 1] -= (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]) * norm[k][l];
		/* The norm of the grad for gamma */
		for (int l = 0; l < numOfColumnGroups; l++)
			for (int k = 0; k < numOfRowGroups; k++)
				grad[numOfRowGroups - 1 + l] -= (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]) * norm[k][l];

		return grad;
	}

	/*
	 * Pack model parameters to params = K - 1 theta parameters + L gamma
	 * parameters
	 */
	protected void packParameters(double[] params, double[] myTheta,
			double[] myGamma) {
		for (int k = 1; k < numOfRowGroups; k++)
			params[k - 1] = myTheta[k];
		for (int l = 0; l < numOfColumnGroups; l++)
			params[numOfRowGroups - 1 + l] = myGamma[l];
	}

	/*
	 * Unpack params = K - 1 theta parameters + L gamma parameters to model
	 * parameters
	 */
	protected void unpackParameters(double[] params, double[] myTheta,
			double[] myGamma) {
		myTheta[0] = 0;
		for (int k = 1; k < numOfRowGroups; k++)
			myTheta[k] = params[k - 1];
		for (int l = 0; l < numOfColumnGroups; l++)
			myGamma[l] = params[numOfRowGroups - 1 + l];
	}

	protected void updateModelHyperParameters() {
		return;
	}

	public double[] computeMisRates(String outputZFilename,
			String outputWFilename) {

		double[] misRates = new double[3];

		try {
			BufferedReader rowReader = new BufferedReader(new FileReader(
					new File(outputZFilename)));
			int numRows = Integer.parseInt(rowReader.readLine());
			int[] trueZs = new int[numRows];
			for (int i = 0; i < numRows; i++)
				trueZs[i] = Integer.parseInt(rowReader.readLine());
			rowReader.close();
			// ArrayMethods.printIntegerArray(trueZs, "\t", System.out);

			TreeMap<Double, Integer> sortedThetas = new TreeMap<Double, Integer>();
			for (int k = 0; k < getNumOfRowGroups(); k++)
				sortedThetas.put(theta[k], k);
			int[] rowLabels = new int[getNumOfRowGroups()];
			int rowRank = 0;
			for (int k : sortedThetas.values())
				rowLabels[k] = rowRank++;
			ArrayMethods.printIntegerArray(rowLabels, "\t", System.out);

			BufferedReader colReader = new BufferedReader(new FileReader(
					new File(outputWFilename)));
			int numCols = Integer.parseInt(colReader.readLine());
			int[] trueWs = new int[numCols];
			for (int j = 0; j < numCols; j++)
				trueWs[j] = Integer.parseInt(colReader.readLine());
			colReader.close();
			// ArrayMethods.printIntegerArray(trueWs, "\t", System.out);

			TreeMap<Double, Integer> sortedGammas = new TreeMap<Double, Integer>();
			for (int l = 0; l < getNumOfColumnGroups(); l++)
				sortedGammas.put(gamma[l], l);
			int[] colLabels = new int[getNumOfColumnGroups()];
			int colRank = 0;
			for (int l : sortedGammas.values())
				colLabels[l] = colRank++;
			ArrayMethods.printIntegerArray(colLabels, "\t", System.out);

			double rowRate = 0.0;
			for (int i = 0; i < graph.getNumOfRows(); i++) {
				int myLabel = rowLabels[whichMax(z[i])];
				if (myLabel != trueZs[i])
					rowRate += 1.0;
			}
			rowRate /= graph.getNumOfRows();

			// w[j][l]

			double colRate = 0.0;
			for (int j = 0; j < graph.getNumOfColumns(); j++) {
				int myLabel = colLabels[whichMax(w[j])];
				if (myLabel != trueWs[j])
					colRate += 1.0;
			}
			colRate /= graph.getNumOfColumns();

			double totalRate = (rowRate * graph.getNumOfRows() + colRate
					* graph.getNumOfColumns())
					/ (graph.getNumOfRows() + graph.getNumOfColumns());

			misRates[0] = rowRate;
			misRates[1] = colRate;
			misRates[2] = totalRate;

		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return misRates;
	}

	public static int whichMax(double[] values) {
		if (values.length < 1)
			return -1;

		int maxIndex = 0;
		double maxValue = values[0];
		for (int k = 1; k < values.length; k++) {
			if (maxValue < values[k]) {
				maxValue = values[k];
				maxIndex = k;
			}
		}
		return maxIndex;

	}

	/*
	 * -2 * LB + ((K-1) + (L-1) + K + L + numRowCovs*numCols +
	 * numColCovs*numRows)*log(numRows*numCols)
	 */
	public double computeBIC() {
		
		double BIC = -2 * getLikelihoodLowerBound();
		
		double mixingParameters = numOfRowGroups - 1.0 + numOfColumnGroups
				- 1.0;
		double modelParamters = numOfRowGroups + numOfColumnGroups;
		double covariateParameters = 0.0;
		BIC += (mixingParameters + modelParamters + covariateParameters)
				* Math.log(graph.getNumOfRows() * graph.getNumOfColumns());
		
		return BIC;
	}
	
}
