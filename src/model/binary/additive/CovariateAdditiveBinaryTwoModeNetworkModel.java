package model.binary.additive;

import java.io.PrintStream;
import java.util.Iterator;

import utils.ArrayMethods;
import exception.UnsupportedMethodException;
import graph.MatrixElement;
import graph.SparseAccess;

public class CovariateAdditiveBinaryTwoModeNetworkModel extends
		AdditiveSparseBinaryTwoModeNetworkModel {

	public double[] getRowBeta() {
		return rowBeta;
	}

	public double[] getColumnBeta() {
		return columnBeta;
	}

	public double[] getEdgeBeta() {
		return edgeBeta;
	}

	// This parameter and maxEvaluations are used for backtracking
	protected boolean usingBackTracking = true;

	protected double[] rowBeta;
	protected double[] prevRowBeta;
	protected double[] columnBeta;
	protected double[] prevColumnBeta;
	protected double[] edgeBeta;
	protected double[] prevEdgeBeta;

	/* The lazy initialization is used here */
	public CovariateAdditiveBinaryTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
		// Create row coefficients
		if (numOfRowCovariates > 0) {
			rowBeta = new double[numOfRowCovariates];
			prevRowBeta = new double[numOfRowCovariates];
		}

		// Create column coefficients
		if (numOfColCovariates > 0) {
			columnBeta = new double[numOfColCovariates];
			prevColumnBeta = new double[numOfColCovariates];
		}

		// Create edge coefficients
		if (numOfEdgeCovariates > 0) {
			edgeBeta = new double[numOfEdgeCovariates];
			prevEdgeBeta = new double[numOfEdgeCovariates];
		}

	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
	}

	public double getLogProb(int y, int i, int j, int k, int l) {

		double linSum = computeLinearCoefficientCombination(i, j, k, l);
		double expLinSum = Math.exp(linSum);

		// The zero part of the probability
		double logProb = -1.0 * Math.log(1 + expLinSum);

		// The non-zero part of the probability
		if (y == 1)
			logProb += linSum;

		return logProb;
	}

	/* Compute the linear combination of coefficients and covariates */
	protected double computeLinearCoefficientCombination(int i, int j, int k,
			int l) {
		double linearSum = theta[k] + gamma[l];
		for (int p = 0; p < numOfRowCovariates; p++)
			linearSum += rowBeta[p] * rowCovariates[i][p];
		for (int p = 0; p < numOfColCovariates; p++)
			linearSum += columnBeta[p] * columnCovariates[j][p];
		for (int p = 0; p < numOfEdgeCovariates; p++)
			linearSum += edgeBeta[p] * edgeCovariates[i][j][p];
		return linearSum;
	}

	protected double computeSumOfLogProbs() {
		double sumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					double sumW = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++)
							sumW += w[j][l]
									* getLogProb((int) graph.getElement(i, j),
											i, j, k, l);
					sumOfLogProbs += z[i][k] * sumW;
				}
		} catch (Exception e) {
			System.out.println("Errors in getLogProb(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return sumOfLogProbs;
	}

	/* This function must use prevW when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		try {
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					rowSumOfLogProbs += prevW[j][l]
							* getLogProb((int) graph.getElement(i, j), i, j, k,
									l);
		} catch (Exception e) {
			System.out.println("Errors in  computeRowSumOfLogProbs(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return rowSumOfLogProbs;
	}

	/* This function must use prevZ when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					columnSumOfLogProbs += prevZ[i][k]
							* getLogProb((int) graph.getElement(i, j), i, j, k,
									l);
		} catch (Exception e) {
			System.out.println("Errors in computeColumnSumOfLogProbs(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return columnSumOfLogProbs;
	}

	public void updateCachedParameters() {
		// This dense version does not cache any parameter.
		return;
	}

	/* Use MM updates to increase the likelihood lower bound in M step */
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleArray(theta, prevTheta);
		ArrayMethods.copyDoubleArray(gamma, prevGamma);
		if (numOfRowCovariates > 0)
			ArrayMethods.copyDoubleArray(rowBeta, prevRowBeta);
		if (numOfColCovariates > 0)
			ArrayMethods.copyDoubleArray(columnBeta, prevColumnBeta);
		if (numOfEdgeCovariates > 0)
			ArrayMethods.copyDoubleArray(edgeBeta, prevEdgeBeta);

		// Theta parameters
		double[] thetaGradient = new double[theta.length];
		ArrayMethods.resetDoubleArray(thetaGradient);
		double[] thetaHessian = new double[theta.length];
		ArrayMethods.resetDoubleArray(thetaHessian);

		// GammprevEdgeBetaa parameters
		double[] gammaGradient = new double[gamma.length];
		ArrayMethods.resetDoubleArray(gammaGradient);
		double[] gammaHessian = new double[gamma.length];
		ArrayMethods.resetDoubleArray(gammaHessian);

		// Row parameters
		double[] rowGradient = null;
		double[] rowHessian = null;
		if (numOfRowCovariates > 0) {
			rowGradient = new double[rowBeta.length];
			ArrayMethods.resetDoubleArray(rowGradient);
			rowHessian = new double[rowBeta.length];
			ArrayMethods.resetDoubleArray(rowHessian);
		}

		// Column parameters
		double[] columnGradient = null;
		double[] columnHessian = null;
		if (numOfColCovariates > 0) {
			columnGradient = new double[columnBeta.length];
			ArrayMethods.resetDoubleArray(columnGradient);
			columnHessian = new double[columnBeta.length];
			ArrayMethods.resetDoubleArray(columnHessian);
		}

		// Edge parameters
		double[] edgeGradient = null;
		double[] edgeHessian = null;
		if (numOfEdgeCovariates > 0) {
			edgeGradient = new double[edgeBeta.length];
			ArrayMethods.resetDoubleArray(edgeGradient);
			edgeHessian = new double[edgeBeta.length];
			ArrayMethods.resetDoubleArray(edgeHessian);
		}

		// For the non-zero part: update gradients
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			for (; it.hasNext();) {

				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();

				// Update theta parameters
				for (int k = 0; k < numOfRowGroups; k++) {
					double sum = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						sum += w[j][l];
					thetaGradient[k] += z[i][k] * sum;
				}

				// Update gamma parameters
				for (int l = 0; l < numOfColumnGroups; l++) {
					double sum = 0;
					for (int k = 0; k < numOfRowGroups; k++)
						sum += z[i][k];
					gammaGradient[l] += w[j][l] * sum;
				}

				// Update covariates. Here we assume that covariates are not
				// dependent on group structures so that them can be factorized
				// out
				// of the loops over K and L
				double sumZxW = 0;
				for (int k = 0; k < numOfRowGroups; k++) {
					double sumW = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						sumW += w[j][l];
					sumZxW += z[i][k] * sumW;
				}

				// Update row covariates
				for (int p = 0; p < numOfRowCovariates; p++)
					rowGradient[p] += rowCovariates[i][p] * sumZxW;

				// Update column covariates
				for (int p = 0; p < numOfColCovariates; p++)
					columnGradient[p] += columnCovariates[j][p] * sumZxW;

				// Update edge covariates
				for (int p = 0; p < numOfEdgeCovariates; p++) {
					edgeGradient[p] += edgeCovariates[i][j][p] * sumZxW;
				}

			}

		// For the zero-part: update gradients and Hessian
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int k = 0; k < numOfRowGroups; k++) {
					for (int l = 0; l < numOfColumnGroups; l++) {

						// Compute the common term
						double linComb = computeLinearCoefficientCombination(i,
								j, k, l);
						double expLinComb = Math.exp(linComb);
						double comTerm = z[i][k] * w[j][l] * expLinComb
								/ (1.0 + expLinComb);

						// Update theta covariates
						thetaGradient[k] -= comTerm;
						thetaHessian[k] -= comTerm;

						// Update gamma covariates
						gammaGradient[l] -= comTerm;
						gammaHessian[l] -= comTerm;

						// Update row covariates
						for (int p = 0; p < numOfRowCovariates; p++) {
							rowGradient[p] -= comTerm * rowCovariates[i][p];
							rowHessian[p] -= comTerm
									* Math.pow(rowCovariates[i][p], 2);
						}

						// Update column covariates
						for (int p = 0; p < numOfColCovariates; p++) {
							columnGradient[p] -= comTerm
									* columnCovariates[i][p];
							columnHessian[p] -= comTerm
									* Math.pow(columnCovariates[i][p], 2);
						}

						// Update edge covariates
						for (int p = 0; p < numOfEdgeCovariates; p++) {
							edgeGradient[p] -= comTerm
									* edgeCovariates[i][j][p];
							edgeHessian[p] -= comTerm
									* Math.pow(edgeCovariates[i][j][p], 2);

						}
					}
				}

		// Multiply Hessian by (2 + row + column + edge)
		int dimension = 2 + numOfRowCovariates + numOfColCovariates
				+ numOfEdgeCovariates;
		for (int k = 0; k < numOfRowGroups; k++)
			thetaHessian[k] *= dimension;
		// Update gamma parameters
		for (int l = 0; l < numOfColumnGroups; l++)
			gammaHessian[l] *= dimension;
		// Update row covariates
		for (int p = 0; p < numOfRowCovariates; p++)
			rowHessian[p] *= dimension;
		// Update column covariates
		for (int p = 0; p < numOfColCovariates; p++)
			columnHessian[p] *= dimension;
		// Update edge covariates
		for (int p = 0; p < numOfEdgeCovariates; p++)
			edgeHessian[p] *= dimension;

		/*
		 * Using one-step Newton method as discussed in Chapter 14 of Lange's
		 * book but with backtracking option make sure the likelihood is
		 * increased
		 */
		if (usingBackTracking) {
			double prevBound = computeSumOfLogProbs();
			double scaling = 1.0;
			int step = 1;
			for (; step <= maxEvaluations; step++) {
				computeEstimates(thetaGradient, thetaHessian, gammaGradient,
						gammaHessian, rowGradient, rowHessian, columnGradient,
						columnHessian, edgeGradient, edgeHessian, scaling);
				double newBound = computeSumOfLogProbs();
				if (newBound > prevBound)
					break;
				else
					prevBound = newBound;
			}
			if (step > maxEvaluations) {
				System.out
						.println("We could not improve the bound in M step!!!");
				ArrayMethods.copyDoubleArray(prevTheta, theta);
				ArrayMethods.copyDoubleArray(prevGamma, gamma);
				if (numOfRowCovariates > 0)
					ArrayMethods.copyDoubleArray(prevRowBeta, rowBeta);
				if (numOfColCovariates > 0)
					ArrayMethods.copyDoubleArray(prevColumnBeta, columnBeta);
				if (numOfEdgeCovariates > 0)
					ArrayMethods.copyDoubleArray(prevEdgeBeta, edgeBeta);
			}
		} else {
			/*
			 * Update parameters using one-step Newton's method if backtracking
			 * is not used.
			 */
			computeEstimates(thetaGradient, thetaHessian, gammaGradient,
					gammaHessian, rowGradient, rowHessian, columnGradient,
					columnHessian, edgeGradient, edgeHessian, 1.0);
		}

	}

	protected void computeEstimates(double[] thetaGradient,
			double[] thetaHessian, double[] gammaGradient,
			double[] gammaHessian, double[] rowGradient, double[] rowHessian,
			double[] columnGradient, double[] columnHessian,
			double[] edgeGradient, double[] edgeHessian, double scaling) {
		// Update theta parameters
		theta[0] = 0;
		// If numOfRowGroups = 1, we only need to estimate gamma
		if (numOfRowGroups > 1)
			for (int k = 1; k < numOfRowGroups; k++)
				theta[k] = prevTheta[k] - scaling * thetaGradient[k]
						/ thetaHessian[k];

		// Update gamma parameters
		for (int l = 0; l < numOfColumnGroups; l++)
			gamma[l] = prevGamma[l] - scaling * gammaGradient[l]
					/ gammaHessian[l];

		// Update row covariates
		for (int p = 0; p < numOfRowCovariates; p++)
			rowBeta[p] = prevRowBeta[p] - scaling * rowGradient[p]
					/ rowHessian[p];

		// Update column covariates
		for (int p = 0; p < numOfColCovariates; p++)
			columnBeta[p] = prevColumnBeta[p] - scaling * columnGradient[p]
					/ columnHessian[p];

		// Update edge covariates
		for (int p = 0; p < numOfEdgeCovariates; p++)
			edgeBeta[p] = prevEdgeBeta[p] - scaling * edgeGradient[p]
					/ edgeHessian[p];
	}

	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {

	}

	protected void updateModelHyperParameters() {
		// This model does not have any hyper-parameters for model parameters.
		// In Bayesian versions, this method must be implemented to update
		// corresponding hyper-parameters.
		return;
	}

	/* Getters and Setters */

	public boolean isUsingBackTracking() {
		return usingBackTracking;
	}

	public void setUsingBackTracking(boolean usingBackTracking) {
		this.usingBackTracking = usingBackTracking;
	}

	/* Reset model paramters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		ArrayMethods.resetDoubleArray(rowBeta);
		ArrayMethods.resetDoubleArray(prevRowBeta);
		ArrayMethods.resetDoubleArray(columnBeta);
		ArrayMethods.resetDoubleArray(prevColumnBeta);
		ArrayMethods.resetDoubleArray(edgeBeta);
		ArrayMethods.resetDoubleArray(prevEdgeBeta);
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		super.printModel(outputStream, formatString);
		if (numOfRowCovariates > 0) {
			System.out.println("rowBeta: ");
			ArrayMethods.printDoubleArray(rowBeta, " ", System.out);
		}
		if (numOfColCovariates > 0) {
			System.out.println("columnBeta: ");
			ArrayMethods.printDoubleArray(columnBeta, " ", System.out);
		}
		if (numOfEdgeCovariates > 0) {
			System.out.println("edgeBeta: ");
			ArrayMethods.printDoubleArray(edgeBeta, " ", System.out);
		}
	}

	@Override
	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int p = 0; p < numOfRowCovariates; p++)
			if (Math.abs((rowBeta[p] - prevRowBeta[p]) / rowBeta[p]) > relativeParameterPrecision)
				return true;

		for (int p = 0; p < numOfColCovariates; p++)
			if (Math.abs((columnBeta[p] - prevColumnBeta[p]) / columnBeta[p]) > relativeParameterPrecision)
				return true;

		for (int p = 0; p < numOfEdgeCovariates; p++)
			if (Math.abs((edgeBeta[p] - prevEdgeBeta[p]) / edgeBeta[p]) > relativeParameterPrecision)
				return true;

		return false;
	}
}
