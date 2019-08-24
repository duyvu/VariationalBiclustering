package model.binary.additive;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;

import utils.ArrayMethods;
import exception.UnsupportedMethodException;
import graph.MatrixElement;
import graph.SparseAccess;

public class RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkModel extends
		AdditiveSparseBinaryTwoModeNetworkModel {

	// If the backtracking option is in use, maxStepHalvingIterations is the
	// number of step-halving iterations will be tried.
	protected boolean usingBackTracking = true;
	protected int maxStepHalvingIterations = 10;

	protected double[][] rowBeta;
	protected double[][] prevRowBeta;
	protected double[][] columnBeta;
	protected double[][] prevColumnBeta;

	/* The lazy initialization is used here */
	public RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
		// Create row coefficients
		if (numOfRowCovariates > 0) {
			rowBeta = new double[graph.getNumOfColumns()][numOfRowCovariates];
			prevRowBeta = new double[graph.getNumOfColumns()][numOfRowCovariates];
		}
		// Create column coefficients
		if (numOfColCovariates > 0) {
			columnBeta = new double[graph.getNumOfRows()][numOfColCovariates];
			prevColumnBeta = new double[graph.getNumOfRows()][numOfColCovariates];
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
			linearSum += rowBeta[j][p] * rowCovariates[i][p];
		for (int p = 0; p < numOfColCovariates; p++)
			linearSum += columnBeta[i][p] * columnCovariates[j][p];
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
			ArrayMethods.copyDoubleMatrix(rowBeta, prevRowBeta);
		if (numOfColCovariates > 0)
			ArrayMethods.copyDoubleMatrix(columnBeta, prevColumnBeta);

		// Theta parameters
		double[] thetaGradient = new double[theta.length];
		ArrayMethods.resetDoubleArray(thetaGradient);
		double[] thetaHessian = new double[theta.length];
		ArrayMethods.resetDoubleArray(thetaHessian);

		// Gamma parameters
		double[] gammaGradient = new double[gamma.length];
		ArrayMethods.resetDoubleArray(gammaGradient);
		double[] gammaHessian = new double[gamma.length];
		ArrayMethods.resetDoubleArray(gammaHessian);

		// Row parameters
		double[][] rowGradient = null;
		double[][] rowHessian = null;
		if (numOfRowCovariates > 0) {
			rowGradient = new double[graph.getNumOfColumns()][rowBeta.length];
			ArrayMethods.resetDoubleMatrix(rowGradient);
			rowHessian = new double[graph.getNumOfColumns()][rowBeta.length];
			ArrayMethods.resetDoubleMatrix(rowHessian);
		}

		// Column parameters
		double[][] columnGradient = null;
		double[][] columnHessian = null;
		if (numOfColCovariates > 0) {
			columnGradient = new double[graph.getNumOfRows()][columnBeta.length];
			ArrayMethods.resetDoubleMatrix(columnGradient);
			columnHessian = new double[graph.getNumOfRows()][columnBeta.length];
			ArrayMethods.resetDoubleMatrix(columnHessian);
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

				// Update covariates.
				// Here we assume that covariates are not dependent on the group
				// structure so that them can be factorized out of the loops
				// over K and L
				double sumZxW = 0;
				for (int k = 0; k < numOfRowGroups; k++) {
					double sumW = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						sumW += w[j][l];
					sumZxW += z[i][k] * sumW;
				}

				// Update row covariates
				for (int p = 0; p < numOfRowCovariates; p++)
					rowGradient[j][p] += rowCovariates[i][p] * sumZxW;

				// Update column covariates
				for (int p = 0; p < numOfColCovariates; p++)
					columnGradient[i][p] += columnCovariates[j][p] * sumZxW;

			}

		double[] thetaTerm = new double[theta.length];
		ArrayMethods.resetDoubleArray(thetaTerm);
		double[] gammaTerm = new double[gamma.length];
		ArrayMethods.resetDoubleArray(gammaTerm);

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
						thetaTerm[k] -= comTerm;

						// Update gamma covariates
						gammaTerm[l] -= comTerm;

						// Update row covariates
						for (int p = 0; p < numOfRowCovariates; p++) {
							rowGradient[j][p] -= comTerm * rowCovariates[i][p];
							rowHessian[j][p] -= comTerm
									* Math.pow(rowCovariates[i][p], 2);
						}

						// Update column covariates
						for (int p = 0; p < numOfColCovariates; p++) {
							columnGradient[i][p] -= comTerm
									* columnCovariates[j][p];
							columnHessian[i][p] -= comTerm
									* Math.pow(columnCovariates[j][p], 2);
						}

					}
				}

		// Multiply Hessian by (2 + row + column + edge)
		int totalCoefficients = 2 + numOfRowCovariates + numOfColCovariates + 0;

		// Update theta parameters
		for (int k = 0; k < numOfRowGroups; k++) {
			thetaGradient[k] += thetaTerm[k];
			thetaHessian[k] = totalCoefficients * thetaTerm[k];
		}

		// Update gamma parameters
		for (int l = 0; l < numOfColumnGroups; l++) {
			gammaGradient[l] += gammaTerm[l];
			gammaHessian[l] = totalCoefficients * gammaTerm[l];
		}

		// Update row covariates
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int p = 0; p < numOfRowCovariates; p++)
				rowHessian[j][p] *= totalCoefficients;

		// Update column covariates
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int p = 0; p < numOfColCovariates; p++)
				columnHessian[i][p] *= totalCoefficients;

		/*
		 * Using one-step Newton method as discussed in Chapter 14 of Kenneth
		 * Lange's book but if the backtracking option is in use, do
		 * step-halving until one of these conditions are satisfied: (a) the
		 * likelihood is increased or (b) maxStepHalvingIterations is reached
		 */
		double scaling = 1.0;
		if (!usingBackTracking) {
			/*
			 * Update parameters using one-step Newton's method with the
			 * step-length 1.0 if the backtracking option is not used.
			 */

			computeEstimates(thetaGradient, thetaHessian, gammaGradient,
					gammaHessian, rowGradient, rowHessian, columnGradient,
					columnHessian, scaling);
		} else {
			/*
			 * Otherwise, we will do backtracking with a shared scaling across
			 * all parameters to save the computation
			 */

			double prevBound = computeSumOfLogProbs();

			computeEstimates(thetaGradient, thetaHessian, gammaGradient,
					gammaHessian, rowGradient, rowHessian, columnGradient,
					columnHessian, scaling);

			double newBound = computeSumOfLogProbs();

			// If new parameters do not improve the bound
			if (newBound <= prevBound) {
				// Do step-halving
				for (int step = 0; step < maxStepHalvingIterations; step++) {
					scaling /= 2.0;
					computeEstimates(thetaGradient, thetaHessian,
							gammaGradient, gammaHessian, rowGradient,
							rowHessian, columnGradient, columnHessian, scaling);
					newBound = computeSumOfLogProbs();
					// Until the bound improved
					if (newBound > prevBound)
						return;
				}

				// Or the maximum number of iterations is reached
				System.out
						.println("We could not improve the bound in M step! Previous estimates are used.");
				ArrayMethods.copyDoubleArray(prevTheta, theta);
				ArrayMethods.copyDoubleArray(prevGamma, gamma);
				if (numOfRowCovariates > 0)
					ArrayMethods.copyDoubleMatrix(prevRowBeta, rowBeta);
				if (numOfColCovariates > 0)
					ArrayMethods.copyDoubleMatrix(prevColumnBeta, columnBeta);
			}
		}

	}

	protected void computeEstimates(double[] thetaGradient,
			double[] thetaHessian, double[] gammaGradient,
			double[] gammaHessian, double[][] rowGradient,
			double[][] rowHessian, double[][] columnGradient,
			double[][] columnHessian, double scaling) {

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
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int p = 0; p < numOfRowCovariates; p++)
				rowBeta[j][p] = prevRowBeta[j][p] - scaling * rowGradient[j][p]
						/ rowHessian[j][p];

		// Update column covariates
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int p = 0; p < numOfColCovariates; p++)
				columnBeta[i][p] = prevColumnBeta[i][p] - scaling
						* columnGradient[i][p] / columnHessian[i][p];

		// System.out.println("z: ");
		// ArrayMethods.printDoubleMatrix(z, " ", 30, System.out);
		//
		// System.out.println("w: ");
		// ArrayMethods.printDoubleMatrix(w, " ", 30, System.out);
		//
		// System.out.println("theta: ");
		// ArrayMethods.printDoubleArray(theta, " ", System.out);
		//
		// System.out.println("gamma: ");
		// ArrayMethods.printDoubleArray(gamma, " ", System.out);
		//
		// if (numOfRowCovariates > 0) {
		// System.out.println("rowBeta: ");
		// ArrayMethods.printDoubleMatrix(rowBeta, " ", 30, System.out);
		// }
		//
		// if (numOfColCovariates > 0) {
		// System.out.println("columnBeta: ");
		// ArrayMethods.printDoubleMatrix(columnBeta, " ", 30, System.out);
		// }

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

	/* Reset model paramters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		ArrayMethods.resetDoubleMatrix(rowBeta);
		ArrayMethods.resetDoubleMatrix(prevRowBeta);
		if (numOfColCovariates > 0) {
			ArrayMethods.resetDoubleMatrix(columnBeta);
			ArrayMethods.resetDoubleMatrix(prevColumnBeta);
		}
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		super.printModel(outputStream, formatString);
		if (numOfRowCovariates > 0) {
			System.out.println("rowBeta: ");
			ArrayMethods.printDoubleMatrix(rowBeta, " ", System.out);
		}
		if (numOfColCovariates > 0) {
			System.out.println("columnBeta: ");
			ArrayMethods.printDoubleMatrix(columnBeta, " ", System.out);
		}
	}

	@Override
	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int p = 0; p < numOfColCovariates; p++)
				if (Math.abs((columnBeta[i][p] - prevColumnBeta[i][p])
						/ columnBeta[i][p]) > relativeParameterPrecision)
					return true;

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int p = 0; p < numOfRowCovariates; p++)
				if (Math.abs((rowBeta[j][p] - prevRowBeta[j][p])
						/ rowBeta[j][p]) > relativeParameterPrecision)
					return true;

		return false;
	}

	/* Getters and Setters */

	public boolean isUsingBackTracking() {
		return usingBackTracking;
	}

	public void setUsingBackTracking(boolean usingBackTracking) {
		this.usingBackTracking = usingBackTracking;
	}

	public int getMaxStepHalvingIterations() {
		return maxStepHalvingIterations;
	}

	public void setMaxStepHalvingIterations(int maxStepHalvingIterations) {
		this.maxStepHalvingIterations = maxStepHalvingIterations;
	}

	public double[] computeRowMSEs(double[] trueRowBeta) {
		double[] MSEs = new double[trueRowBeta.length];
		for (int p = 0; p < trueRowBeta.length; p++) {
			MSEs[p] = 0.0;
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				MSEs[p] += Math.pow(rowBeta[j][p] - trueRowBeta[p], 2);
			MSEs[p] /= graph.getNumOfColumns();
		}
		return MSEs;
	}

	public double[] computeMisRates(String outputZFilename,
			String outputWFilename) {

		double[] misRates = new double[3];

		try {
			if (getNumOfRowGroups() != 2 || getNumOfColumnGroups() != 3) {
				System.out
						.println("This misclassfication evalution only works for 2 x 3");
				System.exit(-1);
			}

			BufferedReader rowReader = new BufferedReader(new FileReader(
					new File(outputZFilename)));
			int numRows = Integer.parseInt(rowReader.readLine());
			int[] trueZs = new int[numRows];
			for (int i = 0; i < numRows; i++)
				trueZs[i] = Integer.parseInt(rowReader.readLine());
			rowReader.close();
			// ArrayMethods.printIntegerArray(trueZs, "\t", System.out);

			int[] rowLabels = new int[2];
			for (int k = 0; k < 2; k++) {
				if (theta[k] == Math.max(theta[0], theta[1]))
					rowLabels[k] = 1;
				else
					rowLabels[k] = 0;
			}
			// ArrayMethods.printIntegerArray(rowLabels, "\t", System.out);

			BufferedReader colReader = new BufferedReader(new FileReader(
					new File(outputWFilename)));
			int numCols = Integer.parseInt(colReader.readLine());
			int[] trueWs = new int[numCols];
			for (int j = 0; j < numCols; j++)
				trueWs[j] = Integer.parseInt(colReader.readLine());
			colReader.close();
			// ArrayMethods.printIntegerArray(trueWs, "\t", System.out);

			int[] colLabels = new int[3];
			for (int l = 0; l < 3; l++) {
				if (gamma[l] == Math
						.max(gamma[0], Math.max(gamma[1], gamma[2])))
					colLabels[l] = 2;
				else if (gamma[l] == Math.min(gamma[0],
						Math.min(gamma[1], gamma[2])))
					colLabels[l] = 0;
				else
					colLabels[l] = 1;
			}
			// ArrayMethods.printIntegerArray(colLabels, "\t", System.out);

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
	 * -2 * LB + ((L-1) + K + L + (K-1) + numRowCovs*numCols +
	 * numColCovs*numRows)*log(numRows*numCols)
	 */
	public double computeBIC() {
		
		double BIC = -2 * getLikelihoodLowerBound();
		
		double mixingParameters = numOfRowGroups - 1.0 + numOfColumnGroups
				- 1.0;
		double modelParamters = numOfRowGroups + numOfColumnGroups;
		double covariateParameters = numOfColCovariates * graph.getNumOfRows()
				+ numOfRowCovariates * graph.getNumOfColumns();
		BIC += (mixingParameters + modelParamters + covariateParameters)
				* Math.log(graph.getNumOfRows() * graph.getNumOfColumns());
		
		return BIC;
	}

}
