package model.gauss.crossing.sameVar;

import java.io.PrintStream;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import umontreal.iro.lecuyer.probdist.NormalDist;
import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class CovariateCrossingDenseGaussTwoModeNetworkModel extends
		CrossingDenseGaussTwoModeNetworkModel {

	public double[] getRowBeta() {
		return rowBeta;
	}

	public double[] getColumnBeta() {
		return columnBeta;
	}

	public double[] getEdgeBeta() {
		return edgeBeta;
	}

	protected double[] rowBeta;
	protected double[] prevRowBeta;
	protected double[] columnBeta;
	protected double[] prevColumnBeta;
	protected double[] edgeBeta;
	protected double[] prevEdgeBeta;

	/* The lazy initialization is used here */
	public CovariateCrossingDenseGaussTwoModeNetworkModel() {
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

	/* Reset model paramters */
	protected void resetModelParamters() {
		super.resetModelParamters();
	}

	@Override
	public double getLogProb(double y, int i, int j, int k, int l) {
		return Math.log(NormalDist.density(computeExpectation(i, j, k, l),
				Math.sqrt(sigma2), y));
	}

	/* Compute the linear combination of coefficients and covariates */
	protected double computeExpectation(int i, int j, int k, int l) {
		double linearSum = theta[k][l];
		for (int p = 0; p < numOfRowCovariates; p++)
			linearSum += rowBeta[p] * rowCovariates[i][p];
		for (int p = 0; p < numOfColCovariates; p++)
			linearSum += columnBeta[p] * columnCovariates[j][p];
		for (int p = 0; p < numOfEdgeCovariates; p++)
			linearSum += edgeBeta[p] * edgeCovariates[i][j][p];
		return linearSum;
	}

	/* Compute the linear combination of coefficients OF covariates */
	protected double computeLinearCoefCovLinearCombination(int i, int j) {
		double linearSum = 0;
		for (int p = 0; p < numOfRowCovariates; p++)
			linearSum += rowBeta[p] * rowCovariates[i][p];
		for (int p = 0; p < numOfColCovariates; p++)
			linearSum += columnBeta[p] * columnCovariates[j][p];
		for (int p = 0; p < numOfEdgeCovariates; p++)
			linearSum += edgeBeta[p] * edgeCovariates[i][j][p];
		return linearSum;
	}

	@Override
	protected double computeSumOfLogProbs() {
		double sumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					double sumW = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++)
							sumW += w[j][l]
									* getLogProb(graph.getElement(i, j), i, j,
											k, l);
					sumOfLogProbs += z[i][k] * sumW;
				}
		} catch (Exception e) {
			System.out.println("Errors in getLogProb(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return sumOfLogProbs;
	}

	@Override
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		try {
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					rowSumOfLogProbs += prevW[j][l]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			System.out.println("Errors in computeRowSumOfLogProbs(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return rowSumOfLogProbs;
	}

	@Override
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					columnSumOfLogProbs += prevZ[i][k]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			System.out.println("Errors in computeColumnSumOfLogProbs(...) "
					+ this.getClass().getName());
			System.exit(-1);
		}

		return columnSumOfLogProbs;
	}

	@Override
	public void updateCachedParameters() {
		// This dense version does not cache any parameter.
		return;
	}

	@Override
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {
		runCoordinate_MStep();
	}

	@Override
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
		runCoordinate_MStep();
	}

	protected void computeCovariates(double[] vector, double[][] matrix, int i,
			int j) {

		// i
		for (int k = 0; k < numOfRowCovariates; k++) {

			vector[k] = rowCovariates[i][k];

			for (int l = 0; l <= k; l++)
				// to itself
				matrix[k][l] = rowCovariates[i][l] * rowCovariates[i][k];

			for (int l = 0; l < numOfColCovariates; l++)
				// j and i
				matrix[numOfRowCovariates + l][k] = columnCovariates[j][l]
						* rowCovariates[i][k];

			for (int l = 0; l < numOfEdgeCovariates; l++)
				// ij and i
				matrix[numOfRowCovariates + numOfColCovariates + l][k] = edgeCovariates[i][j][l]
						* rowCovariates[i][k];

		}

		// j
		for (int k = 0; k < numOfColCovariates; k++) {

			vector[numOfRowCovariates + k] = columnCovariates[j][k];

			for (int l = 0; l <= k; l++)
				// to itself
				matrix[numOfRowCovariates + k][numOfRowCovariates + l] = columnCovariates[j][l]
						* columnCovariates[j][k];
			for (int l = 0; l < numOfEdgeCovariates; l++)
				// ij and j
				matrix[numOfRowCovariates + numOfColCovariates + l][numOfRowCovariates
						+ k] = edgeCovariates[i][j][l] * columnCovariates[j][k];

		}

		// ij
		for (int k = 0; k < numOfEdgeCovariates; k++) {

			vector[numOfRowCovariates + numOfColCovariates + k] = edgeCovariates[i][j][k];

			for (int l = 0; l <= k; l++)
				// to itself
				matrix[numOfRowCovariates + numOfColCovariates + k][numOfRowCovariates
						+ numOfColCovariates + l] = edgeCovariates[i][j][l]
						* edgeCovariates[i][j][k];
		}
	}

	protected void runCoordinate_MStep() {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		prevSigma2 = sigma2;
		if (numOfRowCovariates > 0)
			ArrayMethods.copyDoubleArray(rowBeta, prevRowBeta);
		if (numOfColCovariates > 0)
			ArrayMethods.copyDoubleArray(columnBeta, prevColumnBeta);
		if (numOfEdgeCovariates > 0)
			ArrayMethods.copyDoubleArray(edgeBeta, prevEdgeBeta);

		/* Update THETA */
		try {
			for (int k = 0; k < numOfRowGroups; k++) {
				for (int l = 0; l < numOfColumnGroups; l++) {
					theta[k][l] = 0;
					double denominator = 0;
					for (int i = 0; i < graph.getNumOfRows(); i++) {
						double thetaW = 0;
						double denW = 0;
						for (int j = 0; j < graph.getNumOfColumns(); j++) {
							thetaW += w[j][l]
									* (graph.getElement(i, j) - computeLinearCoefCovLinearCombination(
											i, j));
							denW += w[j][l];
						}
						theta[k][l] += z[i][k] * thetaW;
						denominator += z[i][k] * denW;
					}
					theta[k][l] /= denominator;
				}
			}
		} catch (Exception e) {
			System.out.println("Errors in runCoordinate_MStep(...) ");
			System.exit(-1);
		}

		try {
			/* Update NU */
			int numOfCovariates = numOfRowCovariates + numOfColCovariates
					+ numOfEdgeCovariates;
			double[][] nuMat = new double[numOfCovariates][numOfCovariates];
			ArrayMethods.resetDoubleMatrix(nuMat);
			double[] nuVec = new double[numOfCovariates];
			ArrayMethods.resetDoubleArray(nuVec);
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int j = 0; j < graph.getNumOfColumns(); j++) {
					double[] covVector = new double[numOfCovariates];
					double[][] prodCovMatrix = new double[numOfCovariates][numOfCovariates];
					computeCovariates(covVector, prodCovMatrix, i, j);
					for (int k = 0; k < numOfRowGroups; k++)
						for (int l = 0; l < numOfColumnGroups; l++) {
							double adjustedY_ij = graph.getElement(i, j)
									- theta[k][l];
							for (int p1 = 0; p1 < numOfCovariates; p1++) {
								for (int p2 = 0; p2 <= p1; p2++)
									nuMat[p1][p2] += z[i][k] * w[j][l]
											* prodCovMatrix[p1][p2];
								nuVec[p1] += z[i][k] * w[j][l] * covVector[p1]
										* adjustedY_ij;
							}
						}
				}
			ArrayMethods.reflexDoubleMatrix(nuMat);
			// Compute inverse of nuNum
			RealMatrix iNuMat = new Array2DRowRealMatrix(nuMat);
			DecompositionSolver solver = new LUDecompositionImpl(iNuMat)
					.getSolver();
			iNuMat = solver.getInverse();
			// Compute estimated coefficients (nuMat)^-1 x nuVec
			double[] estimatedCoefficients = iNuMat.operate(nuVec);

			// Copy estimated coefficients into row, column, edge betas
			for (int k = 0; k < numOfRowCovariates; k++)
				rowBeta[k] = estimatedCoefficients[k];
			for (int k = 0; k < numOfColCovariates; k++)
				columnBeta[k] = estimatedCoefficients[numOfRowCovariates + k];
			for (int k = 0; k < numOfEdgeCovariates; k++)
				edgeBeta[k] = estimatedCoefficients[numOfRowCovariates
						+ numOfColCovariates + k];
		} catch (Exception e) {
			System.out.println("Errors in runCoordinate_MStep(...) ");
			System.exit(-1);
		}

		/* Update variance paramter SIGMA2 */
		try {
			sigma2 = 0;
			double denominator = 0;

			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					double sigma2W = 0;
					double denW = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++) {
							sigma2W += w[j][l]
									* Math.pow(graph.getElement(i, j)
											- computeExpectation(i, j, k, l), 2);
							denW += w[j][l];
						}
					sigma2 += z[i][k] * sigma2W;
					denominator += z[i][k] * denW;
				}
			sigma2 /= denominator;
		} catch (Exception e) {
			System.out.println("Errors in runCoordinate_MStep(...) ");
			System.exit(-1);
		}
	}

	@Override
	protected void updateModelHyperParameters() {
		// This model does not have any hyper-parameters for model parameters.
		// In Bayesian versions, this method must be implemented to update
		// corresponding hyper-parameters.
		return;
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
}
