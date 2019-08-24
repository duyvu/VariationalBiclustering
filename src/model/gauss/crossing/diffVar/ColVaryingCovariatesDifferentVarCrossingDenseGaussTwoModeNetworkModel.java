package model.gauss.crossing.diffVar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import umontreal.iro.lecuyer.probdist.NormalDist;
import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel
		extends NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel {

	protected boolean dynamicWeighting = false;
	protected double omega = 1.0;
	protected double[][] nu;

	protected double[][] rowBeta;
	protected double[][] prevRowBeta;

	/* The lazy initialization is used here */
	public ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
		// Create column-varying coefficients for row covariates
		if (numOfRowCovariates > 0) {

			rowBeta = new double[graph.getNumOfColumns()][numOfRowCovariates];
			prevRowBeta = new double[graph.getNumOfColumns()][numOfRowCovariates];

			nu = new double[graph.getNumOfRows()][numOfRowCovariates + 1];
		}
	}

	protected void updateWeights() {
		for (int i = 0; i < graph.getNumOfRows(); i++) {

			double norm = 0;

			for (int p = 0; p < numOfRowCovariates; p++) {
				nu[i][p] = Math.pow(Math.abs(rowCovariates[i][p]) + 1.0, omega);
				norm += nu[i][p];
			}

			nu[i][numOfRowCovariates] = 1.0;
			norm += nu[i][numOfRowCovariates];

			for (int p = 0; p <= numOfRowCovariates; p++)
				nu[i][p] /= norm;
		}
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
	}

	/* Reset model parameters */
	protected void resetModelParamters() {

		super.resetModelParamters();
		ArrayMethods.resetDoubleMatrix(rowBeta);
		ArrayMethods.resetDoubleMatrix(prevRowBeta);

		omega = 1.0;
		updateWeights();

		updateCachedParameters();
	}

	@Override
	/* */
	public void updateCachedParameters() {
		return;
	}

	@Override
	public double getLogProb(double y, int i, int j, int k, int l) {
		double yMean = computeLinearCoefficientCombination(i, j, k, l);
		double prob = NormalDist.density(yMean, Math.sqrt(sigma2[k][l]), y);
		if (prob > 0.0)
			return Math.log(prob);
		else
			return -1000;
	}

	/* Compute the linear combination of coefficients and covariates */
	protected double computeLinearCoefficientCombination(int i, int j, int k,
			int l) {
		double linearSum = theta[k][l];
		for (int p = 0; p < numOfRowCovariates; p++)
			linearSum += rowBeta[j][p] * rowCovariates[i][p];
		return linearSum;
	}

	@Override
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		ArrayMethods.copyDoubleMatrix(sigma2, prevSigma2);
		if (numOfRowCovariates > 0)
			ArrayMethods.copyDoubleMatrix(rowBeta, prevRowBeta);

		if (dynamicWeighting) {
			// When omega goes to 0, all covariates receive the same weights
			omega /= 2.0;
			updateWeights();
		}

		// Theta parameters
		double[][] thetaNum = new double[numOfRowGroups][numOfColumnGroups];
		ArrayMethods.resetDoubleMatrix(thetaNum);
		double[][] thetaDen = new double[numOfRowGroups][numOfColumnGroups];
		ArrayMethods.resetDoubleMatrix(thetaDen);

		// Row parameters
		double[][] rowCovNum = null;
		double[][] rowCovDen = null;
		if (numOfRowCovariates > 0) {
			rowCovNum = new double[graph.getNumOfColumns()][numOfRowCovariates];
			ArrayMethods.resetDoubleMatrix(rowCovNum);
			rowCovDen = new double[graph.getNumOfColumns()][numOfRowCovariates];
			ArrayMethods.resetDoubleMatrix(rowCovDen);
		}

		/* Compute the sum to estimate mean and coefficient parameters */

		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				for (int j = 0; j < graph.getNumOfColumns(); j++) {
					for (int i = 0; i < graph.getNumOfRows(); i++) {

						double yMean = computeLinearCoefficientCombination(i,
								j, k, l);
						double residual = graph.getElement(i, j) - yMean;

						double thetaWeight = (z[i][k] * w[j][l])
								* (1 / nu[i][numOfRowCovariates]);
						thetaNum[k][l] += thetaWeight * residual;
						thetaDen[k][l] += thetaWeight
								* (1 / nu[i][numOfRowCovariates]);

						for (int p = 0; p < numOfRowCovariates; p++) {
							double rowCovWeight = (z[i][k] * w[j][l])
									* (rowCovariates[i][p] / nu[i][p]);
							rowCovNum[j][p] += rowCovWeight * residual;
							rowCovDen[j][p] += rowCovWeight
									* (rowCovariates[i][p] / nu[i][p]);
						}
					}
				}
			}
		}

		/* Estimate mean parameters */
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				theta[k][l] = prevTheta[k][l] + thetaNum[k][l] / thetaDen[k][l];

		/* Estimate coefficient parameters */
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int p = 0; p < numOfRowCovariates; p++)
				rowBeta[j][p] = prevRowBeta[j][p] + rowCovNum[j][p]
						/ rowCovDen[j][p];

		/* Estimate variance parameters */
		estimateVarianceParameters();

	}

	protected void estimateVarianceParameters()
			throws UnsupportedMethodException {

		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				sigma2[k][l] = 0;
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					for (int j = 0; j < graph.getNumOfColumns(); j++) {
						double yMean = computeLinearCoefficientCombination(i,
								j, k, l);
						sigma2[k][l] += z[i][k] * w[j][l]
								* Math.pow(graph.getElement(i, j) - yMean, 2);
						denominator += z[i][k] * w[j][l];
					}
				}
				sigma2[k][l] /= denominator;
			}
		}

	}

	@Override
	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
	}

	@Override
	protected void updateModelHyperParameters() {
		return;
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		super.printModel(outputStream, formatString);
		if (numOfRowCovariates > 0) {
			System.out.println("rowBeta: ");
			ArrayMethods.printDoubleMatrix(rowBeta, " ", System.out);
		}
	}

	@Override
	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int p = 0; p < numOfRowCovariates; p++)
				if (Math.abs((rowBeta[j][p] - prevRowBeta[j][p])
						/ rowBeta[j][p]) > relativeParameterPrecision)
					return true;

		return false;
	}

	/* Getters and Setters */

	public boolean isDynamicWeighting() {
		return dynamicWeighting;
	}

	public void setDynamicWeighting(boolean dynamicWeighting) {
		this.dynamicWeighting = dynamicWeighting;
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

			double[] rowTheta = new double[getNumOfRowGroups()];
			ArrayMethods.resetDoubleArray(rowTheta);
			for (int k = 0; k < getNumOfRowGroups(); k++)
				for (int l = 0; l < getNumOfColumnGroups(); l++)
					rowTheta[k] += theta[k][l];
			int[] rowLabels = new int[2];
			for (int k = 0; k < 2; k++) {
				if (rowTheta[k] == Math.max(rowTheta[0], rowTheta[1]))
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

			double[] colTheta = new double[getNumOfColumnGroups()];
			ArrayMethods.resetDoubleArray(colTheta);
			for (int l = 0; l < getNumOfColumnGroups(); l++)
				for (int k = 0; k < getNumOfRowGroups(); k++)
					colTheta[l] += theta[k][l];
			int[] colLabels = new int[3];
			for (int l = 0; l < 3; l++) {
				if (colTheta[l] == Math.max(colTheta[0],
						Math.max(colTheta[1], colTheta[2])))
					colLabels[l] = 2;
				else if (colTheta[l] == Math.min(colTheta[0],
						Math.min(colTheta[1], colTheta[2])))
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
	 * -2 * LB + ((K-1) + (L-1) + 2 * K * L + numRowCovs*numCols +
	 * numColCovs*numRows)*log(numRows*numCols)
	 */
	public double computeBIC() {

		double BIC = -2 * getLikelihoodLowerBound();

		double mixingParameters = numOfRowGroups - 1.0 + numOfColumnGroups
				- 1.0;
		double modelParamters = 2 * numOfRowGroups * numOfColumnGroups;
		double covariateParameters = numOfRowCovariates
				* graph.getNumOfColumns();
		BIC += (mixingParameters + modelParamters + covariateParameters)
				* Math.log(graph.getNumOfRows() * graph.getNumOfColumns());

		return BIC;
	}

}
