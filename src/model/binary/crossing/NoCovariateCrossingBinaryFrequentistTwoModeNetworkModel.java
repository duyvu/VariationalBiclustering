package model.binary.crossing;

import java.util.Iterator;

import utils.ArrayMethods;

import exception.UnsupportedMethodException;

import graph.ColumnElement;
import graph.MatrixElement;
import graph.RowElement;
import graph.SparseAccess;

public class NoCovariateCrossingBinaryFrequentistTwoModeNetworkModel extends
		CrossingSparseBinaryFrequentistTwoModeNetworkModel {

	protected double[][] logPi0;
	protected double[][] logPi1;

	/* The lazy initialization is used here */
	public NoCovariateCrossingBinaryFrequentistTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
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
		computeLogPis(theta, logPi0, logPi1);
	}

	protected void computeLogPis(double[][] theta, double[][] logPi0,
			double[][] logPi1) {
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				double sum = theta[k][l];
				double norm = -Math.log(1 + Math.exp(sum));
				logPi0[k][l] = norm;
				logPi1[k][l] = sum + norm;
			}
	}

	/* Use MM updates to increase the likelihood lower bound in M step */
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);

		double[][] pi1 = new double[numOfRowGroups][numOfColumnGroups];
		ArrayMethods.resetDoubleMatrix(pi1);

		// compute the numerators
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			while (it.hasNext()) {
				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();
				for (int k = 0; k < numOfRowGroups; k++)
					for (int l = 0; l < numOfColumnGroups; l++)
						pi1[k][l] += z[i][k] * w[j][l];
			}

		// compute pi and convert to theta
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				pi1[k][l] /= (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]);
				theta[k][l] = Math.log(pi1[k][l] / (1.0 - pi1[k][l]));
			}

	}

	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
		runMM_MStep(emIteraction);
	}

	protected void updateModelHyperParameters() {
		return;
	}

}
