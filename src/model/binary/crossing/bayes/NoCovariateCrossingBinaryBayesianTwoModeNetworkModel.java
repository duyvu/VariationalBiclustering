package model.binary.crossing.bayes;

import java.util.Iterator;

import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;

import utils.ArrayMethods;

import exception.UnsupportedMethodException;

import graph.ColumnElement;
import graph.MatrixElement;
import graph.RowElement;
import graph.SparseAccess;

public class NoCovariateCrossingBinaryBayesianTwoModeNetworkModel extends
		CrossingSparseBinaryBayesianTwoModeNetworkModel {

	protected double[][] digamma_alpha_theta;
	protected double[][] digamma_beta_theta;
	protected double[][] digamma_alpha_beta_theta;

	/* The lazy initialization is used here */
	public NoCovariateCrossingBinaryBayesianTwoModeNetworkModel() {
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
		digamma_alpha_theta = new double[numOfRowGroups][numOfColumnGroups];
		digamma_beta_theta = new double[numOfRowGroups][numOfColumnGroups];
		digamma_alpha_beta_theta = new double[numOfRowGroups][numOfColumnGroups];
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		updateCachedParameters();
	}

	public double getLogProb(int y, int i, int j, int k, int l) {
		return (y == 1) ? (digamma_alpha_theta[k][l] - digamma_alpha_beta_theta[k][l])
				: (digamma_beta_theta[k][l] - digamma_alpha_beta_theta[k][l]);
	}

	@Override
	protected double sumModelPriorKL() {
		double sum = 0.0;
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				sum += (alpha_theta[k][l] - alpha_theta_0[k][l])
						* digamma_alpha_theta[k][l];
				sum += (beta_theta[k][l] - beta_theta_0[k][l])
						* digamma_beta_theta[k][l];
				sum += (alpha_theta_0[k][l] - alpha_theta[k][l]
						+ beta_theta_0[k][l] - beta_theta[k][l])
						* digamma_alpha_beta_theta[k][l];
				sum += Beta.logBeta(alpha_theta_0[k][l], beta_theta_0[k][l])
						- Beta.logBeta(alpha_theta[k][l], beta_theta[k][l]);
			}
		return sum;
	}

	@Override
	protected double computeSumOfLogProbs() {
		return computeSumOfLogProbs(digamma_alpha_theta, digamma_beta_theta);
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

				sumOfLogProbs += sumOfZW_Y1
						* (digamma_alpha_theta[k][l] - digamma_beta_theta[k][l]);
				sumOfLogProbs += (alpha[k] - lambda_alpha_0[k])
						* (beta[l] - lambda_beta_0[l])
						* (digamma_beta_theta[k][l] - digamma_alpha_beta_theta[k][l]);
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

			rowSumOfLogProbs += sumOfW_Y1
					* (digamma_alpha_theta[k][l] - digamma_beta_theta[k][l]);
			rowSumOfLogProbs += (beta[l] - lambda_beta_0[l])
					* (digamma_beta_theta[k][l] - digamma_alpha_beta_theta[k][l]);
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

			columnSumOfLogProbs += sumOfZ_Y1
					* (digamma_alpha_theta[k][l] - digamma_beta_theta[k][l]);
			columnSumOfLogProbs += (alpha[k] - lambda_alpha_0[k])
					* (digamma_beta_theta[k][l] - digamma_alpha_beta_theta[k][l]);
		}

		return columnSumOfLogProbs;
	}

	public void updateCachedParameters() {
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				digamma_alpha_theta[k][l] = Gamma.digamma(alpha_theta[k][l]);
				digamma_beta_theta[k][l] = Gamma.digamma(beta_theta[k][l]);
				digamma_alpha_beta_theta[k][l] = Gamma
						.digamma(alpha_theta[k][l] + beta_theta[k][l]);
			}
	}

	/* Use MM updates to increase the likelihood lower bound in M step */
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(alpha_theta, prev_alpha_theta);
		ArrayMethods.copyDoubleMatrix(beta_theta, prev_beta_theta);

		ArrayMethods.copyDoubleMatrix(alpha_theta_0, alpha_theta);
		ArrayMethods.copyDoubleMatrix(beta_theta_0, beta_theta);

		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			while (it.hasNext()) {
				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();
				for (int k = 0; k < numOfRowGroups; k++)
					for (int l = 0; l < numOfColumnGroups; l++) {
						alpha_theta[k][l] += z[i][k] * w[j][l];
						beta_theta[k][l] -= z[i][k] * w[j][l];
					}
			}

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				beta_theta[k][l] += (alpha[k] - lambda_alpha_0[k])
						* (beta[l] - lambda_beta_0[l]);

	}

	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
		runMM_MStep(emIteraction);
	}

	/* We don't have hyper parameters to optimize */
	protected void updateModelHyperParameters() {
		return;
	}

}
