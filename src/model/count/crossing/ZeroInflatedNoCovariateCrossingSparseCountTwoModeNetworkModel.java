package model.count.crossing;

import java.io.PrintStream;
import java.util.Iterator;

import umontreal.iro.lecuyer.probdist.PoissonDist;

import graph.SparseAccess;
import graph.MatrixElement;
import graph.RowElement;
import graph.ColumnElement;

import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class ZeroInflatedNoCovariateCrossingSparseCountTwoModeNetworkModel
		extends NoCovariateCrossingSparseCountTwoModeNetworkModel {

	/* Inflation-admixture parameters */
	protected double[][] gamma;
	protected double[][] prevGamma;
	protected double[][] nu;

	/* The lazy initialization is used here */
	public ZeroInflatedNoCovariateCrossingSparseCountTwoModeNetworkModel() {
	}

	/* This is the last time to finalize model data and parameters */
	@Override
	public void finalizeModel() {
		super.finalizeModel();
		gamma = new double[numOfRowGroups][numOfColumnGroups];
		prevGamma = new double[numOfRowGroups][numOfColumnGroups];
		nu = new double[numOfRowGroups][numOfColumnGroups];
	}

	@Override
	public double getLogProb(int y, int i, int j, int k, int l) {
		if (y > 0) {
			return (Math.log(gamma[k][l]) + getLogTruncatedPoisson(y, i, j, k,
					l));
		} else if (y == 0) {
			return Math.log((1 - gamma[k][l]) + gamma[k][l]
					* PoissonDist.prob(theta[k][l], y));
		} else
			return Double.NaN;
	}

	public double getLogTruncatedPoisson(int y, int i, int j, int k, int l) {
		double prob = PoissonDist.prob(theta[k][l], y);
		if (prob > 0.0)
			return Math.log(prob);
		else
			return -1000;
	}

	@Override
	protected double computeSumOfLogProbs() {

		double sumOfLogProbs = 0;

		try {

			Iterator<MatrixElement> it = ((SparseAccess) graph)
					.elementIterator();
			if (it != null)
				for (; it.hasNext();) {
					MatrixElement element = it.next();
					int i = element.getRow();
					int j = element.getCol();
					int value = (int) element.getValue();
					for (int k = 0; k < numOfRowGroups; k++)
						for (int l = 0; l < numOfColumnGroups; l++) {
							sumOfLogProbs += z[i][k] * w[j][l]
									* getLogProb(value, i, j, k, l);

							sumOfLogProbs -= z[i][k] * w[j][l]
									* getLogProb(0, i, j, k, l);
						}
				}

			for (int k = 0; k < numOfRowGroups; k++)
				for (int l = 0; l < numOfColumnGroups; l++) {
					sumOfLogProbs += (graph.getNumOfRows() * alpha[k])
							* (graph.getNumOfColumns() * beta[l])
							* getLogProb(0, 0, 0, k, l);
				}

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return sumOfLogProbs;
	}

	@Override
	/* This function must use prevW when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		try {

			Iterator<ColumnElement> it = ((SparseAccess) graph).rowIterator(i);
			if (it != null)
				for (; it.hasNext();) {
					ColumnElement element = it.next();
					int j = element.getCol();
					int value = (int) element.getValue();
					for (int l = 0; l < numOfColumnGroups; l++) {
						rowSumOfLogProbs += prevW[j][l]
								* getLogProb(value, i, j, k, l);
						rowSumOfLogProbs -= prevW[j][l]
								* getLogProb(0, i, j, k, l);
					}
				}

			for (int l = 0; l < numOfColumnGroups; l++)
				rowSumOfLogProbs += (graph.getNumOfColumns() * beta[l])
						* getLogProb(0, 0, 0, k, l);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return rowSumOfLogProbs;
	}

	@Override
	/* This function must use prevZ when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		try {

			Iterator<RowElement> it = ((SparseAccess) graph).columnIterator(j);
			if (it != null)
				for (; it.hasNext();) {
					RowElement element = it.next();
					int i = element.getRow();
					int value = (int) element.getValue();
					for (int k = 0; k < numOfRowGroups; k++) {
						columnSumOfLogProbs += prevZ[i][k]
								* getLogProb(value, i, j, k, l);
						columnSumOfLogProbs -= prevZ[i][k]
								* getLogProb(0, i, j, k, l);
					}
				}

			for (int k = 0; k < numOfRowGroups; k++)
				columnSumOfLogProbs += (graph.getNumOfRows() * alpha[k])
						* getLogProb(0, 0, 0, k, l);

		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return columnSumOfLogProbs;
	}

	protected void estimateMeans() throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		ArrayMethods.copyDoubleMatrix(gamma, prevGamma);

		/* Compute nu at t step */
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				nu[k][l] = (1.0 - prevGamma[k][l])
						/ ((1.0 - prevGamma[k][l]) + prevGamma[k][l]
								* PoissonDist.prob(prevTheta[k][l], 0));

		/* Estimate inflation-admixture and mean parameters at t+1 step */
		ArrayMethods.resetDoubleMatrix(gamma);
		ArrayMethods.resetDoubleMatrix(theta);
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			for (; it.hasNext();) {
				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();
				int value = (int) element.getValue();
				for (int k = 0; k < numOfRowGroups; k++)
					for (int l = 0; l < numOfColumnGroups; l++) {
						theta[k][l] += z[i][k] * w[j][l] * value;
						gamma[k][l] += z[i][k] * w[j][l];
						gamma[k][l] -= z[i][k] * w[j][l] * (1 - nu[k][l]);
					}
			}
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				gamma[k][l] += (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]) * (1 - nu[k][l]);

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				theta[k][l] /= gamma[k][l];
				double denominator = (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]);
				gamma[k][l] /= denominator;
			}

		System.out.println("z: ");
		ArrayMethods.printDoubleMatrix(z, " ", 30, System.out);

		System.out.println("w: ");
		ArrayMethods.printDoubleMatrix(w, " ", 30, System.out);

		System.out.println("theta: ");
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);

		System.out.println("gamma: ");
		ArrayMethods.printDoubleMatrix(gamma, " ", System.out);

	}

	/* Reset model parameters */
	@Override
	protected void resetModelParamters() {
		super.resetModelParamters();
		ArrayMethods.resetDoubleMatrix(gamma);
		ArrayMethods.resetDoubleMatrix(prevGamma);
		ArrayMethods.resetDoubleMatrix(nu);
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		super.printModel(outputStream, formatString);
		System.out.println("gamma: ");
		ArrayMethods.printDoubleMatrix(gamma, " ", System.out);
	}

	@Override
	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs((gamma[k][l] - prevGamma[k][l]) / gamma[k][l]) > relativeParameterPrecision)
					return true;

		return false;
	}

}
