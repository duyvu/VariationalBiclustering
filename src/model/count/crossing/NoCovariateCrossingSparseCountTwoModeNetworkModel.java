package model.count.crossing;

import java.util.Iterator;

import umontreal.iro.lecuyer.probdist.PoissonDist;

import graph.SparseAccess;
import graph.MatrixElement;
import graph.RowElement;
import graph.ColumnElement;

import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class NoCovariateCrossingSparseCountTwoModeNetworkModel extends
		CrossingSparseCountTwoModeNetworkModel {

	/* The lazy initialization is used here */
	public NoCovariateCrossingSparseCountTwoModeNetworkModel() {
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
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		updateCachedParameters();
	}

	@Override
	public void updateCachedParameters() {
		return;
	}

	@Override
	public double getLogProb(int y, int i, int j, int k, int l) {
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

	@Override
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {
		/* Use the closed-form estimators */
		estimateMeans();
	}

	@Override
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
		/* Use the closed-form estimators */
		estimateMeans();
	}

	protected void estimateMeans() throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);

		/* Estimate mean parameters */
		ArrayMethods.resetDoubleMatrix(theta);
		Iterator<MatrixElement> it = ((SparseAccess) graph).elementIterator();
		if (it != null)
			for (; it.hasNext();) {
				MatrixElement element = it.next();
				int i = element.getRow();
				int j = element.getCol();
				int value = (int) element.getValue();
				for (int k = 0; k < numOfRowGroups; k++)
					for (int l = 0; l < numOfColumnGroups; l++)
						theta[k][l] += z[i][k] * w[j][l] * value;
			}

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				double denominator = (graph.getNumOfRows() * alpha[k])
						* (graph.getNumOfColumns() * beta[l]);
				theta[k][l] /= denominator;
			}

		System.out.println("z: ");
		ArrayMethods.printDoubleMatrix(z, " ", 30, System.out);

		System.out.println("w: ");
		ArrayMethods.printDoubleMatrix(w, " ", 30, System.out);

		System.out.println("theta: ");
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);

	}

	@Override
	protected void updateModelHyperParameters() {
		return;
	}

}
