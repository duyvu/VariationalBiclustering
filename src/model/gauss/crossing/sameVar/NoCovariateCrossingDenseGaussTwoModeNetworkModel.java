package model.gauss.crossing.sameVar;

import umontreal.iro.lecuyer.probdist.NormalDist;
import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class NoCovariateCrossingDenseGaussTwoModeNetworkModel extends
		CrossingDenseGaussTwoModeNetworkModel {

	/* The lazy initialization is used here */
	public NoCovariateCrossingDenseGaussTwoModeNetworkModel() {
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
	public double getLogProb(double y, int i, int j, int k, int l) {

		// System.out.println("y =" + y);
		// System.out.println("k =" + k);
		// System.out.println("l =" + l);
		// System.out.println("theta[k][l] =" + theta[k][l]);
		// System.out.println("Math.sqrt(sigma2) =" + Math.sqrt(sigma2));
		// System.out
		// .println("Math.log(P(...)) ="
		// + Math.log(NormalDist.density(theta[k][l],
		// Math.sqrt(sigma2), y)));

		return Math.log(NormalDist.density(theta[k][l], Math.sqrt(sigma2), y));
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
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					rowSumOfLogProbs += prevW[j][l]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
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
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					columnSumOfLogProbs += prevZ[i][k]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
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
		estimateMeansVariances();
	}

	@Override
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
		/* Use the closed-form estimators */
		estimateMeansVariances();
	}

	protected void estimateMeansVariances() throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		prevSigma2 = sigma2;

		/* Estimate mean parameters */
		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				theta[k][l] = 0;
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					double thetaW = 0;
					double denW = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++) {
						thetaW += w[j][l] * graph.getElement(i, j);
						denW += w[j][l];
					}
					theta[k][l] += z[i][k] * thetaW;
					denominator += z[i][k] * denW;
				}
				theta[k][l] /= denominator;
			}
		}

		/* Estimate variance parameters */
		sigma2 = 0;
		double denominator = 0;
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {
				double sigma2W = 0;
				double denW = 0;
				for (int j = 0; j < graph.getNumOfColumns(); j++)
					for (int l = 0; l < numOfColumnGroups; l++) {
						sigma2W += w[j][l]
								* Math.pow(
										graph.getElement(i, j) - theta[k][l], 2);
						denW += w[j][l];
					}
				sigma2 += z[i][k] * sigma2W;
				denominator += z[i][k] * denW;
			}
		sigma2 /= denominator;

		System.out.println("z: ");
		ArrayMethods.printDoubleMatrix(z, " ", 30, System.out);

		System.out.println("w: ");
		ArrayMethods.printDoubleMatrix(w, " ", 30, System.out);

		System.out.println("theta: ");
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);

		System.out.println("sigma2: " + sigma2);

	}

	@Override
	protected void updateModelHyperParameters() {
		return;
	}

}
