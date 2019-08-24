package model;

import optimization.QPSolver;
import utils.ArrayMethods;

public abstract class ContinuousFrequentistTwoModeNetworkModel extends FrequentistTwoModeNetworkModel {

	public abstract double getLogProb(double y, int i, int j, int k, int l);

	protected void runFP_EStep(int emIteraction, int eStep) {

		/* Store current membership variables to previous versions */
		ArrayMethods.copyDoubleMatrix(z, prevZ);
		ArrayMethods.copyDoubleMatrix(w, prevW);

		/* Update latent membership variables for rows */
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {
				z[i][k] = Math.log(alpha[k]);
				z[i][k] += computeRowSumOfLogProbs(i, k);
			}
		/*
		 * Convert from log scale to probability scale then normalize latent
		 * membership variables for rows
		 */
		normalizeLogMembership2Membership(z);

		/* Update latent membership variables for columns */
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				w[j][l] = Math.log(beta[l]);
				w[j][l] += computeColumnSumOfLogProbs(j, l);
			}
		/*
		 * Convert from log scale to probability scale then normalize latent
		 * membership variables for columns
		 */
		normalizeLogMembership2Membership(w);

	}

	/* This function must use prevW when computing the sum */
	protected abstract double computeRowSumOfLogProbs(int i, int k);

	/* This function must use prevZ when computing the sum */
	protected abstract double computeColumnSumOfLogProbs(int j, int l);

	protected void runMM_EStep(int emIteraction, int eStep) {

		/* Store current membership variables to previous versions */
		ArrayMethods.copyDoubleMatrix(z, prevZ);
		ArrayMethods.copyDoubleMatrix(w, prevW);

		/* Collect coefficients for QP problems */
		/* Quadratic coefficients for rows */
		double[][] zA = new double[graph.getNumOfRows()][numOfRowGroups];
		/* Linear coefficients for rows */
		double[][] zS = new double[graph.getNumOfRows()][numOfRowGroups];
		/* Quadratic coefficients for columns */
		double[][] wA = new double[graph.getNumOfColumns()][numOfColumnGroups];
		/* Linear coefficients for columns */
		double[][] wS = new double[graph.getNumOfColumns()][numOfColumnGroups];

		// Collect coefficients by the row membership distribution
		double[] logAlpha = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			logAlpha[k] = Math.log(alpha[k]);
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {
				zA[i][k] = 1 / prevZ[i][k];
				zS[i][k] = 1 + logAlpha[k] - Math.log(prevZ[i][k]);
			}

		// Collect coefficients by the column membership distribution
		double[] logBeta = new double[numOfColumnGroups];
		for (int k = 0; k < numOfColumnGroups; k++)
			logBeta[k] = Math.log(beta[k]);
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				wA[j][l] = 1 / prevW[j][l];
				wS[j][l] = 1 + logBeta[l] - Math.log(prevW[j][l]);
			}

		try {
			// Collect coefficients by the interaction terms
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int j = 0; j < graph.getNumOfColumns(); j++)
					for (int k = 0; k < numOfRowGroups; k++)
						for (int l = 0; l < numOfColumnGroups; l++) {
							double logProb = getLogProb(graph.getElement(i, j),
									i, j, k, l);
							if (logProb < 0) {
								zA[i][k] += -.5 * (prevW[j][l] / prevZ[i][k])
										* logProb;
								wA[j][l] += -.5 * (prevZ[i][k] / prevW[j][l])
										* logProb;
							} else {
								zA[i][k] += .5 * logProb;
								wA[j][l] += .5 * logProb;
								double linearCoef = (prevZ[i][k] + prevW[j][l])
										* logProb;
								zS[i][k] += linearCoef;
								wS[j][l] += linearCoef;
							}
						}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		// Solve quadratic programming problems for row and column membership
		// variables separately.

		/* Update latent membership variables for rows */
		/* Solve QP problems */
		QPSolver.solve01Bounds(zA, zS, z, minMembership);
		/* Normalize latent membership variables for rows */
		normalizeMembership(z);

		/* Update latent membership variables for columns */
		/* Solve QP problems */
		QPSolver.solve01Bounds(wA, wS, w, minMembership);
		/* Normalize latent membership variables for columns */
		normalizeMembership(w);

	}

}
