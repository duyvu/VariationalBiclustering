package model;

import optimization.QPSolver;

import org.apache.commons.math.special.Gamma;

import utils.ArrayMethods;

public abstract class DiscreteBayesianTwoModeNetworkModel extends
		BayesianTwoModeNetworkModel {

	public abstract double getLogProb(int y, int i, int j, int k, int l);

	protected void runFP_EStep(int emIteraction, int eStep) {

		/* Store current membership variables to previous versions */
		ArrayMethods.copyDoubleMatrix(z, prevZ);
		ArrayMethods.copyDoubleMatrix(w, prevW);

		/* Update latent membership variables for rows */
		double digammaSumAlpha = Gamma.digamma(ArrayMethods
				.sumDoubleArray(alpha));
		double[] digammmaAlphaDifference = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			digammmaAlphaDifference[k] = Gamma.digamma(alpha[k])
					- digammaSumAlpha;

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {
				z[i][k] = digammmaAlphaDifference[k];
				z[i][k] += computeRowSumOfLogProbs(i, k);
			}

		/*
		 * Convert from log scale to probability scale then normalize latent
		 * membership variables for rows
		 */
		normalizeLogMembership2Membership(z);

		double digammaSumBeta = Gamma
				.digamma(ArrayMethods.sumDoubleArray(beta));
		double[] digammmaBetaDifference = new double[numOfColumnGroups];
		for (int l = 0; l < numOfColumnGroups; l++)
			digammmaBetaDifference[l] = Gamma.digamma(beta[l]) - digammaSumBeta;

		/* Update latent membership variables for columns */
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				w[j][l] = digammmaBetaDifference[l];
				w[j][l] += computeColumnSumOfLogProbs(j, l);
			}

		/*
		 * Convert from log scale to probability scale then normalize latent
		 * membership variables for columns
		 */
		normalizeLogMembership2Membership(w);

	}

	protected void runMM_EStep(int emIteraction, int eStep) {

		/* Store current membership variables to previous versions */
		ArrayMethods.copyDoubleMatrix(z, prevZ);
		ArrayMethods.copyDoubleMatrix(w, prevW);

		/* Update latent membership variables for rows */
		runMMEStep4Rows();

		/* Update latent membership variables for columns */
		runMMEStep4Columns();

	}

	protected void runMMEStep4Rows() {
		/* Update latent membership variables for rows */
		double[][] zA = new double[graph.getNumOfRows()][numOfRowGroups];
		double[][] zS = new double[graph.getNumOfRows()][numOfRowGroups];

		// Calculate quadratic coefficients
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {

				zA[i][k] = computeRowSumOfLogProbs(i, k);

				// A[i][k] must be equal or less than 0.
				if (zA[i][k] > 0) {
					// However, A[i][k] can be greater than 0 because of
					// numerical precision. We cut it off to 0 in this case
					System.out.println("zA(" + i + ", " + k + ") is "
							+ zA[i][k] + " > 0.");
					zA[i][k] = 0;
				}

				zA[i][k] = 1 - zA[i][k] / 2;
				zA[i][k] /= prevZ[i][k];
			}

		// Calculate linear coefficients
		double digammaSumAlpha = Gamma.digamma(ArrayMethods
				.sumDoubleArray(alpha));
		double[] digammmaAlphaDifference = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			digammmaAlphaDifference[k] = Gamma.digamma(alpha[k])
					- digammaSumAlpha;

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				zS[i][k] = 1 + digammmaAlphaDifference[k]
						- Math.log(prevZ[i][k]);

		// Solve QP problems
		QPSolver.solve01Bounds(zA, zS, z, minMembership);

		/* Normalize latent membership variables for rows */
		normalizeMembership(z);
	}

	protected void runMMEStep4Columns() {
		/* Update latent membership variables for columns */
		double[][] wA = new double[graph.getNumOfColumns()][numOfColumnGroups];
		double[][] wS = new double[graph.getNumOfColumns()][numOfColumnGroups];

		// Calculate quadratic coefficients
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++) {

				wA[j][l] = computeColumnSumOfLogProbs(j, l);

				// A[j][l] must be equal or less than 0.
				if (wA[j][l] > 0) {
					// However, A[j][l] can be greater than 0 because of
					// numerical precision. We cut it off to 0 in this case
					System.out.println("wA(" + j + ", " + l + ") is "
							+ wA[j][l] + " > 0.");
					wA[j][l] = 0;
				}

				wA[j][l] = 1 - wA[j][l] / 2;
				wA[j][l] /= prevW[j][l];
			}

		// Calculate linear coefficients
		double digammaSumBeta = Gamma
				.digamma(ArrayMethods.sumDoubleArray(beta));
		double[] digammmaBetaDifference = new double[numOfColumnGroups];
		for (int l = 0; l < numOfColumnGroups; l++)
			digammmaBetaDifference[l] = Gamma.digamma(beta[l]) - digammaSumBeta;
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				wS[j][l] = 1 + digammmaBetaDifference[l]
						- Math.log(prevW[j][l]);

		// Solve QP problems
		QPSolver.solve01Bounds(wA, wS, w, minMembership);

		/* Normalize latent membership variables for columns */
		normalizeMembership(w);
	}

	/* This function must use prevW when computing the sum */
	protected abstract double computeRowSumOfLogProbs(int i, int k);

	/* This function must use prevZ when computing the sum */
	protected abstract double computeColumnSumOfLogProbs(int j, int l);

}
