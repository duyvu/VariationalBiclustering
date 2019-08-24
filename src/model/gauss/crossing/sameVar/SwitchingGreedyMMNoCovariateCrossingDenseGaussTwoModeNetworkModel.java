package model.gauss.crossing.sameVar;

import optimization.QPSolver;
import utils.ArrayMethods;

public class SwitchingGreedyMMNoCovariateCrossingDenseGaussTwoModeNetworkModel
		extends CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel {

	protected int thresholdIterations = 1000;
	protected double greedyStep = 2.0;

	/* The lazy initialization is used here */
	public SwitchingGreedyMMNoCovariateCrossingDenseGaussTwoModeNetworkModel() {
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
	}

	/* Run a full EM iteration by cycling through block updates */
	@Override
	public void runEMIteration(int emIteration, int maxESteps,
			double membershipPrecision, boolean isMMinEStep, boolean isMMinMStep) {

		if (emIteration < thresholdIterations) {

			// E Step
			for (int eStep = 0; eStep < maxESteps; eStep++) {
				runEStep(emIteration, eStep, isMMinEStep);
				if (isMMinEStep
						|| !areMembershipVariablesChangedSignificantly(membershipPrecision))
					break;
			}

			// M Step
			runMStep(emIteration, isMMinMStep);

		} else {

			// Store previous values to test for the improvement of the EM
			// update but all later steps in this procedure will use the current
			// versions: z, w, alpha, beta, theta, gamma
			ArrayMethods.copyDoubleMatrix(z, prevZ);
			ArrayMethods.copyDoubleMatrix(w, prevW);
			ArrayMethods.copyDoubleArray(alpha, prevAlpha);
			ArrayMethods.copyDoubleArray(beta, prevBeta);
			ArrayMethods.copyDoubleMatrix(theta, prevTheta);
			prevSigma2 = sigma2;

			// Compute cached parameters
			for (int k = 0; k < numOfRowGroups; k++)
				alphaSum[k] = graph.getNumOfRows() * alpha[k];
			for (int l = 0; l < numOfColumnGroups; l++)
				betaSum[l] = graph.getNumOfColumns() * beta[l];

			if (isInterleaving)
				interleavingUpdates(maxESteps, membershipPrecision,
						isMMinEStep, isMMinMStep);
			else
				noInterleavingUpdates(maxESteps, membershipPrecision,
						isMMinEStep, isMMinMStep);
		}

		int convergedMembershipCounts = 0;
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				if (z[i][k] > .90) {
					convergedMembershipCounts++;
					break;
				}
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (w[j][l] > .90) {
					convergedMembershipCounts++;
					break;
				}

		System.out.println("The percentage of converged memberships is " + 1.0
				* convergedMembershipCounts
				/ (graph.getNumOfRows() + graph.getNumOfColumns()));

	}

	/*
	 * Update latent variables of rows indexed from fromRow (inclusively) to
	 * toRow (exclusively)
	 */
	protected void updateLatentRows_MM(int fromRow, int toRow) {

		System.out.println("Running Greedy MM Rows Cyclic");

		/* Update latent membership variables for rows */
		double[][] zA = new double[toRow - fromRow][numOfRowGroups];
		double[][] zS = new double[toRow - fromRow][numOfRowGroups];

		double[] logAlpha = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			logAlpha[k] = Math.log(alpha[k]);

		try {
			// Collect coefficients by the interaction terms
			for (int i = fromRow, index = 0; i < toRow; i++, index++)
				for (int k = 0; k < numOfRowGroups; k++)
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++) {
							double logProb = getLogProb(graph.getElement(i, j),
									i, j, k, l);
							if (logProb < 0) {
								zA[index][k] += -.5 * (w[j][l] / z[i][k])
										* logProb;

							} else {
								zA[index][k] += .5 * logProb;
								double linearCoef = (z[i][k] + w[j][l])
										* logProb;
								zS[index][k] += linearCoef;
							}
						}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		for (int i = fromRow, index = 0; i < toRow; i++, index++)
			for (int k = 0; k < numOfRowGroups; k++) {
				zA[index][k] += 1 / z[i][k];
				zS[index][k] = 1 + logAlpha[k] - Math.log(z[i][k]);
			}

		// Solve QP problems
		QPSolver.solve01Bounds(zA, zS, z, fromRow, toRow, minMembership);

		double oneStepLB = getLikelihoodLowerBound();

		// Apply a more GREEDY step length but save the one step update first
		double[][] oneStepZ = new double[toRow - fromRow][numOfRowGroups];
		for (int rowIndex = fromRow, i = 0; rowIndex < toRow; rowIndex++, i++)
			for (int k = 0; k < numOfRowGroups; k++)
				oneStepZ[i][k] = z[rowIndex][k];

		while (greedyStep > 1.0) {

			for (int rowIndex = fromRow; rowIndex < toRow; rowIndex++)
				for (int k = 0; k < numOfRowGroups; k++)
					z[rowIndex][k] = prevZ[rowIndex][k] + greedyStep
							* (z[rowIndex][k] - prevZ[rowIndex][k]);

			double greedyStepLB = getLikelihoodLowerBound();

			// Decrease the greedy step if the current greedy step fails
			if (greedyStepLB <= oneStepLB) {
				greedyStep /= 2.0;
				for (int rowIndex = fromRow, i = 0; rowIndex < toRow; rowIndex++, i++)
					for (int k = 0; k < numOfRowGroups; k++)
						z[rowIndex][k] = oneStepZ[i][k];
			} else
				break;
		}

		/* Normalize latent membership variables for rows */
		normalizeMembership(z, fromRow, toRow);
	}

	protected void updateLatentColumns_MM(int fromCol, int toCol) {

		System.out.println("Running Greedy MM Columns Cyclic");

		/* Update latent membership variables for columns */
		double[][] wA = new double[toCol - fromCol][numOfColumnGroups];
		double[][] wS = new double[toCol - fromCol][numOfColumnGroups];

		double[] logBeta = new double[numOfColumnGroups];
		for (int k = 0; k < numOfColumnGroups; k++)
			logBeta[k] = Math.log(beta[k]);

		try {
			// Collect coefficients by the interaction terms
			for (int j = fromCol, index = 0; j < toCol; j++, index++)
				for (int l = 0; l < numOfColumnGroups; l++)
					for (int i = 0; i < graph.getNumOfRows(); i++)
						for (int k = 0; k < numOfRowGroups; k++) {
							double logProb = getLogProb(graph.getElement(i, j),
									i, j, k, l);
							if (logProb < 0) {
								wA[index][l] += -.5 * (z[i][k] / w[j][l])
										* logProb;
							} else {
								wA[index][l] += .5 * logProb;
								double linearCoef = (z[i][k] + w[j][l])
										* logProb;
								wS[index][l] += linearCoef;
							}
						}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		for (int j = fromCol, index = 0; j < toCol; j++, index++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				wA[index][l] += 1 / w[j][l];
				wS[index][l] = 1 + logBeta[l] - Math.log(w[j][l]);
			}

		// Solve QP problems
		QPSolver.solve01Bounds(wA, wS, w, fromCol, toCol, minMembership);

		double oneStepLB = getLikelihoodLowerBound();

		// Apply a more GREEDY step length but save the one step update first
		double[][] oneStepW = new double[toCol - fromCol][numOfColumnGroups];
		for (int rowIndex = fromCol, i = 0; rowIndex < toCol; rowIndex++, i++)
			for (int k = 0; k < numOfColumnGroups; k++)
				oneStepW[i][k] = w[rowIndex][k];

		while (greedyStep > 1.0) {

			for (int rowIndex = fromCol; rowIndex < toCol; rowIndex++)
				for (int k = 0; k < numOfColumnGroups; k++)
					w[rowIndex][k] = prevW[rowIndex][k] + greedyStep
							* (w[rowIndex][k] - prevW[rowIndex][k]);

			double greedyStepLB = getLikelihoodLowerBound();

			// Decrease the greedy step if the current greedy step fails
			if (greedyStepLB <= oneStepLB) {
				greedyStep /= 2.0;
				for (int rowIndex = fromCol, i = 0; rowIndex < toCol; rowIndex++, i++)
					for (int k = 0; k < numOfColumnGroups; k++)
						w[rowIndex][k] = oneStepW[i][k];
			} else
				break;
		}

		/* Normalize latent membership variables for columns */
		normalizeMembership(w, fromCol, toCol);
	}

	/* Getters and Setters */

	public int getThresholdIterations() {
		return thresholdIterations;
	}

	public void setThresholdIterations(int thresholdIterations) {
		this.thresholdIterations = thresholdIterations;
	}

	public double getGreedyStep() {
		return greedyStep;
	}

	public void setGreedyStep(double greedyStep) {
		this.greedyStep = greedyStep;
	}
}
