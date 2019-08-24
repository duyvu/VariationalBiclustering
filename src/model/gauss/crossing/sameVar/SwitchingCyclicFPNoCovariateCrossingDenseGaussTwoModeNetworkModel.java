package model.gauss.crossing.sameVar;

import utils.ArrayMethods;

public class SwitchingCyclicFPNoCovariateCrossingDenseGaussTwoModeNetworkModel
		extends CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel {

	protected int thresholdIterations = 1000;

	/* The lazy initialization is used here */
	public SwitchingCyclicFPNoCovariateCrossingDenseGaussTwoModeNetworkModel() {
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
				interleavingUpdates(maxESteps, membershipPrecision, false,
						isMMinMStep);
			else
				noInterleavingUpdates(maxESteps, membershipPrecision, false,
						isMMinMStep);
		}

	}

	/* Getters and Setters */

	public int getThresholdIterations() {
		return thresholdIterations;
	}

	public void setThresholdIterations(int thresholdIterations) {
		this.thresholdIterations = thresholdIterations;
	}

}
