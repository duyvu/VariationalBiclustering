package model;

import java.util.HashMap;

import utils.ArrayMethods;

public abstract class FrequentistTwoModeNetworkModel extends
		TwoModeNetworkModel {

	public void setHyperParameters(HashMap<String, String> configurationMap) {
	}

	/* Compute the likelihood lower bound */
	public double getLikelihoodLowerBound() {

		updateCachedParameters();

		double lowerBound = 0;

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				lowerBound += z[i][k] * Math.log(alpha[k] / z[i][k]);

		// System.out.println(lowerBound);

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				lowerBound += w[j][l] * Math.log(beta[l] / w[j][l]);

		// System.out.println(lowerBound);

		double sumOfLogProbs = computeSumOfLogProbs();
		lowerBound += sumOfLogProbs;

		// System.out.println(lowerBound);
		// System.out.println(sumOfLogProbs);

		return lowerBound;
	}

	public void updateLatentHyperParameters() {

		/*
		 * Store current mixing probabilities for rows and columns to previous
		 * versions
		 */
		ArrayMethods.copyDoubleArray(alpha, prevAlpha);
		ArrayMethods.copyDoubleArray(beta, prevBeta);

		/* Update mixing probabilities for rows */
		ArrayMethods.resetDoubleArray(alpha);
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				alpha[k] += z[i][k];
		for (int k = 0; k < numOfRowGroups; k++)
			alpha[k] = alpha[k] / graph.getNumOfRows();

		/* Update mixing probabilities for columns */
		ArrayMethods.resetDoubleArray(beta);
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				beta[l] += w[j][l];
		for (int l = 0; l < numOfColumnGroups; l++)
			beta[l] = beta[l] / graph.getNumOfColumns();

	}

}
