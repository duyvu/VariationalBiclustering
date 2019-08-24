package model;

import java.util.HashMap;

import org.apache.commons.math3.special.Gamma;

import utils.ArrayMethods;

public abstract class BayesianTwoModeNetworkModel extends TwoModeNetworkModel {

	protected double _lambda_alpha_0;
	protected double[] lambda_alpha_0 = null;

	protected double _lambda_beta_0;
	protected double[] lambda_beta_0 = null;

	public void setHyperParameters(HashMap<String, String> configurationMap) {
		if (configurationMap.get("lambda_alpha_0") != null)
			_lambda_alpha_0 = Double.parseDouble(configurationMap
					.get("lambda_alpha_0"));
		else
			_lambda_alpha_0 = 1.0;
		System.out.println("Setting lambda_alpha_0 = " + _lambda_alpha_0);

		if (configurationMap.get("lambda_beta_0") != null)
			_lambda_beta_0 = Double.parseDouble(configurationMap
					.get("lambda_beta_0"));
		else
			_lambda_beta_0 = 1.0;
		System.out.println("Setting lambda_beta_0 = " + _lambda_beta_0);

	}

	/* Compute the likelihood lower bound */
	public double getLikelihoodLowerBound() {

		updateCachedParameters();

		double lowerBound = 0;

		// Rows

		// Hyper

		lowerBound -= Gamma.logGamma(ArrayMethods.sumDoubleArray(alpha))
				- Gamma.logGamma(ArrayMethods.sumDoubleArray(lambda_alpha_0));

		for (int k = 0; k < numOfRowGroups; k++)
			lowerBound -= Gamma.logGamma(lambda_alpha_0[k])
					- Gamma.logGamma(alpha[k]);

		double digammaSumAlpha = Gamma.digamma(ArrayMethods
				.sumDoubleArray(alpha));

		double[] digammaAlpha = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			digammaAlpha[k] = Gamma.digamma(alpha[k]);

		for (int k = 0; k < numOfRowGroups; k++)
			lowerBound -= (alpha[k] - lambda_alpha_0[k])
					* (digammaAlpha[k] - digammaSumAlpha);

		lowerBound -= digammaSumAlpha * graph.getNumOfRows();

		// Membership

		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				lowerBound -= z[i][k] * Math.log(z[i][k]);

		for (int k = 0; k < numOfRowGroups; k++) {
			double sumZ_k = 0.0;
			for (int i = 0; i < graph.getNumOfRows(); i++)
				sumZ_k += z[i][k];
			lowerBound += digammaAlpha[k] * sumZ_k;
		}

		// Columns

		// Hyper

		lowerBound -= Gamma.logGamma(ArrayMethods.sumDoubleArray(beta))
				- Gamma.logGamma(ArrayMethods.sumDoubleArray(lambda_beta_0));

		for (int l = 0; l < numOfColumnGroups; l++)
			lowerBound -= Gamma.logGamma(lambda_beta_0[l])
					- Gamma.logGamma(beta[l]);

		double[] digammaBeta = new double[numOfColumnGroups];
		for (int l = 0; l < numOfColumnGroups; l++)
			digammaBeta[l] = Gamma.digamma(beta[l]);

		double digammaSumBeta = Gamma
				.digamma(ArrayMethods.sumDoubleArray(beta));

		for (int l = 0; l < numOfColumnGroups; l++)
			lowerBound -= (beta[l] - lambda_beta_0[l])
					* (digammaBeta[l] - digammaSumBeta);

		lowerBound -= digammaSumBeta * graph.getNumOfColumns();

		// Membership

		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				lowerBound -= w[j][l] * Math.log(w[j][l]);

		for (int l = 0; l < numOfColumnGroups; l++) {
			double sumW_l = 0.0;
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				sumW_l += w[j][l];
			lowerBound += Gamma.digamma(beta[l]) * sumW_l;
		}

		// Model

		double sumModelPriorKL = sumModelPriorKL();
		lowerBound -= sumModelPriorKL;

		double sumOfLogProbs = computeSumOfLogProbs();
		lowerBound += sumOfLogProbs;

		return lowerBound;
	}

	/* Sum over K x L KL divergence of model parameters */
	protected abstract double sumModelPriorKL();

	public void updateLatentHyperParameters() {

		/*
		 * Store current mixing probabilities for rows and columns to previous
		 * versions
		 */
		ArrayMethods.copyDoubleArray(alpha, prevAlpha);
		ArrayMethods.copyDoubleArray(beta, prevBeta);

		/* Update mixing probabilities for rows */
		ArrayMethods.copyDoubleArray(lambda_alpha_0, alpha);
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++)
				alpha[k] += z[i][k];

		/* Update mixing probabilities for columns */
		ArrayMethods.copyDoubleArray(lambda_beta_0, beta);
		for (int j = 0; j < graph.getNumOfColumns(); j++)
			for (int l = 0; l < numOfColumnGroups; l++)
				beta[l] += w[j][l];

	}

	/* Setters and Getters */
	public double[] getLambdaAlpha() {
		return lambda_alpha_0;
	}

	public void setLambdaAlpha_0(double lambda_alpha) {
		this.lambda_alpha_0 = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			this.lambda_alpha_0[k] = lambda_alpha;
	}

	public void setLambdaAlpha(double[] lambda_alpha) {
		this.lambda_alpha_0 = new double[numOfRowGroups];
		ArrayMethods.copyDoubleArray(lambda_alpha, this.lambda_alpha_0);
	}

	public double[] getLambdaBeta() {
		return lambda_beta_0;
	}

	public void setLambdaBeta_0(double lambda_beta) {
		this.lambda_beta_0 = new double[numOfColumnGroups];
		for (int l = 0; l < numOfColumnGroups; l++)
			this.lambda_beta_0[l] = lambda_beta;
	}

	public void setLambdaBeta(double[] lambda_beta) {
		this.lambda_beta_0 = new double[numOfColumnGroups];
		ArrayMethods.copyDoubleArray(lambda_beta, this.lambda_beta_0);
	}

}
