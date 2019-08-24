package model.binary.crossing.bayes;

import java.io.PrintStream;
import java.util.HashMap;

import model.binary.SparseBinaryBayesianTwoModeNetworkModel;

import utils.ArrayMethods;

public abstract class CrossingSparseBinaryBayesianTwoModeNetworkModel extends
		SparseBinaryBayesianTwoModeNetworkModel {

	protected double _alpha_theta_0;
	protected double[][] alpha_theta_0;

	protected double _beta_theta_0;
	protected double[][] beta_theta_0;

	protected double[][] alpha_theta;
	protected double[][] prev_alpha_theta;

	protected double[][] beta_theta;
	protected double[][] prev_beta_theta;

	public void setHyperParameters(HashMap<String, String> configurationMap) {
		super.setHyperParameters(configurationMap);

		if (configurationMap.get("alpha_theta_0") != null)
			_alpha_theta_0 = Double.parseDouble(configurationMap
					.get("alpha_theta_0"));
		else
			_alpha_theta_0 = 1.0;
		System.out.println("Setting alpha_theta_0 = " + _alpha_theta_0);

		if (configurationMap.get("beta_theta_0") != null)
			_beta_theta_0 = Double.parseDouble(configurationMap
					.get("beta_theta_0"));
		else
			_beta_theta_0 = 1.0;
		System.out.println("Setting beta_theta_0 = " + _beta_theta_0);

	}

	public void setAlpha_Theta_0(double alpha_theta_0) {
		this.alpha_theta_0 = new double[numOfRowGroups][numOfColumnGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				this.alpha_theta_0[k][l] = alpha_theta_0;
	}

	public void setBeta_Theta_0(double beta_theta_0) {
		this.beta_theta_0 = new double[numOfRowGroups][numOfColumnGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				this.beta_theta_0[k][l] = beta_theta_0;
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		alpha_theta = new double[numOfRowGroups][numOfColumnGroups];
		prev_alpha_theta = new double[numOfRowGroups][numOfColumnGroups];
		beta_theta = new double[numOfRowGroups][numOfColumnGroups];
		prev_beta_theta = new double[numOfRowGroups][numOfColumnGroups];

		// Distribute singular hyper-parameters to their matrix forms

		setLambdaAlpha_0(_lambda_alpha_0);
		setLambdaBeta_0(_lambda_beta_0);

		setAlpha_Theta_0(_alpha_theta_0);
		setBeta_Theta_0(_beta_theta_0);

	}

	/* Reset model parameters and every data structure for a new EM run */
	public void reset(boolean isRandom, boolean isMMinMStep) {
		super.reset(isRandom, isMMinMStep);
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		ArrayMethods.resetDoubleMatrix(alpha_theta);
		ArrayMethods.resetDoubleMatrix(prev_alpha_theta);
		ArrayMethods.resetDoubleMatrix(beta_theta);
		ArrayMethods.resetDoubleMatrix(prev_beta_theta);
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		System.out.println("z: ");
		ArrayMethods.printDoubleMatrix(z, " ", System.out);
		System.out.println("lambda_alpha: ");
		ArrayMethods.printDoubleArray(alpha, " ", System.out);
		System.out.println("w: ");
		ArrayMethods.printDoubleMatrix(w, " ", System.out);
		System.out.println("lambda_beta: ");
		ArrayMethods.printDoubleArray(beta, " ", System.out);
		System.out.println("alpha_theta: ");
		ArrayMethods.printDoubleMatrix(alpha_theta, " ", System.out);
		System.out.println("beta_theta: ");
		ArrayMethods.printDoubleMatrix(beta_theta, " ", System.out);
	}

	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs((alpha_theta[k][l] - prev_alpha_theta[k][l])
						/ alpha_theta[k][l]) > relativeParameterPrecision)
					return true;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs((beta_theta[k][l] - prev_beta_theta[k][l])
						/ beta_theta[k][l]) > relativeParameterPrecision)
					return true;

		System.out.print("alpha_theta: ");
		ArrayMethods.printDoubleMatrix(alpha_theta, " ", System.out);
		System.out.print("prev_alpha_theta: ");
		ArrayMethods.printDoubleMatrix(prev_alpha_theta, " ", System.out);

		System.out.print("beta_theta: ");
		ArrayMethods.printDoubleMatrix(beta_theta, " ", System.out);
		System.out.print("prev_beta_theta: ");
		ArrayMethods.printDoubleMatrix(prev_beta_theta, " ", System.out);

		System.out
				.println("Parameters alpha_theta or beta_theta did not pass the test!");

		System.out.println("relativeParameterPrecision: "
				+ relativeParameterPrecision);

		return false;
	}
}
