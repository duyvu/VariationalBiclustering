package model.binary.additive;

import java.io.PrintStream;

import model.binary.SparseBinaryFrequentistTwoModeNetworkModel;

import utils.ArrayMethods;

public abstract class AdditiveSparseBinaryTwoModeNetworkModel extends
		SparseBinaryFrequentistTwoModeNetworkModel {

	protected double[] theta;
	protected double[] prevTheta;
	protected double[] gamma;
	protected double[] prevGamma;

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		theta = new double[numOfRowGroups];
		prevTheta = new double[numOfRowGroups];
		gamma = new double[numOfColumnGroups];
		prevGamma = new double[numOfColumnGroups];
	}

	/* Reset model parameters and every data structure for a new EM run */
	public void reset(boolean isRandom, boolean isMMinMStep) {
		super.reset(isRandom, isMMinMStep);
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		ArrayMethods.resetDoubleArray(theta);
		ArrayMethods.resetDoubleArray(prevTheta);
		ArrayMethods.resetDoubleArray(gamma);
		ArrayMethods.resetDoubleArray(prevGamma);
	}

	@Override
	public void printModel(PrintStream outputStream, String formatString) {
		System.out.println("z: ");
		ArrayMethods.printDoubleMatrix(z, " ", System.out);
		System.out.println("alpha: ");
		ArrayMethods.printDoubleArray(alpha, " ", System.out);
		System.out.println("w: ");
		ArrayMethods.printDoubleMatrix(w, " ", System.out);
		System.out.println("beta: ");
		ArrayMethods.printDoubleArray(beta, " ", System.out);
		System.out.println("theta: ");
		ArrayMethods.printDoubleArray(theta, " ", System.out);
		System.out.println("gamma: ");
		ArrayMethods.printDoubleArray(gamma, " ", System.out);
	}

	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int k = 0; k < numOfRowGroups; k++)
			if (Math.abs((theta[k] - prevTheta[k]) / theta[k]) > relativeParameterPrecision)
				return true;

		for (int l = 0; l < numOfColumnGroups; l++)
			if (Math.abs((gamma[l] - prevGamma[l]) / gamma[l]) > relativeParameterPrecision)
				return true;

		System.out.print("theta: ");
		ArrayMethods.printDoubleArray(theta, " ", System.out);
		System.out.print("prevTheta: ");
		ArrayMethods.printDoubleArray(prevTheta, " ", System.out);
		System.out.print("gamma: ");
		ArrayMethods.printDoubleArray(gamma, " ", System.out);
		System.out.print("prevGamma: ");
		ArrayMethods.printDoubleArray(prevGamma, " ", System.out);

		System.out.println("relativeParameterPrecision: "
				+ relativeParameterPrecision);

		System.out.println("Theta va Gamma did not pass the test!");

		return false;
	}
}
