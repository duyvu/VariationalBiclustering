package model.binary.crossing;

import java.io.PrintStream;

import model.binary.SparseBinaryFrequentistTwoModeNetworkModel;

import utils.ArrayMethods;

public abstract class CrossingSparseBinaryFrequentistTwoModeNetworkModel extends
		SparseBinaryFrequentistTwoModeNetworkModel {

	protected double[][] theta;
	protected double[][] prevTheta;

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		theta = new double[numOfRowGroups][numOfColumnGroups];
		prevTheta = new double[numOfRowGroups][numOfColumnGroups];
	}

	/* Reset model parameters and every data structure for a new EM run */
	public void reset(boolean isRandom, boolean isMMinMStep) {
		super.reset(isRandom, isMMinMStep);
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		ArrayMethods.resetDoubleMatrix(theta);
		ArrayMethods.resetDoubleMatrix(prevTheta);
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
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);
	}

	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs((theta[k][l] - prevTheta[k][l]) / theta[k][l]) > relativeParameterPrecision)
					return true;

		System.out.print("theta: ");
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);
		System.out.print("prevTheta: ");
		ArrayMethods.printDoubleMatrix(prevTheta, " ", System.out);
		System.out.print("gamma: ");

		System.out.println("relativeParameterPrecision: "
				+ relativeParameterPrecision);

		System.out.println("Paramters theta did not pass the test!");

		return false;
	}
}
