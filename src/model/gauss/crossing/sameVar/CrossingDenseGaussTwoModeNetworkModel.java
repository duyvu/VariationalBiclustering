package model.gauss.crossing.sameVar;

import java.io.PrintStream;

import model.gauss.DenseGaussTwoModeNetworkModel;

import utils.ArrayMethods;

public abstract class CrossingDenseGaussTwoModeNetworkModel extends
		DenseGaussTwoModeNetworkModel {

	/* Both means and variances are varying across groups */
	protected double[][] theta;
	protected double[][] prevTheta;
	protected double sigma2;
	protected double prevSigma2;

	/* This is the last time to finalize model data and parameters */
	@Override
	public void finalizeModel() {
		super.finalizeModel();
		theta = new double[numOfRowGroups][numOfColumnGroups];
		prevTheta = new double[numOfRowGroups][numOfColumnGroups];
	}

	/* Reset model parameters and every data structure for a new EM run */
	@Override
	public void reset(boolean isRandom, boolean isMMinMStep) {
		super.reset(isRandom, isMMinMStep);
	}

	/* Reset model parameters */
	@Override
	protected void resetModelParamters() {
		ArrayMethods.resetDoubleMatrix(theta);
		ArrayMethods.resetDoubleMatrix(prevTheta);
		sigma2 = minScale;
		prevSigma2 = sigma2;
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
		System.out.println("sigma2: " + sigma2);
	}

	@Override
	public boolean areOtherModelParametersChangedSignificantly(
			double relativeParameterPrecision) {

		if (super
				.areOtherModelParametersChangedSignificantly(relativeParameterPrecision))
			return true;

		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				if (Math.abs((theta[k][l] - prevTheta[k][l]) / theta[k][l]) > relativeParameterPrecision)
					return true;

		if (Math.abs((sigma2 - prevSigma2) / sigma2) > relativeParameterPrecision)
			return true;

		// System.out.print("theta: ");
		// ArrayMethods.printDoubleMatrix(theta, " ", System.out);
		// System.out.print("prevTheta: ");
		// ArrayMethods.printDoubleMatrix(prevTheta, " ", System.out);
		// System.out.print("sigma2: " + sigma2);
		// System.out.print("prevSigma2: " + prevSigma2);
		// System.out.println("relativeParameterPrecision: "
		// + relativeParameterPrecision);
		// System.out.println("Theta and Sigma2 did not pass the test!");

		return false;
	}
}
