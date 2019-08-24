package model.count.crossing;

import java.io.PrintStream;

import utils.ArrayMethods;
import model.count.SparseCountTwoModeNetworkModel;

public abstract class CrossingSparseCountTwoModeNetworkModel extends
		SparseCountTwoModeNetworkModel {

	/* Mean parameters */
	protected double[][] theta;
	protected double[][] prevTheta;

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

		return false;
	}

}
