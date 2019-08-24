package model.binary.additive;

import graph.ColumnElement;
import graph.RowElement;
import graph.SparseAccess;

import java.util.Iterator;

import optimization.QPSolver;
import utils.ArrayMethods;

public class CyclicNoCovariateAdditiveBinaryTwoModeNetworkModel extends
		NoCovariateAdditiveBinaryTwoModeNetworkModel {

	protected boolean isInterleaving = false;
	protected int rowBlockSize = 100;
	protected int columnBlockSize = 100;

	protected double[] alphaSum;
	protected double[] betaSum;

	/* The lazy initialization is used here */
	public CyclicNoCovariateAdditiveBinaryTwoModeNetworkModel() {
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		alphaSum = new double[numOfRowGroups];
		betaSum = new double[numOfColumnGroups];
	}

	/* This function must use w when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		for (int l = 0; l < numOfColumnGroups; l++) {

			double sumOfW_Y1 = 0;
			Iterator<ColumnElement> it = ((SparseAccess) graph).rowIterator(i);
			if (it != null)
				for (; it.hasNext();) {
					ColumnElement element = it.next();
					int j = element.getCol();
					sumOfW_Y1 += w[j][l];
				}

			rowSumOfLogProbs += sumOfW_Y1 * logPi1[k][l];
			rowSumOfLogProbs += (graph.getNumOfColumns() * beta[l] - sumOfW_Y1)
					* logPi0[k][l];
		}

		return rowSumOfLogProbs;
	}

	/* This function must use z when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		for (int k = 0; k < numOfRowGroups; k++) {

			double sumOfZ_Y1 = 0;
			Iterator<RowElement> it = ((SparseAccess) graph).columnIterator(j);
			if (it != null)
				for (; it.hasNext();) {
					RowElement element = it.next();
					int i = element.getRow();
					sumOfZ_Y1 += z[i][k];
				}

			columnSumOfLogProbs += sumOfZ_Y1 * logPi1[k][l];
			columnSumOfLogProbs += (graph.getNumOfRows() * alpha[k] - sumOfZ_Y1)
					* logPi0[k][l];
		}

		return columnSumOfLogProbs;
	}

	/* Run a full EM iteration by cycling through block updates */
	@Override
	public void runEMIteration(int emIteration, int maxESteps,
			double membershipPrecision, boolean isMMinEStep, boolean isMMinMStep) {

		// Store previous values to test for the improvement of the EM update
		// but all later steps in this procedure will use the current versions:
		// z, w, alpha, beta, theta, gamma
		ArrayMethods.copyDoubleMatrix(z, prevZ);
		ArrayMethods.copyDoubleMatrix(w, prevW);
		ArrayMethods.copyDoubleArray(alpha, prevAlpha);
		ArrayMethods.copyDoubleArray(beta, prevBeta);
		ArrayMethods.copyDoubleArray(theta, prevTheta);
		ArrayMethods.copyDoubleArray(gamma, prevGamma);

		// Compute cached parameters
		for (int k = 0; k < numOfRowGroups; k++)
			alphaSum[k] = graph.getNumOfRows() * alpha[k];
		for (int l = 0; l < numOfColumnGroups; l++)
			betaSum[l] = graph.getNumOfColumns() * beta[l];

		if (isInterleaving)
			interleavingUpdates(maxESteps, membershipPrecision, isMMinEStep,
					isMMinMStep);
		else
			noInterleavingUpdates(maxESteps, membershipPrecision, isMMinEStep,
					isMMinMStep);

	}

	protected void noInterleavingUpdates(int maxESteps,
			double membershipPrecision, boolean isMMinEStep, boolean isMMinMStep) {
		// Update row blocks
		for (int blockIndex = 0; blockIndex < Math
				.ceil((graph.getNumOfRows() / rowBlockSize)); blockIndex++) {
			updateRowBlock(blockIndex, maxESteps, membershipPrecision,
					isMMinEStep);
		}

		// Update column blocks
		for (int blockIndex = 0; blockIndex < Math.ceil((graph
				.getNumOfColumns() / columnBlockSize)); blockIndex++) {
			updateColumnBlock(blockIndex, maxESteps, membershipPrecision,
					isMMinEStep);
		}
	}

	protected void interleavingUpdates(int maxESteps,
			double membershipPrecision, boolean isMMinEStep, boolean isMMinMStep) {

		int rowBlockIndex = 0;
		int numRowBlocks = (int) Math
				.ceil((graph.getNumOfRows() / rowBlockSize));

		int columnBlockIndex = 0;
		int numColumnBlocks = (int) Math
				.ceil((graph.getNumOfColumns() / columnBlockSize));

		while ((rowBlockIndex < numRowBlocks)
				&& (columnBlockIndex < numColumnBlocks)) {
			// Update a row block
			updateRowBlock(rowBlockIndex++, maxESteps, membershipPrecision,
					isMMinEStep);
			// Update a column block
			updateColumnBlock(columnBlockIndex++, maxESteps,
					membershipPrecision, isMMinEStep);
		}

		// Update the tail of row blocks
		for (; rowBlockIndex < numRowBlocks; rowBlockIndex++)
			updateRowBlock(rowBlockIndex, maxESteps, membershipPrecision,
					isMMinEStep);

		// Update the tail of column blocks
		for (; columnBlockIndex < numColumnBlocks; columnBlockIndex++)
			updateColumnBlock(columnBlockIndex, maxESteps, membershipPrecision,
					isMMinEStep);

	}

	protected void updateRowBlock(int rowBlockIndex, int maxESteps,
			double membershipPrecision, boolean isMMinEStep) {

		int fromRow = rowBlockSize * rowBlockIndex;
		int toRow = Math.min(rowBlockSize * (rowBlockIndex + 1),
				graph.getNumOfRows());

		// Update a block of latent row variables
		if (isMMinEStep)
			updateLatentRows_MM(fromRow, toRow);
		else
			updateLatentRows_FP(fromRow, toRow, maxESteps, membershipPrecision);

		// Update alpha
		updateAlpha(fromRow, toRow);

		// Update parameters of the model
		updateModelParamters();

		// Update cached parameters
		updateCachedParameters();

	}

	protected void updateColumnBlock(int columnBlockIndex, int maxESteps,
			double membershipPrecision, boolean isMMinEStep) {

		int fromCol = columnBlockSize * columnBlockIndex;
		int toCol = Math.min(columnBlockSize * (columnBlockIndex + 1),
				graph.getNumOfColumns());

		// Update a block of latent column variables
		if (isMMinEStep)
			updateLatentColumns_MM(fromCol, toCol);
		else
			updateLatentColumns_FP(fromCol, toCol, maxESteps,
					membershipPrecision);

		// Update beta
		updateBeta(fromCol, toCol);

		// Update parameters of the model
		updateModelParamters();

		// Update cached parameters
		updateCachedParameters();

	}

	/*
	 * Update latent variables of rows indexed from fromRow (inclusively) to
	 * toRow (exclusively)
	 */
	protected void updateLatentRows_MM(int fromRow, int toRow) {

		System.out.println("Running MM Rows Cyclic");

		/* Update latent membership variables for rows */
		double[][] zA = new double[toRow - fromRow][numOfRowGroups];
		double[][] zS = new double[toRow - fromRow][numOfRowGroups];

		// Calculate quadratic coefficients
		for (int i = fromRow, index = 0; i < toRow; i++, index++)
			for (int k = 0; k < numOfRowGroups; k++) {
				zA[index][k] = computeRowSumOfLogProbs(i, k);
				// A[i][k] must be equal or less than 0.
				if (zA[index][k] > 0) {
					// However, A[i][k] can be greater than 0 because of
					// numerical precision. We cut it off to 0 in this case
					System.out.println("zA(" + i + ", " + k + ") is "
							+ zA[index][k] + " > 0.");
					zA[index][k] = 0;
				}
				zA[index][k] = 1 - zA[index][k] / 2;
				zA[index][k] /= z[i][k];
			}

		// Calculate linear coefficients
		double[] logAlpha = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			logAlpha[k] = Math.log(alpha[k]);
		for (int i = fromRow, index = 0; i < toRow; i++, index++)
			for (int k = 0; k < numOfRowGroups; k++)
				zS[index][k] = 1 + logAlpha[k] - Math.log(z[i][k]);

		// Solve QP problems
		QPSolver.solve01Bounds(zA, zS, z, fromRow, toRow, minMembership);

		/* Normalize latent membership variables for rows */
		normalizeMembership(z, fromRow, toRow);
	}

	protected void updateLatentRows_FP(int fromRow, int toRow, int maxESteps,
			double membershipPrecision) {

		System.out.println("Running FP Rows Cyclic");

		for (int eStep = 0; eStep < maxESteps; eStep++) {

			/* Update latent membership variables for rows */
			for (int i = fromRow; i < toRow; i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					z[i][k] = Math.log(alpha[k]);
					z[i][k] += computeRowSumOfLogProbs(i, k);
				}

			/*
			 * Convert from log scale to probability scale then normalize latent
			 * membership variables for rows
			 */
			normalizeLogMembership2Membership(z, fromRow, toRow);

			if (!areMembershipVariablesChangedSignificantly(membershipPrecision))
				break;

		}

	}

	protected void updateLatentColumns_MM(int fromCol, int toCol) {

		System.out.println("Running MM Columns Cyclic");

		/* Update latent membership variables for columns */
		double[][] wA = new double[toCol - fromCol][numOfColumnGroups];
		double[][] wS = new double[toCol - fromCol][numOfColumnGroups];

		// Calculate quadratic coefficients
		for (int j = fromCol, index = 0; j < toCol; j++, index++)
			for (int l = 0; l < numOfColumnGroups; l++) {
				wA[index][l] = computeColumnSumOfLogProbs(j, l);
				// A[j][l] must be equal or less than 0.
				if (wA[index][l] > 0) {
					// However, A[j][l] can be greater than 0 because of
					// numerical precision. We cut it off to 0 in this case
					System.out.println("wA(" + j + ", " + l + ") is "
							+ wA[index][l] + " > 0.");
					wA[index][l] = 0;
				}

				wA[index][l] = 1 - wA[index][l] / 2;
				wA[index][l] /= w[j][l];
			}

		// Calculate linear coefficients
		double[] logBeta = new double[numOfColumnGroups];
		for (int l = 0; l < numOfColumnGroups; l++)
			logBeta[l] = Math.log(beta[l]);
		for (int j = fromCol, index = 0; j < toCol; j++, index++)
			for (int l = 0; l < numOfColumnGroups; l++)
				wS[index][l] = 1 + logBeta[l] - Math.log(w[j][l]);

		// Solve QP problems
		QPSolver.solve01Bounds(wA, wS, w, fromCol, toCol, minMembership);

		/* Normalize latent membership variables for columns */
		normalizeMembership(w, fromCol, toCol);
	}

	protected void updateLatentColumns_FP(int fromCol, int toCol,
			int maxESteps, double membershipPrecision) {

		System.out.println("Running FP Columns Cyclic");

		for (int eStep = 0; eStep < maxESteps; eStep++) {

			/* Update latent membership variables for columns */
			for (int j = fromCol; j < toCol; j++)
				for (int l = 0; l < numOfColumnGroups; l++) {
					w[j][l] = Math.log(beta[l]);
					w[j][l] += computeColumnSumOfLogProbs(j, l);
				}

			/*
			 * Convert from log scale to probability scale then normalize latent
			 * membership variables for columns
			 */
			normalizeLogMembership2Membership(w, fromCol, toCol);

			if (!areMembershipVariablesChangedSignificantly(membershipPrecision))
				break;

		}

	}

	protected void updateAlpha(int fromRow, int toRow) {
		/* Update mixing probabilities for rows */
		for (int i = fromRow; i < toRow; i++)
			for (int k = 0; k < numOfRowGroups; k++)
				alphaSum[k] += (z[i][k] - prevZ[i][k]);
		for (int k = 0; k < numOfRowGroups; k++)
			alpha[k] = alphaSum[k] / graph.getNumOfRows();
	}

	protected void updateBeta(int fromCol, int toCol) {
		/* Update mixing probabilities for columns */
		for (int j = fromCol; j < toCol; j++)
			// numOfRowGroups >= 1 AND numOfColumnGroups >= 1
			for (int l = 0; l < numOfColumnGroups; l++)
				betaSum[l] += (w[j][l] - prevW[j][l]);
		for (int l = 0; l < numOfColumnGroups; l++)
			beta[l] = betaSum[l] / graph.getNumOfColumns();
	}

	protected void updateModelParamters() {

		double[] newTheta = new double[numOfRowGroups];
		double[] newGamma = new double[numOfColumnGroups];

		double[][] lambda = new double[numOfRowGroups][numOfColumnGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			for (int l = 0; l < numOfColumnGroups; l++)
				lambda[k][l] = 1 + Math.exp(theta[k] + gamma[l]);

		// Estimate theta: to make the model identifiable, theta[0] = 0 while
		// other parameters are estimated.
		newTheta[0] = 0;
		// If numOfRowGroups = 1, we only need to estimate gamma
		if (numOfRowGroups > 1) {
			for (int k = 1; k < numOfRowGroups; k++) {
				// Compute the numerator
				double numerator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++)
					numerator += z[i][k] * rowDegrees[i];
				// Compute the denominator
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					double term = 0;
					for (int l = 0; l < numOfColumnGroups; l++)
						term += (Math.exp(gamma[l]) / lambda[k][l])
								* (graph.getNumOfColumns() * beta[l]);
					denominator += (z[i][k] / Math.exp(theta[k])) * term;
				}
				newTheta[k] = .5 * Math.log(numerator / denominator);
			}
		}

		// Estimate gamma
		for (int l = 0; l < numOfColumnGroups; l++) {
			// Compute the numerator
			double numerator = 0;
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				numerator += w[j][l] * columnDegrees[j];
			// Compute the denominator
			double denominator = 0;
			for (int j = 0; j < graph.getNumOfColumns(); j++) {
				double term = 0;
				for (int k = 0; k < numOfRowGroups; k++) {
					term += (Math.exp(theta[k]) / lambda[k][l])
							* (graph.getNumOfRows() * alpha[k]);
				}
				denominator += (w[j][l] / Math.exp(gamma[l])) * term;
			}
			newGamma[l] = .5 * Math.log(numerator / denominator);
		}

		ArrayMethods.copyDoubleArray(newTheta, theta);
		ArrayMethods.copyDoubleArray(newGamma, gamma);
	}

	/* Getters and Setters */

	public boolean isInterleaving() {
		return isInterleaving;
	}

	public void setInterleaving(boolean isInterleaving) {
		this.isInterleaving = isInterleaving;
	}

	public int getRowBlockSize() {
		return rowBlockSize;
	}

	public void setRowBlockSize(int rowBlockSize) {
		this.rowBlockSize = rowBlockSize;
	}

	public int getColumnBlockSize() {
		return columnBlockSize;
	}

	public void setColumnBlockSize(int columnBlockSize) {
		this.columnBlockSize = columnBlockSize;
	}

}
