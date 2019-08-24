package model.gauss.crossing.sameVar;

import optimization.QPSolver;

import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel extends
		NoCovariateCrossingDenseGaussTwoModeNetworkModel {

	protected boolean isInterleaving = false;
	protected int rowBlockSize = 100;
	protected int columnBlockSize = 100;

	protected double[] alphaSum;
	protected double[] betaSum;

	/* The lazy initialization is used here */
	public CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel() {
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
		alphaSum = new double[numOfRowGroups];
		betaSum = new double[numOfColumnGroups];
	}

	@Override
	/* This function must use w when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		try {
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					rowSumOfLogProbs += w[j][l]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		return rowSumOfLogProbs;
	}

	@Override
	/* This function must use z when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					columnSumOfLogProbs += z[i][k]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
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
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		prevSigma2 = sigma2;

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
		try {
			updateModelParamters();
		} catch (UnsupportedMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

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
		try {
			updateModelParamters();
		} catch (UnsupportedMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// Update cached parameters
		updateCachedParameters();

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

	/*
	 * Update latent variables of rows indexed from fromRow (inclusively) to
	 * toRow (exclusively)
	 */
	protected void updateLatentRows_MM(int fromRow, int toRow) {

		System.out.println("Running MM Rows Cyclic");

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

	protected void updateModelParamters() throws UnsupportedMethodException {
		/* Estimate mean parameters */
		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				theta[k][l] = 0;
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					double thetaW = 0;
					double denW = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++) {
						thetaW += w[j][l] * graph.getElement(i, j);
						denW += w[j][l];
					}
					theta[k][l] += z[i][k] * thetaW;
					denominator += z[i][k] * denW;
				}
				theta[k][l] /= denominator;
			}
		}

		/* Estimate variance parameters */
		sigma2 = 0;
		double denominator = 0;
		for (int i = 0; i < graph.getNumOfRows(); i++)
			for (int k = 0; k < numOfRowGroups; k++) {
				double sigma2W = 0;
				double denW = 0;
				for (int j = 0; j < graph.getNumOfColumns(); j++)
					for (int l = 0; l < numOfColumnGroups; l++) {
						sigma2W += w[j][l]
								* Math.pow(
										graph.getElement(i, j) - theta[k][l], 2);
						denW += w[j][l];
					}
				sigma2 += z[i][k] * sigma2W;
				denominator += z[i][k] * denW;
			}
		sigma2 /= denominator;

		System.out.println("theta: ");
		ArrayMethods.printDoubleMatrix(theta, " ", System.out);

		System.out.println("sigma2: " + sigma2);

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
