package model.gauss.crossing.sameVar;

import optimization.QPSolver;
import exception.UnsupportedMethodException;
import utils.ArrayMethods;

public class SwitchingBlockRelaxationNoCovariateCrossingDenseGaussTwoModeNetworkModel
		extends CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel {

	protected int thresholdIterations = 1000;

	/* The lazy initialization is used here */
	public SwitchingBlockRelaxationNoCovariateCrossingDenseGaussTwoModeNetworkModel() {
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
				interleavingUpdatesBR(membershipPrecision, isMMinMStep);
			else
				noInterleavingUpdatesBR(membershipPrecision, isMMinMStep);
		}

	}

	protected void interleavingUpdatesBR(double membershipPrecision,
			boolean isMMinMStep) {

		int rowBlockIndex = 0;
		int numRowBlocks = (int) Math
				.ceil((graph.getNumOfRows() / rowBlockSize));

		int columnBlockIndex = 0;
		int numColumnBlocks = (int) Math
				.ceil((graph.getNumOfColumns() / columnBlockSize));

		while ((rowBlockIndex < numRowBlocks)
				&& (columnBlockIndex < numColumnBlocks)) {
			// Update a row block
			updateRowBlock(rowBlockIndex++, membershipPrecision);
			// Update a column block
			updateColumnBlock(columnBlockIndex++, membershipPrecision);
		}

		// Update the tail of row blocks
		for (; rowBlockIndex < numRowBlocks; rowBlockIndex++)
			updateRowBlock(rowBlockIndex, membershipPrecision);

		// Update the tail of column blocks
		for (; columnBlockIndex < numColumnBlocks; columnBlockIndex++)
			updateColumnBlock(columnBlockIndex, membershipPrecision);
	}

	protected void noInterleavingUpdatesBR(double membershipPrecision,
			boolean isMMinMStep) {
		// Update row blocks
		for (int blockIndex = 0; blockIndex < Math
				.ceil((graph.getNumOfRows() / rowBlockSize)); blockIndex++) {
			updateRowBlock(blockIndex, membershipPrecision);
		}

		// Update column blocks
		for (int blockIndex = 0; blockIndex < Math.ceil((graph
				.getNumOfColumns() / columnBlockSize)); blockIndex++) {
			updateColumnBlock(blockIndex, membershipPrecision);
		}
	}

	protected void updateRowBlock(int rowBlockIndex, double membershipPrecision) {
		
		int fromRow = rowBlockSize * rowBlockIndex;
		int toRow = Math.min(rowBlockSize * (rowBlockIndex + 1),
				graph.getNumOfRows());

		// Update a block of latent row variables
		updateLatentRows_BR(fromRow, toRow);

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

	protected void updateColumnBlock(int columnBlockIndex,
			double membershipPrecision) {

		int fromCol = columnBlockSize * columnBlockIndex;
		int toCol = Math.min(columnBlockSize * (columnBlockIndex + 1),
				graph.getNumOfColumns());

		// Update a block of latent column variables
		updateLatentColumns_BR(fromCol, toCol);

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

	/*
	 * Update latent variables of rows indexed from fromRow (inclusively) to
	 * toRow (exclusively)
	 */
	protected void updateLatentRows_BR(int fromRow, int toRow) {

		System.out.println("Running BR Rows Cyclic");

		/* Update latent membership variables for rows */
		double[][] zM = new double[toRow - fromRow][numOfRowGroups];
		double[][] zS = new double[toRow - fromRow][numOfRowGroups];

		double[] logAlpha = new double[numOfRowGroups];
		for (int k = 0; k < numOfRowGroups; k++)
			logAlpha[k] = Math.log(alpha[k]);

		try {
			// Collect coefficients by the interaction terms
			for (int i = fromRow, index = 0; i < toRow; i++, index++)
				for (int k = 0; k < numOfRowGroups; k++) {
					double logProb = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++)
							logProb += w[j][k]
									* getLogProb(graph.getElement(i, j), i, j,
											k, l);
					zM[index][k] = 0.5 / z[i][k];
					zS[index][k] = logAlpha[k] + logProb - Math.log(z[i][k])
							- 1.0;
				}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		// Solve QP problems
		QPSolver.solveBounds(zM, zS, z, fromRow, toRow, minMembership);

		/* Normalize latent membership variables for rows */
		normalizeMembership(z, fromRow, toRow);
	}

	protected void updateLatentColumns_BR(int fromCol, int toCol) {

		System.out.println("Running BR Columns Cyclic");

		/* Update latent membership variables for columns */
		double[][] wM = new double[toCol - fromCol][numOfColumnGroups];
		double[][] wS = new double[toCol - fromCol][numOfColumnGroups];

		double[] logBeta = new double[numOfColumnGroups];
		for (int k = 0; k < numOfColumnGroups; k++)
			logBeta[k] = Math.log(beta[k]);

		try {
			// Collect coefficients by the interaction terms
			for (int j = fromCol, index = 0; j < toCol; j++, index++)
				for (int l = 0; l < numOfColumnGroups; l++) {
					double logProb = 0;
					for (int i = 0; i < graph.getNumOfRows(); i++)
						for (int k = 0; k < numOfRowGroups; k++)
							logProb += z[i][k]
									* getLogProb(graph.getElement(i, j), i, j,
											k, l);
					wM[index][l] = 0.5 / w[j][l];
					wS[index][l] = logBeta[l] + logProb - Math.log(w[j][l]) - 1;
				}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}

		// Solve QP problems
		QPSolver.solveBounds(wM, wS, w, fromCol, toCol, minMembership);

		/* Normalize latent membership variables for columns */
		normalizeMembership(w, fromCol, toCol);
	}

	/* Getters and Setters */

	public int getThresholdIterations() {
		return thresholdIterations;
	}

	public void setThresholdIterations(int thresholdIterations) {
		this.thresholdIterations = thresholdIterations;
	}

}
