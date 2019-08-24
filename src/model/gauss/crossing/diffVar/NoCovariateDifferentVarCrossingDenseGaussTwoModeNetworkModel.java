package model.gauss.crossing.diffVar;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;

import umontreal.iro.lecuyer.probdist.NormalDist;
import utils.ArrayMethods;
import exception.UnsupportedMethodException;

public class NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel
		extends DifferentVarCrossingDenseGaussTwoModeNetworkModel {

	/* The lazy initialization is used here */
	public NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel() {
	}

	public void setData(String networkDataFile, String rowCovFile,
			String colCovFile, String edgeCovFile, int dataFormat,
			String dataDelim) {
		super.setData(networkDataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
	}

	/* This is the last time to finalize model data and parameters */
	public void finalizeModel() {
		super.finalizeModel();
	}

	/* Reset model parameters */
	protected void resetModelParamters() {
		super.resetModelParamters();
		updateCachedParameters();
	}

	@Override
	public void updateCachedParameters() {
		return;
	}

	@Override
	public double getLogProb(double y, int i, int j, int k, int l) {
		double prob = NormalDist.density(theta[k][l], Math.sqrt(sigma2[k][l]),
				y);
		if (prob > 0.0)
			return Math.log(prob);
		else
			return -1000;
	}

	@Override
	protected double computeSumOfLogProbs() {
		double sumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++) {
					double value = 0;
					for (int j = 0; j < graph.getNumOfColumns(); j++)
						for (int l = 0; l < numOfColumnGroups; l++)
							value += w[j][l]
									* getLogProb(graph.getElement(i, j), i, j,
											k, l);
					sumOfLogProbs += z[i][k] * value;
				}
		} catch (Exception e) {
			System.out.println("Errors in getLogProb(...) "
					+ this.getClass().getName());
			e.printStackTrace();
			System.exit(-1);
		}

		return sumOfLogProbs;
	}

	@Override
	/* This function must use prevW when computing the sum */
	protected double computeRowSumOfLogProbs(int i, int k) {
		double rowSumOfLogProbs = 0;

		try {
			for (int j = 0; j < graph.getNumOfColumns(); j++)
				for (int l = 0; l < numOfColumnGroups; l++)
					rowSumOfLogProbs += prevW[j][l]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			System.out.println("Errors in  computeRowSumOfLogProbs(...) "
					+ this.getClass().getName());
			e.printStackTrace();
			System.exit(-1);
		}

		return rowSumOfLogProbs;
	}

	@Override
	/* This function must use prevZ when computing the sum */
	protected double computeColumnSumOfLogProbs(int j, int l) {
		double columnSumOfLogProbs = 0;

		try {
			for (int i = 0; i < graph.getNumOfRows(); i++)
				for (int k = 0; k < numOfRowGroups; k++)
					columnSumOfLogProbs += prevZ[i][k]
							* getLogProb(graph.getElement(i, j), i, j, k, l);
		} catch (Exception e) {
			System.out.println("Errors in  computeColumnSumOfLogProbs(...) "
					+ this.getClass().getName());
			e.printStackTrace();
			System.exit(-1);
		}

		return columnSumOfLogProbs;
	}

	@Override
	protected void runMM_MStep(int emIteraction)
			throws UnsupportedMethodException {
		/* Use the closed-form estimators */
		estimateMeansVariances();
	}

	@Override
	/* Use Gradient updates to maximize the likelihood lower bound in M step */
	protected void runGradient_MStep(int emIteraction)
			throws UnsupportedMethodException {
	}

	/*
	 * TO-DO: Exchange loops to factorize z[i][k] out to improve the
	 * performance, i.e., reduce the number of multiply operations
	 */
	protected void estimateMeansVariances() throws UnsupportedMethodException {

		/* Store current model parameters to previous versions */
		ArrayMethods.copyDoubleMatrix(theta, prevTheta);
		ArrayMethods.copyDoubleMatrix(sigma2, prevSigma2);

		/* Estimate mean parameters */
		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				theta[k][l] = 0;
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					for (int j = 0; j < graph.getNumOfColumns(); j++) {
						theta[k][l] += z[i][k] * w[j][l]
								* graph.getElement(i, j);
						denominator += z[i][k] * w[j][l];
					}
				}
				theta[k][l] /= denominator;
			}
		}

		/* Estimate variance parameters */
		for (int k = 0; k < numOfRowGroups; k++) {
			for (int l = 0; l < numOfColumnGroups; l++) {
				sigma2[k][l] = 0;
				double denominator = 0;
				for (int i = 0; i < graph.getNumOfRows(); i++) {
					for (int j = 0; j < graph.getNumOfColumns(); j++) {
						sigma2[k][l] += z[i][k]
								* w[j][l]
								* Math.pow(
										graph.getElement(i, j) - theta[k][l], 2);
						denominator += z[i][k] * w[j][l];
					}
				}
				sigma2[k][l] /= denominator;
			}
		}

	}

	@Override
	protected void updateModelHyperParameters() {
		return;
	}

	public double[] computeMisRates(String outputZFilename,
			String outputWFilename) {

		double[] misRates = new double[3];

		try {

			BufferedReader rowReader = new BufferedReader(new FileReader(
					new File(outputZFilename)));
			int numRows = Integer.parseInt(rowReader.readLine());
			int[] trueZs = new int[numRows];
			for (int i = 0; i < numRows; i++)
				trueZs[i] = Integer.parseInt(rowReader.readLine());
			rowReader.close();
			// ArrayMethods.printIntegerArray(trueZs, "\t", System.out);

			double[] rowTheta = new double[getNumOfRowGroups()];
			ArrayMethods.resetDoubleArray(rowTheta);
			for (int k = 0; k < getNumOfRowGroups(); k++)
				for (int l = 0; l < getNumOfColumnGroups(); l++)
					rowTheta[k] += theta[k][l];

			TreeMap<Double, Integer> sortedRowThetas = new TreeMap<Double, Integer>();
			for (int k = 0; k < getNumOfRowGroups(); k++)
				sortedRowThetas.put(rowTheta[k], k);
			int[] rowLabels = new int[getNumOfRowGroups()];
			int rowRank = 0;
			for (int k : sortedRowThetas.values())
				rowLabels[k] = rowRank++;
			ArrayMethods.printIntegerArray(rowLabels, "\t", System.out);

			BufferedReader colReader = new BufferedReader(new FileReader(
					new File(outputWFilename)));
			int numCols = Integer.parseInt(colReader.readLine());
			int[] trueWs = new int[numCols];
			for (int j = 0; j < numCols; j++)
				trueWs[j] = Integer.parseInt(colReader.readLine());
			colReader.close();
			// ArrayMethods.printIntegerArray(trueWs, "\t", System.out);

			double[] colTheta = new double[getNumOfColumnGroups()];
			ArrayMethods.resetDoubleArray(colTheta);
			for (int l = 0; l < getNumOfColumnGroups(); l++)
				for (int k = 0; k < getNumOfRowGroups(); k++)
					colTheta[l] += theta[k][l];

			TreeMap<Double, Integer> sortedColThetas = new TreeMap<Double, Integer>();
			for (int l = 0; l < getNumOfColumnGroups(); l++)
				sortedColThetas.put(colTheta[l], l);
			int[] colLabels = new int[getNumOfColumnGroups()];
			int colRank = 0;
			for (int l : sortedColThetas.values())
				colLabels[l] = colRank++;
			ArrayMethods.printIntegerArray(colLabels, "\t", System.out);

			double rowRate = 0.0;
			for (int i = 0; i < graph.getNumOfRows(); i++) {
				int myLabel = rowLabels[whichMax(z[i])];
				if (myLabel != trueZs[i])
					rowRate += 1.0;
			}
			rowRate /= graph.getNumOfRows();

			double colRate = 0.0;
			for (int j = 0; j < graph.getNumOfColumns(); j++) {
				int myLabel = colLabels[whichMax(w[j])];
				if (myLabel != trueWs[j])
					colRate += 1.0;
			}
			colRate /= graph.getNumOfColumns();

			double totalRate = (rowRate * graph.getNumOfRows() + colRate
					* graph.getNumOfColumns())
					/ (graph.getNumOfRows() + graph.getNumOfColumns());

			misRates[0] = rowRate;
			misRates[1] = colRate;
			misRates[2] = totalRate;

		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return misRates;
	}

	public static int whichMax(double[] values) {
		if (values.length < 1)
			return -1;

		int maxIndex = 0;
		double maxValue = values[0];
		for (int k = 1; k < values.length; k++) {
			if (maxValue < values[k]) {
				maxValue = values[k];
				maxIndex = k;
			}
		}
		return maxIndex;

	}

	/*
	 * -2 * LB + ((K-1) + (L-1) + 2 * K * L + numRowCovs*numCols +
	 * numColCovs*numRows)*log(numRows*numCols)
	 */
	public double computeBIC() {

		double BIC = -2 * getLikelihoodLowerBound();

		double mixingParameters = numOfRowGroups - 1.0 + numOfColumnGroups
				- 1.0;
		double modelParamters = 2 * numOfRowGroups * numOfColumnGroups;
		double covariateParameters = 0.0;
		BIC += (mixingParameters + modelParamters + covariateParameters)
				* Math.log(graph.getNumOfRows() * graph.getNumOfColumns());

		return BIC;
	}

}
