package model;

import org.apache.commons.math.optimization.general.ConjugateGradientFormula;

import umontreal.iro.lecuyer.rng.RandomStream;

public abstract class NetworkModel {

	public final static int SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK = 0;

	public final static int DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 1;

	public final static int COVARIATE_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK = 2;

	public final static int COVARIATE_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 3;

	public final static int ROW_COL_VARYING_COVARIATES_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK = 4;

	public final static int CYCLIC_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK = 5;

	public final static int CYCLIC_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 6;

	public final static int SWITCHING_FP_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 7;

	public final static int SWITCHING_BLOCK_RELAXATION_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 8;

	public final static int SWITCHING_GREEDY_MM_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 9;

	public final static int SPARSE_CROSSING_POISSON_TWO_MODE_NETWORK = 10;

	public final static int SPARSE_CROSSING_INFLATED_ZERO_POISSON_TWO_MODE_NETWORK = 11;

	public final static int DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 12;

	public final static int CYCLIC_DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 13;

	public final static int COL_VARYING_COVARIATES_DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK = 14;

	public final static int SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK = 15;

	public final static int COVARIATE_SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK = 16;

	public final static int BAYESIAN_SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK = 17;
	
	/* The minimum value of a latent membership variable */
	protected double minMembership = 1e-100;

	/* The minimum value of a location parameter if relevant */
	protected double minLocation = 1e-100;

	/* The minimum value of a scale parameter if relevant */
	protected double minScale = 1e-100;

	/* The random generator */
	protected RandomStream randomGenerator;

	/* Some parameters for gradient and Newton-Raphson optimization methods */
	protected int maxIterations = 100;
	protected int maxEvaluations = 100;
	protected double relativeThreshold = 1e-6;
	/*
	 * In order to perform only relative checks, the absolute tolerance must be
	 * set to a negative value
	 */
	protected double absoluteThreshold = -1.0;
	protected ConjugateGradientFormula gradientFormula = ConjugateGradientFormula.FLETCHER_REEVES;

	/* Finalize all model parameters and cached data structures */
	public abstract void finalizeModel();

	/* Getters and Setters */

	public double getMinMembership() {
		return minMembership;
	}

	public void setMinMembership(double minMembership) {
		this.minMembership = minMembership;
	}

	public double getMinLocation() {
		return minLocation;
	}

	public void setMinLocation(double minLocation) {
		this.minLocation = minLocation;
	}

	public double getMinScale() {
		return minScale;
	}

	public void setMinScale(double minScale) {
		this.minScale = minScale;
	}

	public RandomStream getRandomGenerator() {
		return randomGenerator;
	}

	public void setRandomGenerator(RandomStream randomGenerator) {
		this.randomGenerator = randomGenerator;
	}

	public void setRowBlockSize(int rowBlockSize) {
	}

	public void setColumnBlockSize(int columnBlockSize) {
	}

	public void setInterleaving(boolean isInterleaving) {
	}

	public void setThresholdIterations(int thresholdIterations) {
	}

	public void setGreedyStep(double greedyStep) {
	}

}
