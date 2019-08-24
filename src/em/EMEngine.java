package em;

import java.io.PrintStream;

import umontreal.iro.lecuyer.rng.LFSR113;
import umontreal.iro.lecuyer.rng.RandomStream;
import utils.ArrayMethods;

public abstract class EMEngine {

	public static final int MATRIX_DATA_FORMAT = 0;
	public static final int SPARSE_DATA_FORMAT = 1;

	public static final int SPACE_DELIM = 0;
	public static final int TAB_DELIM = 1;

	/* The number of random starting EM trials to runF */
	protected int numOfEMRuns;

	/* The number of EM iterations in each EM run */
	protected int maxEMIterations;

	/*
	 * The number of E step to run if membership latent variables have not been
	 * converged. This is only relevant for the fixed-point equation method.
	 */
	protected int maxESteps;

	/* Which method updates, MM or FP, are used in E step */
	protected boolean isMMinEStep;

	/* Which method updates, MM or FP, are used in M step */
	protected boolean isMMinMStep;

	/*
	 * The E step is stopped if no difference in latent variables between two
	 * E-step updates is greater than this value
	 */
	protected double membershipPrecision;

	/*
	 * The M step is stopped if no relative change in model parameters is
	 * greater than this value
	 */
	protected double relativeParameterPrecision;

	/*
	 * The current EM run is stopped if the relative change in the lower bound
	 * is smaller than this value
	 */
	protected double relativeLowerBoundPrecision;

	/* The minimum value of a latent membership variable */
	protected double minMembership;

	/* The minimum value of a location parameter if relevant */
	protected double minLocation;

	/* The minimum value of a scale parameter if relevant */
	protected double minScale;

	/* The seed for the random generator */
	protected int[] randomSeeds;

	/* The random generator */
	protected RandomStream randomGenerator;

	public EMEngine() {
		randomSeeds = new int[4];
		for (int k = 0; k < randomSeeds.length; k++)
			randomSeeds[k] = 12345;
		randomGenerator = new LFSR113();
	}

	/* Read the configuration file and initialize the EM engine and the model */
	public abstract void configure(String configurationFilePath, String delim);

	/* Finalize model data and parameters */
	public abstract void finalizeModel();

	/*
	 * Run the EM procedure without logging any information of the run
	 */
	public abstract void runEM();

	/*
	 * Run the EM procedure and log all information of the run
	 */
	public abstract void runEM(PrintStream outputStream, String formatString);

	/*
	 * Run the EM procedure and only print out the final likelihood lower bound
	 */
	public abstract void findLikelihoodLowerBound(PrintStream outputStream,
			String formatString);

	/* Getters and Setters */

	public int getNumOfEMRuns() {
		return numOfEMRuns;
	}

	public void setNumOfEMRuns(int numOfEMRuns) {
		this.numOfEMRuns = numOfEMRuns;
	}

	public int getMaxEMIterations() {
		return maxEMIterations;
	}

	public void setMaxEMIterations(int maxEMIterations) {
		this.maxEMIterations = maxEMIterations;
	}

	public int getMaxESteps() {
		return maxESteps;
	}

	public void setMaxESteps(int maxESteps) {
		this.maxESteps = maxESteps;
	}

	public boolean isMMinEStep() {
		return isMMinEStep;
	}

	public void setMMinEStep(boolean isMMinEStep) {
		this.isMMinEStep = isMMinEStep;
	}

	public boolean isMMinMStep() {
		return isMMinMStep;
	}

	public void setMMinMStep(boolean isMMinMStep) {
		this.isMMinMStep = isMMinMStep;
	}

	public double getMembershipPrecision() {
		return membershipPrecision;
	}

	public void setMembershipPrecision(double membershipPrecision) {
		this.membershipPrecision = membershipPrecision;
	}

	public double getRelativeParameterPrecision() {
		return relativeParameterPrecision;
	}

	public void setRelativeParameterPrecision(double relativeParameterPrecision) {
		this.relativeParameterPrecision = relativeParameterPrecision;
	}

	public double getRelativeLowerBoundPrecision() {
		return relativeLowerBoundPrecision;
	}

	public void setRelativeLowerBoundPrecision(
			double relativeLowerBoundPrecision) {
		this.relativeLowerBoundPrecision = relativeLowerBoundPrecision;
	}

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

	public int[] getRandomSeeds() {
		return randomSeeds;
	}

	public void setRandomSeeds(int[] randomSeeds) {
		ArrayMethods.copyIntegerArray(randomSeeds, this.randomSeeds);
		((LFSR113) randomGenerator).setSeed(this.randomSeeds);
	}

	public void setRandomSeed(int seedValue, int seedIndex) {
		if (seedIndex < 4) {
			randomSeeds[seedIndex] = seedValue;
			((LFSR113) randomGenerator).setSeed(randomSeeds);
		} else
			System.out.println("Seed index must be from 0 to 3!");
	}

}
