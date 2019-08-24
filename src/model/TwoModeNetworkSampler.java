package model;

import exception.UnsupportedMethodException;
import umontreal.iro.lecuyer.rng.RandomStream;

public interface TwoModeNetworkSampler {

	/* Read model parameters from a file */
	public void setParameters(String modelFileName, String delim);

	/*
	 * Print model parameters to a file which is a way to cross-check results
	 * later
	 */
	public void printModelParameters2File(String confirmedModelFilename,
			String delim);

	/*
	 * Generate a sample and save it to files where there is no covariate on
	 * rows, columns, or edges
	 */
	public void generateSampleWITHOUTCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputYFilename, String delim);

	/*
	 * Generate a sample and save it to files where there are covariates on
	 * rows, columns, or edges. These covariates are randomly generated
	 */
	public void generateSampleWITHRandomCovariates(RandomStream sampler,
			int numOfRows, int numOfColumns, String outputZFilename,
			String outputWFilename, String outputRowCovFilename,
			String outputColCovFilename, String outputEdgeCovFilename,
			String outputYFilename, String delim);

	/*
	 * Generate a sample and save it to files where fixed covariates on rows,
	 * columns, or edges are given.
	 */
	public void generateSampleWITHFixedCovariates(RandomStream sampler,
			String outputZFilename, String outputWFilename,
			String outputYFilename, String delim)
			throws UnsupportedMethodException;

}
