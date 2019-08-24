package em;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.TimeUnit;

import model.BayesianTwoModeNetworkModel;
import model.NetworkModel;
import model.TwoModeNetworkModel;
import model.binary.additive.CovariateAdditiveBinaryTwoModeNetworkModel;
import model.binary.additive.CyclicNoCovariateAdditiveBinaryTwoModeNetworkModel;
import model.binary.additive.NoCovariateAdditiveBinaryTwoModeNetworkModel;
import model.binary.additive.RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkModel;
import model.binary.crossing.CovariateCrossingBinaryTwoModeNetworkModel;
import model.binary.crossing.NoCovariateCrossingBinaryFrequentistTwoModeNetworkModel;
import model.binary.crossing.bayes.NoCovariateCrossingBinaryBayesianTwoModeNetworkModel;
import model.count.crossing.NoCovariateCrossingSparseCountTwoModeNetworkModel;
import model.count.crossing.ZeroInflatedNoCovariateCrossingSparseCountTwoModeNetworkModel;
import model.gauss.crossing.diffVar.ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.diffVar.CyclicNoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.diffVar.NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.CovariateCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.NoCovariateCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.SwitchingBlockRelaxationNoCovariateCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.SwitchingCyclicFPNoCovariateCrossingDenseGaussTwoModeNetworkModel;
import model.gauss.crossing.sameVar.SwitchingGreedyMMNoCovariateCrossingDenseGaussTwoModeNetworkModel;

public class TwoModeEMEngine extends EMEngine {

	/* Do we need a random initialization */
	protected boolean isRandomInit = true;

	/* The model requests to run combined EM mode */
	protected boolean isCyclicEM = false;

	/* The number of latent row groups to detect */
	protected int numOfRowGroups;

	/* The number of latent column groups to detect */
	protected int numOfColumnGroups;

	/* The 2-mode model */
	protected TwoModeNetworkModel model;

	public TwoModeEMEngine() {
		super();
	}

	@Override
	public void configure(String configurationFilePath, String delim) {

		/*
		 * Read the configuration file and put pairs of parameters and values
		 * into a Map
		 */
		System.out
				.println("The configuration file is " + configurationFilePath);

		HashMap<String, String> configurationMap = new HashMap<String, String>();
		String line = null;
		try {
			BufferedReader configurationReader = new BufferedReader(
					new FileReader(new File(configurationFilePath)));
			while ((line = configurationReader.readLine()) != null) {
				StringTokenizer tokenizer = new StringTokenizer(line, delim);
				if (tokenizer.countTokens() >= 2) {
					String key = tokenizer.nextToken();
					String value = tokenizer.nextToken();
					configurationMap.put(key, value);
				}
			}
			configurationReader.close();
		} catch (Exception e) {
			System.out
					.println("I/O errors when reading the configuration file "
							+ configurationFilePath);
			e.printStackTrace();
		}

		/* Print out the configuration map for verification */
		for (Iterator<Map.Entry<String, String>> it = configurationMap
				.entrySet().iterator(); it.hasNext();) {
			Map.Entry<String, String> entry = it.next();
			System.out.println(entry.getKey() + " : " + entry.getValue());
		}

		// Read modelType and create the corresponding model
		if (configurationMap.get("modelType") != null) {
			int modelType = Integer.parseInt(configurationMap.get("modelType")
					.trim());
			switch (modelType) {
			case NetworkModel.SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating SparseAdditiveBinaryTwoModeNetwork");
				model = new NoCovariateAdditiveBinaryTwoModeNetworkModel();
				break;
			case NetworkModel.DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating NoCovariateCrossedDenseGaussTwoModeNetworkModel");
				model = new NoCovariateCrossingDenseGaussTwoModeNetworkModel();
				break;
			case NetworkModel.COVARIATE_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating CovariateAdditiveBinaryTwoModeNetworkModel");
				model = new CovariateAdditiveBinaryTwoModeNetworkModel();
				break;
			case NetworkModel.COVARIATE_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating CovariateCrossingDenseGaussTwoModeNetworkModel");
				model = new CovariateCrossingDenseGaussTwoModeNetworkModel();
				break;
			case NetworkModel.ROW_COL_VARYING_COVARIATES_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating RowColCovariatesAdditiveBinaryTwoModeNetworkModel");
				model = new RowColVaryingCovariatesAdditiveBinaryTwoModeNetworkModel();
				break;
			case NetworkModel.CYCLIC_SPARSE_ADDITIVE_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating CyclicNoCovariateAdditiveBinaryTwoModeNetworkModel");
				model = new CyclicNoCovariateAdditiveBinaryTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.CYCLIC_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel");
				model = new CyclicNoCovariateCrossingDenseGaussTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.SWITCHING_FP_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating SwitchingCyclicFPNoCovariateCrossingDenseGaussTwoModeNetworkModel");
				model = new SwitchingCyclicFPNoCovariateCrossingDenseGaussTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.SWITCHING_BLOCK_RELAXATION_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating SwitchingBlockRelaxationNoCovariateCrossingDenseGaussTwoModeNetworkModel");
				model = new SwitchingBlockRelaxationNoCovariateCrossingDenseGaussTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.SWITCHING_GREEDY_MM_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating SwitchingGreedyMMNoCovariateCrossingDenseGaussTwoModeNetworkModel");
				model = new SwitchingGreedyMMNoCovariateCrossingDenseGaussTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.SPARSE_CROSSING_POISSON_TWO_MODE_NETWORK:
				System.out
						.println("Creating NoCovariateCrossingSparseCountTwoModeNetworkModel");
				model = new NoCovariateCrossingSparseCountTwoModeNetworkModel();
				break;
			case NetworkModel.SPARSE_CROSSING_INFLATED_ZERO_POISSON_TWO_MODE_NETWORK:
				System.out
						.println("Creating ZeroInflatedNoCovariateCrossingSparseCountTwoModeNetworkModel");
				model = new ZeroInflatedNoCovariateCrossingSparseCountTwoModeNetworkModel();
				break;
			case NetworkModel.DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel");
				model = new NoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel();
				break;
			case NetworkModel.CYCLIC_DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating CyclicNoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel");
				model = new CyclicNoCovariateDifferentVarCrossingDenseGaussTwoModeNetworkModel();
				isCyclicEM = true;
				break;
			case NetworkModel.COL_VARYING_COVARIATES_DIFFERENT_VAR_DENSE_CROSSING_GAUSS_TWO_MODE_NETWORK:
				System.out
						.println("Creating ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel");
				model = new ColVaryingCovariatesDifferentVarCrossingDenseGaussTwoModeNetworkModel();
				break;
			case NetworkModel.SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating NoCovariateCrossingBinaryTwoModeNetworkModel");
				model = new NoCovariateCrossingBinaryFrequentistTwoModeNetworkModel();
				break;
			case NetworkModel.COVARIATE_SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating CovariateCrossingBinaryTwoModeNetworkModel");
				model = new CovariateCrossingBinaryTwoModeNetworkModel();
				break;
			case NetworkModel.BAYESIAN_SPARSE_CROSSING_BINARY_TWO_MODE_NETWORK:
				System.out
						.println("Creating NoCovariateCrossingBinaryBayesianTwoModeNetworkModel");
				model = new NoCovariateCrossingBinaryBayesianTwoModeNetworkModel();
				break;
			default:
				System.out.println("You must specify the network model type!");
				System.exit(-1);
			}
		} else {
			System.out.println("You must specify the network model type!");
			System.exit(-1);
		}

		/* EM parameters */

		/* The number of random starting EM trials to runF */
		if (configurationMap.get("numOfEMRuns") != null)
			setNumOfEMRuns(Integer
					.parseInt(configurationMap.get("numOfEMRuns")));

		/* The number of EM iterations in each EM run */
		if (configurationMap.get("maxEMIterations") != null)
			setMaxEMIterations(Integer.parseInt(configurationMap
					.get("maxEMIterations")));

		/*
		 * The number of E step to run if membership latent variables have not
		 * been converged. This is only relevanumOfColCovariatesnt for the
		 * fixed-point equation method.
		 */
		if (configurationMap.get("maxESteps") != null)
			setMaxESteps(Integer.parseInt(configurationMap.get("maxESteps")));

		/* Which method updates, MM or FP, are used in E step */
		if (configurationMap.get("isMMinEStep") != null)
			if (Integer.parseInt(configurationMap.get("isMMinEStep")) != 0)
				isMMinEStep = true;
			else
				isMMinEStep = false;

		/* Which method updates, MM or FP, are used in M step */
		if (configurationMap.get("isMMinMStep") != null)
			if (Integer.parseInt(configurationMap.get("isMMinMStep")) != 0)
				isMMinMStep = true;
			else
				isMMinMStep = false;

		/*
		 * The E step is stopped if no difference in latent variables between
		 * two E-step updates is greater than this value
		 */
		if (configurationMap.get("membershipPrecision") != null)
			setMembershipPrecision(Double.parseDouble(configurationMap
					.get("membershipPrecision")));

		/*
		 * The M step is stopped if no relative change in model parameters is
		 * greater than this value
		 */
		if (configurationMap.get("relativeParameterPrecision") != null)
			setRelativeParameterPrecision(Double.parseDouble(configurationMap
					.get("relativeParameterPrecision")));

		/*
		 * The current EM run is stopped if the relative change in the lower
		 * bound is smaller than this valuenumOfColCovariates
		 */
		if (configurationMap.get("relativeLowerBoundPrecision") != null)
			setRelativeLowerBoundPrecision(Double.parseDouble(configurationMap
					.get("relativeLowerBoundPrecision")));

		/* The seed for the random generator */
		int[] myRandomSeeds = new int[4];
		if (configurationMap.get("randomSeed") != null) {
			StringTokenizer seedNumbers = new StringTokenizer(
					configurationMap.get("randomSeed"), ";");
			for (int k = 0; k < myRandomSeeds.length; k++)
				myRandomSeeds[k] = Integer.parseInt(seedNumbers.nextToken()
						.trim());
			setRandomSeeds(myRandomSeeds);
			model.setRandomGenerator(randomGenerator);
		}

		/* Model parameters */

		if (configurationMap.get("dataFile") != null) {

			String dataFile = configurationMap.get("dataFile").trim();

			int dataFormat = SPARSE_DATA_FORMAT; // Default 1 for sparse and 0
													// for matrix
			if (configurationMap.get("dataFormat") != null)
				dataFormat = Integer.parseInt(configurationMap
						.get("dataFormat").trim());

			String dataDelim = "\t"; // Default 1 for tab and 0 for space
			if (configurationMap.get("dataDelim") != null
					&& (Integer.parseInt(configurationMap.get("dataDelim")
							.trim()) == SPACE_DELIM))
				dataDelim = " ";

			String rowCovFile = configurationMap.get("rowCovFile");

			String colCovFile = configurationMap.get("colCovFile");

			String edgeCovFile = configurationMap.get("edgeCovFile");
			model.setData(dataFile, rowCovFile, colCovFile, edgeCovFile,
					dataFormat, dataDelim);

		}

		/* The number of latent row groups to detect */
		if (configurationMap.get("numOfRowGroups") != null) {
			setNumOfRowGroups(Integer.parseInt(configurationMap
					.get("numOfRowGroups")));
		}

		/* The number of latent column groups to detect */
		if (configurationMap.get("numOfColumnGroups") != null) {
			setNumOfColumnGroups(Integer.parseInt(configurationMap
					.get("numOfColumnGroups")));
		}

		/* The minimum value of a latent membership variable */
		if (configurationMap.get("minMembership") != null)
			setMinMembership(Double.parseDouble(configurationMap
					.get("minMembership")));

		/* The minimum value of a location parameter if relevant */
		if (configurationMap.get("minLocation") != null)
			setMinLocation(Double.parseDouble(configurationMap
					.get("minLocation")));

		/* The minimum value of a scale parameter if relevant */
		if (configurationMap.get("minScale") != null)
			setMinScale(Double.parseDouble(configurationMap.get("minScale")));

		if (isCyclicEM) {
			if (configurationMap.get("rowBlockSize") != null) {
				int rowBlockSize = Integer.parseInt(configurationMap
						.get("rowBlockSize"));
				model.setRowBlockSize(rowBlockSize);
				System.out.println("rowBlockSize = " + rowBlockSize);
			}
			if (configurationMap.get("columnBlockSize") != null) {
				int columnBlockSize = Integer.parseInt(configurationMap
						.get("columnBlockSize"));
				model.setColumnBlockSize(columnBlockSize);
				System.out.println("columnBlockSize = " + columnBlockSize);
			}
			if (configurationMap.get("isInterleaving") != null) {
				int isInterleaving = Integer.parseInt(configurationMap
						.get("isInterleaving"));
				model.setInterleaving((isInterleaving != 0));
				System.out.println("isInterleaving = " + (isInterleaving != 0));
			}
			if (configurationMap.get("thresholdIterations") != null) {
				int thresholdIterations = Integer.parseInt(configurationMap
						.get("thresholdIterations"));
				model.setThresholdIterations(thresholdIterations);
				System.out.println("thresholdIterations = "
						+ thresholdIterations);
			}

			if (configurationMap.get("greedyStep") != null) {
				double greedyStep = Double.parseDouble(configurationMap
						.get("greedyStep"));
				model.setGreedyStep(greedyStep);
				System.out.println("greedyStep = " + greedyStep);
			}

		}

		/* The minimum value of a scale parameter if relevant */
		if (configurationMap.get("isRandomInit") != null) {
			isRandomInit = (Integer.parseInt(configurationMap
					.get("isRandomInit")) != 0);
			System.out.println("isRandomInit = " + isRandomInit);
		}

		if (model instanceof BayesianTwoModeNetworkModel)
			((BayesianTwoModeNetworkModel) model)
					.setHyperParameters(configurationMap);

	}

	public void setDataFiles(int dataFormat, String dataDelim, String dataFile,
			String rowCovFile, String colCovFile, String edgeCovFile) {
		model.setData(dataFile, rowCovFile, colCovFile, edgeCovFile,
				dataFormat, dataDelim);
	}

	@Override
	/* Finalize model data and parameters */
	public void finalizeModel() {
		System.out.println(this.getClass().getName() + "::finalizeModel()");
		model.finalizeModel();
	}

	@Override
	/*
	 * Run the EM procedure without logging any information of the run
	 */
	public void runEM() {
		System.out.println(this.getClass().getName()
				+ "::runEM() WITHOU logging information");

		// Assuming that model has been created
		// data and all parameters have beerunEMIterationn configured
		model.reset(isRandomInit, isMMinMStep);

		long startTime = System.nanoTime();

		System.out.println("Starting a EM run");
		for (int emIteration = 0; emIteration < maxEMIterations; emIteration++) {

			System.out.print("EM Iteration " + emIteration + "\t");
			System.out.println(TimeUnit.NANOSECONDS.toSeconds(System.nanoTime()
					- startTime));

			// E Step
			for (int eStep = 0; eStep < maxESteps; eStep++) {
				System.out.println("E Step " + eStep);
				model.runEStep(emIteration, eStep, isMMinEStep);
				if (!model
						.areMembershipVariablesChangedSignificantly(membershipPrecision)) {
					System.out
							.println("Membership parameters are not significantly changed!");
					break;
				}
			}

			// M Step
			model.runMStep(emIteration, isMMinMStep);
			if (!model
					.areOtherModelParametersChangedSignificantly(relativeParameterPrecision)) {
				System.out
						.println("Model parameters are not significantly changed!");
				break;
			}
		}
	}

	@Override
	/*
	 * Run the EM procedure and log all information of the run
	 */
	public void runEM(PrintStream outputStream, String formatString) {
		if (isCyclicEM) {
			System.out.println("Running in the combined EM mode!!!");
			runCombinedEM(outputStream, formatString);
		} else {
			System.out.println("Running in the separate EM mode!!!");
			runSeparateEM(outputStream, formatString);
		}
	}

	/* The non-combined EM method where E and M steps are run separately */
	public void runSeparateEM(PrintStream outputStream, String formatString) {

		outputStream.println(this.getClass().getName()
				+ "::runDetailedEM() WITH logging information");

		// Assuming that model has been created
		// Data and all parameters have been configured
		model.reset(isRandomInit, isMMinMStep);

		double initialLLB = model.getLikelihoodLowerBound();
		long startTime = System.nanoTime();
		outputStream.println("LowerBound:");
		outputStream.print(0 + "\t");
		outputStream.print(0 + "\t");
		outputStream.format(formatString, initialLLB);
		outputStream.print("\t");
		outputStream.format(formatString, 0.0);
		outputStream.println();

		outputStream.println("Starting a EM run");
		double lastLLHLowerBound = 0;
		for (int emIteration = 0; emIteration < maxEMIterations; emIteration++) {

			outputStream.println("EM Iteration " + emIteration);

			// E Step
			for (int eStep = 0; eStep < maxESteps; eStep++) {
				outputStream.println("E Step " + eStep);
				model.runEStep(emIteration, eStep, isMMinEStep);
				if (isMMinEStep
						|| !model
								.areMembershipVariablesChangedSignificantly(membershipPrecision)) {
					outputStream
							.println("Membership parameters are not significantly changed!");
					break;
				}
			}

			// System.out.println("E-step Lower Bound\t"
			// + model.getLikelihoodLowerBound());

			// M Step
			model.runMStep(emIteration, isMMinMStep);
			if (!model
					.areOtherModelParametersChangedSignificantly(relativeParameterPrecision)) {
				outputStream
						.println("Model parameters are not significantly changed!");
				break;
			}

			// System.out.println("M-step Lower Bound\t"
			// + model.getLikelihoodLowerBound());

			// Print out some convergence details
			double newLLHLowerBound = model.getLikelihoodLowerBound();
			double largestTauChange = model.getLargestMembershipChange();

			outputStream.println("LowerBound:");
			outputStream.print((emIteration + 1) + "\t");
			outputStream.print(TimeUnit.NANOSECONDS.toMillis(System.nanoTime()
					- startTime)
					+ "\t");
			outputStream.format(formatString, newLLHLowerBound);
			outputStream.print("\t");
			outputStream.format(formatString, largestTauChange);
			outputStream.println();

			// Check the convergence criterion
			if ((emIteration > 0)
					&& (Math.abs((newLLHLowerBound - lastLLHLowerBound)
							/ newLLHLowerBound) < relativeLowerBoundPrecision)) {
				outputStream
						.println("The EM run stopped since the convergence criterion is satistfied!");
				break;
			} else
				lastLLHLowerBound = newLLHLowerBound;
		}

		model.printModel(System.out, formatString);
	}

	public double[] computeRowMSEs(double[] rowBeta) {
		return model.computeRowMSEs(rowBeta);
	}

	public double[] computeMisRates(String outputZFilename,
			String outputWFilename) {
		return model.computeMisRates(outputZFilename, outputWFilename);
	}
	
	public double computeBIC() {
		return model.computeBIC();
	}

	/* The non-combined EM method where E and M steps are run separately */
	public void runCombinedEM(PrintStream outputStream, String formatString) {

		outputStream.println(this.getClass().getName()
				+ "::runDetailedEM() WITH logging information");

		// Assuming that model has been created
		// Data and all parameters have been configured
		model.reset(isRandomInit, isMMinMStep);

		double initialLLB = model.getLikelihoodLowerBound();
		long startTime = System.nanoTime();
		outputStream.println("LowerBound:");
		outputStream.print(0 + "\t");
		outputStream.print(0 + "\t");
		outputStream.format(formatString, initialLLB);
		outputStream.print("\t");
		outputStream.format(formatString, 0.0);
		outputStream.println();

		outputStream.println("Starting a EM run");
		double lastLLHLowerBound = 0;
		for (int emIteration = 0; emIteration < maxEMIterations; emIteration++) {

			outputStream.println("EM Iteration " + emIteration);

			model.runEMIteration(emIteration, maxESteps, membershipPrecision,
					isMMinEStep, isMMinMStep);

			// Check the convergence criterion
			if (!model
					.areOtherModelParametersChangedSignificantly(relativeParameterPrecision)) {
				outputStream
						.println("Model parameters are not significantly changed!");
				break;
			}

			// Print out some convergence details
			double newLLHLowerBound = model.getLikelihoodLowerBound();
			double largestTauChange = model.getLargestMembershipChange();

			outputStream.println("LowerBound:");
			outputStream.print((emIteration + 1) + "\t");
			outputStream.print(TimeUnit.NANOSECONDS.toMillis(System.nanoTime()
					- startTime)
					+ "\t");
			outputStream.format(formatString, newLLHLowerBound);
			outputStream.print("\t");
			outputStream.format(formatString, largestTauChange);
			outputStream.println();

			// Check the convergence criterion
			if ((emIteration > 0)
					&& (Math.abs((newLLHLowerBound - lastLLHLowerBound)
							/ newLLHLowerBound) < relativeLowerBoundPrecision)) {
				outputStream
						.println("The EM run stopped since the convergence criterion is satistfied!");
				break;
			} else
				lastLLHLowerBound = newLLHLowerBound;

		}

		model.printModel(System.out, formatString);
	}

	@Override
	/*
	 * Run the EM procedure and only print out the final likelihood lower bound
	 */
	public void findLikelihoodLowerBound(PrintStream outputStream,
			String formatString) {
		// Run the EM procedure WITHOUT logging
		runEM();
		// Print out the final likelihood lower bound
		double lowerBound = model.getLikelihoodLowerBound();
		outputStream.print("LowerBound: ");
		outputStream.format(formatString, lowerBound);
	}

	/* Getters and Setters */

	public int getNumOfRowGroups() {
		return numOfRowGroups;
	}

	public void setNumOfRowGroups(int numOfRowGroups) {
		System.out.println("Number of Row Groups is " + numOfRowGroups);
		this.numOfRowGroups = numOfRowGroups;
		if (model != null)
			model.setNumOfRowGroups(this.numOfRowGroups);
	}

	public int getNumOfColumnGroups() {
		return numOfColumnGroups;
	}

	public void setNumOfColumnGroups(int numOfColumnGroups) {
		System.out.println("Number of Column Groups is " + numOfColumnGroups);
		this.numOfColumnGroups = numOfColumnGroups;
		if (model != null)
			model.setNumOfColumnGroups(this.numOfColumnGroups);
	}

	public void setMinMembership(double minMembership) {
		super.setMinMembership(minMembership);
		if (model != null)
			model.setMinMembership(minMembership);
	}

	public void setMinLocation(double minLocation) {
		this.minLocation = minLocation;
		if (model != null)
			model.setMinLocation(minLocation);
	}

	public void setMinScale(double minScale) {
		this.minScale = minScale;
		if (model != null)
			model.setMinScale(minScale);
	}

}
