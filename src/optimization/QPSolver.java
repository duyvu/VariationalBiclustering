package optimization;

public class QPSolver {

	public static void solve01Bounds(double[][] A, double[][] S,
			double[][] tau, double precision) {

		if (A.length <= 0)
			return;

		int n = A[0].length;

		boolean[] J_k = new boolean[n];
		boolean[] J_lambda_a = new boolean[n];

		boolean[] J_lambda_k = new boolean[n];
		boolean[] J_lambda_k_a = new boolean[n];
		boolean[] J_lambda_k_b = new boolean[n];

		double lambda_k;

		for (int rowIndex = 0; rowIndex < A.length; rowIndex++) {
			// Step 1:
			for (int j = 0; j < n; j++) {
				J_k[j] = true;
				J_lambda_a[j] = false;
			}

			boolean done = false;
			// We never loop over n times.
			int count = 0;
			while (!done) {
				count++;

				// Step 2:
				double value1 = 0, value2 = 0;
				for (int j = 0; j < n; j++)
					if (J_k[j]) {
						value1 += 1 / A[rowIndex][j];
						value2 += S[rowIndex][j] / A[rowIndex][j];
					}
				lambda_k = (value2 - 2) / value1;

				// Step 3:
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						if (lambda_k >= S[rowIndex][j]) {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = true;
							J_lambda_k_b[j] = false;
						} else if (lambda_k > (-2 * A[rowIndex][j] + S[rowIndex][j])) {
							J_lambda_k[j] = true;
							J_lambda_k_a[j] = false;
							J_lambda_k_b[j] = false;
						} else {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = false;
							J_lambda_k_b[j] = true;
						}
					} else {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = false;
					}
				}

				// Step 4:
				double delta = 0;
				double value3 = 0, value4 = 0;
				boolean J_lambda_k_empty = true;
				for (int j = 0; j < n; j++) {
					delta += (J_lambda_k_b[j]) ? 1 : 0;
					if (J_lambda_k[j]) {
						value3 += S[rowIndex][j] / A[rowIndex][j];
						value4 += 1 / A[rowIndex][j];
						J_lambda_k_empty = false;
					}
				}
				delta += (value3 / 2 - lambda_k * value4 / 2 - 1);

				// Step 5:
				if (Math.abs(delta) < precision || J_lambda_k_empty
						|| count >= n) {
					// Step 8:
					for (int j = 0; j < n; j++) {
						if (J_lambda_a[j] || J_lambda_k_a[j])
							tau[rowIndex][j] = 0;
						else if (J_lambda_k_b[j])
							tau[rowIndex][j] = 1;
						else
							tau[rowIndex][j] = (S[rowIndex][j] - lambda_k)
									/ (2 * A[rowIndex][j]);
					}
					done = true;
				} else if (delta > 0) { // Step 6:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_a[j]) {
							J_lambda_a[j] = true;
							J_k[j] = false;
						}
					}
				} else { // Step 7:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_b[j])
							tau[rowIndex][j] = 1;
						else
							tau[rowIndex][j] = 0;
					}
					// the solution is reached since one variable is set to 1
					done = true;
				}
			}
		}
	}

	public static void solve01Bounds(double[][] A, double[][] S,
			double[][] tau, int fromRow, int toRow, double precision) {

		if (A.length <= 0 || ((toRow - fromRow) != A.length))
			return;

		int n = A[0].length;

		boolean[] J_k = new boolean[n];
		boolean[] J_lambda_a = new boolean[n];

		boolean[] J_lambda_k = new boolean[n];
		boolean[] J_lambda_k_a = new boolean[n];
		boolean[] J_lambda_k_b = new boolean[n];

		double lambda_k;

		for (int rowIndex = 0, i = fromRow; rowIndex < A.length; rowIndex++, i++) {
			// Step 1:
			for (int j = 0; j < n; j++) {
				J_k[j] = true;
				J_lambda_a[j] = false;
			}

			boolean done = false;
			// We never loop over n times.
			int count = 0;
			while (!done) {
				count++;

				// Step 2:
				double value1 = 0, value2 = 0;
				for (int j = 0; j < n; j++)
					if (J_k[j]) {
						value1 += 1 / A[rowIndex][j];
						value2 += S[rowIndex][j] / A[rowIndex][j];
					}
				lambda_k = (value2 - 2) / value1;

				// Step 3:
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						if (lambda_k >= S[rowIndex][j]) {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = true;
							J_lambda_k_b[j] = false;
						} else if (lambda_k > (-2 * A[rowIndex][j] + S[rowIndex][j])) {
							J_lambda_k[j] = true;
							J_lambda_k_a[j] = false;
							J_lambda_k_b[j] = false;
						} else {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = false;
							J_lambda_k_b[j] = true;
						}
					} else {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = false;
					}
				}

				// Step 4:
				double delta = 0;
				double value3 = 0, value4 = 0;
				boolean J_lambda_k_empty = true;
				for (int j = 0; j < n; j++) {
					delta += (J_lambda_k_b[j]) ? 1 : 0;
					if (J_lambda_k[j]) {
						value3 += S[rowIndex][j] / A[rowIndex][j];
						value4 += 1 / A[rowIndex][j];
						J_lambda_k_empty = false;
					}
				}
				delta += (value3 / 2 - lambda_k * value4 / 2 - 1);

				if (Math.abs(delta) < precision || J_lambda_k_empty
						|| count >= n) {
					for (int j = 0; j < n; j++) {
						if (J_lambda_a[j] || J_lambda_k_a[j])
							tau[i][j] = 0;
						else if (J_lambda_k_b[j])
							tau[i][j] = 1;
						else
							tau[i][j] = (S[rowIndex][j] - lambda_k)
									/ (2 * A[rowIndex][j]);
					}
					done = true;
				} else if (delta > 0) {
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_a[j]) {
							J_lambda_a[j] = true;
							J_k[j] = false;
						}
					}
				} else {
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_b[j])
							tau[i][j] = 1;
						else
							tau[i][j] = 0;
					}
					done = true;
				}
			}
		}
	}

	public static void solveBounds(double[][] M, double[][] S, double[][] tau,
			double precision) {

		if (M.length <= 0)
			return;

		int n = M[0].length;

		boolean[] J_k = new boolean[n];
		boolean[] J_lambda_a = new boolean[n];
		boolean[] J_lambda_b = new boolean[n];

		boolean[] J_lambda_k = new boolean[n];
		boolean[] J_lambda_k_a = new boolean[n];
		boolean[] J_lambda_k_b = new boolean[n];

		double lambda_k;

		for (int rowIndex = 0; rowIndex < M.length; rowIndex++) {

			// Step 0:
			double[] a = new double[n];
			double[] b = new double[n];
			double[] delta_tau = new double[n];
			for (int k = 0; k < n; k++) {
				a[k] = -tau[rowIndex][k];
				b[k] = 1 - tau[rowIndex][k];
			}
			double alpha_k = 0;

			// Step 1:
			for (int j = 0; j < n; j++) {
				J_k[j] = true;
				J_lambda_a[j] = false;
				J_lambda_b[j] = false;
			}
			// The test sum_j d_j a_j <= alpha <= sum_j d_j b_j is always
			// satisfied since the left-hand side is always non-positive
			// while the right-hand side is always non-negative and alpha = 0

			boolean done = false;
			// We never loop over n times.
			int count = 0;
			while (!done) {
				count++;

				// Step 2:
				double numerator = 0, denominartor = 0;
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						numerator += 1 / M[rowIndex][j];
						denominartor += S[rowIndex][j] / M[rowIndex][j];
					}
				}
				double value_a = 0;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k_a[j])
						value_a += a[j];
				}
				double value_b = 0;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k_b[j])
						value_b += b[j];
				}
				denominartor += 2 * value_a + 2 * value_b - 2 * alpha_k;
				lambda_k = denominartor / numerator;

				// Step 3:
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						// 2.32
						double check_A = -2 * M[rowIndex][j] * a[j]
								+ S[rowIndex][j];
						if (lambda_k >= check_A) {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = true;
							J_lambda_k_b[j] = false;
						} else {
							double check_B = -2 * M[rowIndex][j] * b[j]
									+ S[rowIndex][j];
							// 2.34
							if (lambda_k > check_B) {
								J_lambda_k[j] = true;
								J_lambda_k_a[j] = false;
								J_lambda_k_b[j] = false;
							} else { // 2.33
								J_lambda_k[j] = false;
								J_lambda_k_a[j] = false;
								J_lambda_k_b[j] = true;
							}
						}
					} else {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = false;
					}
				}

				// Step 4:
				double sdm = 0, d2m = 0;
				double abAlpha = -alpha_k;
				boolean J_lambda_k_empty = true;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k[j]) {
						sdm += S[rowIndex][j] / M[rowIndex][j];
						d2m += 1 / M[rowIndex][j];
						J_lambda_k_empty = false;
					}
					if (J_lambda_k_a[j])
						abAlpha += a[j];
					if (J_lambda_k_b[j])
						abAlpha += b[j];
				}
				double delta = abAlpha + sdm / 2.0 - lambda_k * d2m / 2.0;

				// Step 5:
				if (Math.abs(delta) < precision || J_lambda_k_empty
						|| count >= n) {
					// Step 8:
					for (int j = 0; j < n; j++) {
						if (J_lambda_a[j] || J_lambda_k_a[j])
							delta_tau[j] = a[j];
						else if (J_lambda_b[j] || J_lambda_k_b[j])
							delta_tau[j] = b[j];
						else
							delta_tau[j] = (S[rowIndex][j] - lambda_k)
									/ (2 * M[rowIndex][j]);
					}
					done = true;
				} else if (delta > 0) { // Step 6:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_a[j]) {
							delta_tau[j] = a[j];
							alpha_k -= a[j];
							J_lambda_a[j] = true;
							J_k[j] = false;
						}
					}
				} else { // Step 7:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_b[j]) {
							delta_tau[j] = b[j];
							alpha_k -= b[j];
							J_lambda_b[j] = true;
							J_k[j] = false;
						}
					}
				}
			}

			// Add the increment
			for (int k = 0; k < n; k++)
				tau[rowIndex][k] += delta_tau[k];
		}

	}

	public static void solveBounds(double[][] M, double[][] S, double[][] tau,
			int fromRow, int toRow, double precision) {

		if (M.length <= 0)
			return;

		int n = M[0].length;

		boolean[] J_k = new boolean[n];
		boolean[] J_lambda_a = new boolean[n];
		boolean[] J_lambda_b = new boolean[n];

		boolean[] J_lambda_k = new boolean[n];
		boolean[] J_lambda_k_a = new boolean[n];
		boolean[] J_lambda_k_b = new boolean[n];

		double lambda_k;

		for (int rowIndex = 0, i = fromRow; rowIndex < M.length; rowIndex++, i++) {

			// Step 0:
			double[] a = new double[n];
			double[] b = new double[n];
			double[] delta_tau = new double[n];
			for (int k = 0; k < n; k++) {
				a[k] = -tau[i][k];
				b[k] = 1 - tau[i][k];
			}
			double alpha_k = 0;

			// Step 1:
			for (int j = 0; j < n; j++) {
				J_k[j] = true;
				J_lambda_a[j] = false;
				J_lambda_b[j] = false;
			}
			// The test sum_j d_j a_j <= alpha <= sum_j d_j b_j is always
			// satisfied since the left-hand side is always non-positive
			// while the right-hand side is always non-negative and alpha = 0

			boolean done = false;
			// We never loop over n times.
			int count = 0;
			while (!done) {
				count++;

				// Step 2:
				double numerator = 0, denominartor = 0;
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						numerator += 1 / M[rowIndex][j];
						denominartor += S[rowIndex][j] / M[rowIndex][j];
					}
				}
				double value_a = 0;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k_a[j])
						value_a += a[j];
				}
				double value_b = 0;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k_b[j])
						value_b += b[j];
				}
				denominartor += 2 * value_a + 2 * value_b - 2 * alpha_k;
				lambda_k = denominartor / numerator;

				// Step 3:
				for (int j = 0; j < n; j++) {
					if (J_k[j]) {
						// 2.32
						double check_A = -2 * M[rowIndex][j] * a[j]
								+ S[rowIndex][j];
						if (lambda_k >= check_A) {
							J_lambda_k[j] = false;
							J_lambda_k_a[j] = true;
							J_lambda_k_b[j] = false;
						} else {
							double check_B = -2 * M[rowIndex][j] * b[j]
									+ S[rowIndex][j];
							// 2.34
							if (lambda_k > check_B) {
								J_lambda_k[j] = true;
								J_lambda_k_a[j] = false;
								J_lambda_k_b[j] = false;
							} else { // 2.33
								J_lambda_k[j] = false;
								J_lambda_k_a[j] = false;
								J_lambda_k_b[j] = true;
							}
						}
					} else {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = false;
					}
				}

				// Step 4:
				double sdm = 0, d2m = 0;
				double abAlpha = -alpha_k;
				boolean J_lambda_k_empty = true;
				for (int j = 0; j < n; j++) {
					if (J_lambda_k[j]) {
						sdm += S[rowIndex][j] / M[rowIndex][j];
						d2m += 1 / M[rowIndex][j];
						J_lambda_k_empty = false;
					}
					if (J_lambda_k_a[j])
						abAlpha += a[j];
					if (J_lambda_k_b[j])
						abAlpha += b[j];
				}
				double delta = abAlpha + sdm / 2.0 - lambda_k * d2m / 2.0;

				// Step 5:
				if (Math.abs(delta) < precision || J_lambda_k_empty
						|| count >= n) {
					// Step 8:
					for (int j = 0; j < n; j++) {
						if (J_lambda_a[j] || J_lambda_k_a[j])
							delta_tau[j] = a[j];
						else if (J_lambda_b[j] || J_lambda_k_b[j])
							delta_tau[j] = b[j];
						else
							delta_tau[j] = (S[rowIndex][j] - lambda_k)
									/ (2 * M[rowIndex][j]);
					}
					done = true;
				} else if (delta > 0) { // Step 6:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_a[j]) {
							delta_tau[j] = a[j];
							alpha_k -= a[j];
							J_lambda_a[j] = true;
							J_k[j] = false;
						}
					}
				} else { // Step 7:
					for (int j = 0; j < n; j++) {
						if (J_lambda_k_b[j]) {
							delta_tau[j] = b[j];
							alpha_k -= b[j];
							J_lambda_b[j] = true;
							J_k[j] = false;
						}
					}
				}
			}

			// Add the increment
			for (int k = 0; k < n; k++)
				tau[i][k] += delta_tau[k];
		}

	}

}
