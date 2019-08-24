package exp;

import java.util.Random;

public class TestSumMatrix {

	public static void main(String[] args) {
		double[][] matrix = new double[39][];
		Random sampler = new Random(12345);
		double sum = 0;
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[0].length; j++) {
				matrix[i][j] = sampler.nextDouble();
				sum += matrix[i][j];
			}
		System.out.println("The final sum is " + sum);
	}

}
