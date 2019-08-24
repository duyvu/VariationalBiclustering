package utils;

import java.io.PrintStream;

public class PrintingMethods {

	public static void printDoubleVector(double[] vector,
			PrintStream outputStream, String delim) {
		for (int k = 0; k < vector.length; k++)
			outputStream.print(vector[k] + delim);
		outputStream.println();
	}

	public static void printDoubleMatrix(double[][] matrix,
			PrintStream outputStream, String delim) {
		for (int k = 0; k < matrix.length; k++) {
			for (int l = 0; l < matrix[0].length; l++)
				outputStream.print(matrix[k][l] + delim);
			outputStream.println();
		}
	}

}
