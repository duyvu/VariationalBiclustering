package graph;

import utils.ArrayMethods;

import exception.UnsupportedMethodException;

public class DenseGraph extends Graph {

	protected double[][] edges;

	protected void createInternalData(int _numOfRows, int _numOfColumns) {
		numOfRows = _numOfRows;
		numOfColumns = _numOfColumns;
		edges = new double[numOfRows][numOfColumns];
	}

	public DenseGraph(int _numOfRows, int _numOfColumns) {
		createInternalData(_numOfRows, _numOfColumns);
	}

	public DenseGraph(int _numOfRows, int _numOfColumns, double[][] _edges) {
		createInternalData(_numOfRows, _numOfColumns);
		ArrayMethods.copyDoubleMatrix(_edges, edges);
	}

	@Override
	public void getMatrix(double[][] matrix) throws UnsupportedMethodException {
		ArrayMethods.copyDoubleMatrix(edges, matrix);
	}

	@Override
	public void setMatrix(double[][] matrix) throws UnsupportedMethodException {
		ArrayMethods.copyDoubleMatrix(matrix, edges);
	}

	@Override
	public double getElement(int i, int j) throws UnsupportedMethodException {
		return edges[i][j];
	}

	@Override
	public void setElement(int i, int j, double value)
			throws UnsupportedMethodException {
		edges[i][j] = value;
	}

}
