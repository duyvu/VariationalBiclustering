package graph;

import exception.UnsupportedMethodException;

public abstract class Graph {

	protected int numOfRows;

	protected int numOfColumns;

	public Graph() {
	}

	public abstract void getMatrix(double[][] matrix)
			throws UnsupportedMethodException;

	public abstract void setMatrix(double[][] matrix)
			throws UnsupportedMethodException;

	public abstract double getElement(int i, int j)
			throws UnsupportedMethodException;

	public abstract void setElement(int i, int j, double value)
			throws UnsupportedMethodException;

	public int getNumOfRows() {
		return numOfRows;
	}

	public void setNumOfRows(int numOfRows) {
		this.numOfRows = numOfRows;
	}

	public int getNumOfColumns() {
		return numOfColumns;
	}

	public void setNumOfColumns(int numOfColumns) {
		this.numOfColumns = numOfColumns;
	}

}
