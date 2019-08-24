package graph;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import exception.UnsupportedMethodException;

public class SparseGraph extends Graph implements SparseAccess {

	LinkedList<MatrixElement> edges;
	HashMap<Integer, LinkedList<ColumnElement>> rows;
	HashMap<Integer, LinkedList<RowElement>> columns;

	protected void createInternalData(int _numOfRows, int _numOfColumns) {
		numOfRows = _numOfRows;
		numOfColumns = _numOfColumns;
		edges = new LinkedList<MatrixElement>();
		rows = new HashMap<Integer, LinkedList<ColumnElement>>();
		columns = new HashMap<Integer, LinkedList<RowElement>>();
	}

	public SparseGraph(int _numOfRows, int _numOfColumns) {
		createInternalData(_numOfRows, _numOfColumns);
	}

	public SparseGraph(int _numOfRows, int _numOfColumns,
			LinkedList<MatrixElement> _edges) {
		createInternalData(_numOfRows, _numOfColumns);
		try {
			for (MatrixElement element : _edges)
				setElement(element.row, element.col, element.value);
		} catch (Exception e) {
			System.out.println("This exception could not happen since"
					+ this.getClass().getName() + " supports setElement()");
		}
	}

	@Override
	public Iterator<MatrixElement> elementIterator() {
		if (edges != null)
			return edges.iterator();
		else
			return null;
	}

	@Override
	public Iterator<ColumnElement> rowIterator(int row) {
		if (rows.get(row) != null)
			return rows.get(row).iterator();
		else
			return null;
	}

	@Override
	public Iterator<RowElement> columnIterator(int col) {
		if (columns.get(col) != null)
			return columns.get(col).iterator();
		else
			return null;
	}

	@Override
	public void getMatrix(double[][] matrix) throws UnsupportedMethodException {
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[0].length; j++)
				matrix[i][j] = 0;
		for (MatrixElement element : edges)
			matrix[element.row][element.col] = element.value;
	}

	@Override
	public void setMatrix(double[][] matrix) throws UnsupportedMethodException {
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[0].length; j++)
				if (matrix[i][j] != 0)
					setElement(i, j, matrix[i][j]);
	}

	@Override
	public double getElement(int i, int j) throws UnsupportedMethodException {
		if (rowIterator(i) != null)
			for (Iterator<ColumnElement> it = rowIterator(i); it.hasNext();) {
				ColumnElement element = it.next();
				if (element.col == j)
					return element.value;
			}
		return 0;
	}

	@Override
	public void setElement(int i, int j, double value)
			throws UnsupportedMethodException {

		edges.add(new MatrixElement(i, j, value));

		if (rows.containsKey(i))
			rows.get(i).add(new ColumnElement(j, value));
		else {
			LinkedList<ColumnElement> newRow = new LinkedList<ColumnElement>();
			newRow.add(new ColumnElement(j, value));
			rows.put(i, newRow);
		}

		if (columns.containsKey(j))
			columns.get(j).add(new RowElement(i, value));
		else {
			LinkedList<RowElement> newColumn = new LinkedList<RowElement>();
			newColumn.add(new RowElement(i, value));
			columns.put(j, newColumn);
		}

	}

	@Override
	public int getNumberOfNonzeroEdges() {
		return edges.size();
	}

}
