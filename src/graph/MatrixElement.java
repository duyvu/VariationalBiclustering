package graph;

public class MatrixElement {

	protected int row;
	protected int col;
	protected double value;

	public int getRow() {
		return row;
	}

	public int getCol() {
		return col;
	}

	public double getValue() {
		return value;
	}

	public MatrixElement(int _row, int _col, double _value) {
		row = _row;
		col = _col;
		value = _value;
	}

	public boolean equals(Object object) {
		MatrixElement anotherMatrixElement = (MatrixElement) object;
		if (row == anotherMatrixElement.row && col == anotherMatrixElement.col
				&& value == anotherMatrixElement.value)
			return true;
		else
			return false;
	}

	public int hashCode() {
		return (new Integer(row)).hashCode() + (new Integer(col)).hashCode()
				+ (new Double(value).hashCode());
	}

}
