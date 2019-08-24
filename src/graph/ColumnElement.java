package graph;

public class ColumnElement {

	protected int col;
	protected double value;

	public int getCol() {
		return col;
	}

	public double getValue() {
		return value;
	}

	public ColumnElement(int _col, double _value) {
		col = _col;
		value = _value;
	}

	public boolean equals(Object object) {
		ColumnElement anotherMatrixElement = (ColumnElement) object;
		if (col == anotherMatrixElement.col
				&& value == anotherMatrixElement.value)
			return true;
		else
			return false;
	}

	public int hashCode() {
		return (new Integer(col)).hashCode() + (new Double(value).hashCode());
	}

}
