package graph;

public class RowElement {

	protected int row;
	protected double value;

	public int getRow() {
		return row;
	}

	public double getValue() {
		return value;
	}

	public RowElement(int _row, double _value) {
		row = _row;
		value = _value;
	}

	public boolean equals(Object object) {
		RowElement anotherMatrixElement = (RowElement) object;
		if (row == anotherMatrixElement.row
				&& value == anotherMatrixElement.value)
			return true;
		else
			return false;
	}

	public int hashCode() {
		return (new Integer(row)).hashCode() + (new Double(value).hashCode());
	}

}
