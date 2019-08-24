package graph;

import java.util.Iterator;

public interface SparseAccess {

	public Iterator<MatrixElement> elementIterator();

	public Iterator<ColumnElement> rowIterator(int row);

	public Iterator<RowElement> columnIterator(int col);

	public int getNumberOfNonzeroEdges();

}
