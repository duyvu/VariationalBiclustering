package model.binary;

import java.util.LinkedList;

import graph.MatrixElement;
import graph.SparseGraph;
import model.DiscreteFrequentistTwoModeNetworkModel;

public abstract class SparseBinaryFrequentistTwoModeNetworkModel extends
		DiscreteFrequentistTwoModeNetworkModel {

	/*
	 * Create the data graph which depends on the subclasses. If the subclasses
	 * are sparse models, sparse graphs are created
	 */
	protected void createGraph(int numOfRows, int numOfColumns,
			LinkedList<MatrixElement> edges) {
		graph = new SparseGraph(numOfRows, numOfColumns, edges);
	}

	/*
	 * Create the data graph which depends on the subclasses. If the subclasses
	 * are sparse models, sparse graphs are created
	 */
	protected void createGraph(int numOfRows, int numOfColumns, double[][] edges) {
		graph = new SparseGraph(numOfRows, numOfColumns);
		try {
			graph.setMatrix(edges);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

}
