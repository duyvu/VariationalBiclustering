package model.count;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.StringTokenizer;

import graph.MatrixElement;
import graph.SparseGraph;
import model.DiscreteFrequentistTwoModeNetworkModel;

public abstract class SparseCountTwoModeNetworkModel extends
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

	/* Read the sparse data */
	protected void readEdgeList(BufferedReader dataReader, String dataDelim,
			int numOfEdges, LinkedList<MatrixElement> edges) {
		try {
			for (int e = 0; e < numOfEdges; e++) {
				StringTokenizer tokenizer = new StringTokenizer(dataReader
						.readLine().trim(), dataDelim);
				int row = Integer.parseInt(tokenizer.nextToken().trim());
				int col = Integer.parseInt(tokenizer.nextToken().trim());
				int value = Integer.parseInt(tokenizer.nextToken().trim());
				edges.add(new MatrixElement(row, col, value));
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/* Read the matrix data */
	protected void readMatrix(BufferedReader dataReader, String dataDelim,
			double[][] edges) {
		try {
			for (int i = 0; i < edges.length; i++) {
				StringTokenizer tokenizer = new StringTokenizer(dataReader
						.readLine().trim(), dataDelim);
				for (int j = 0; j < edges[0].length; j++)
					edges[i][j] = Integer
							.parseInt(tokenizer.nextToken().trim());
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
