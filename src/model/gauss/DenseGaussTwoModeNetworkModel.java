package model.gauss;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.StringTokenizer;

import graph.DenseGraph;
import graph.MatrixElement;
import model.ContinuousFrequentistTwoModeNetworkModel;

public abstract class DenseGaussTwoModeNetworkModel extends
		ContinuousFrequentistTwoModeNetworkModel {

	@Override
	protected void createGraph(int numOfRows, int numOfColumns,
			LinkedList<MatrixElement> edges) {
		graph = new DenseGraph(numOfRows, numOfColumns);
		try {
			for (MatrixElement element : edges)
				graph.setElement(element.getRow(), element.getCol(),
						element.getValue());
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Override
	protected void createGraph(int numOfRows, int numOfColumns, double[][] edges) {
		graph = new DenseGraph(numOfRows, numOfColumns, edges);
	}

	@Override
	protected void readEdgeList(BufferedReader dataReader, String dataDelim,
			int numOfEdges, LinkedList<MatrixElement> edges) {
		try {
			for (int e = 0; e < numOfEdges; e++) {
				StringTokenizer tokenizer = new StringTokenizer(dataReader
						.readLine().trim(), dataDelim);
				int row = Integer.parseInt(tokenizer.nextToken().trim());
				int col = Integer.parseInt(tokenizer.nextToken().trim());
				double value = Double.parseDouble(tokenizer.nextToken().trim());
				edges.add(new MatrixElement(row, col, value));
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	@Override
	protected void readMatrix(BufferedReader dataReader, String dataDelim,
			double[][] edges) {
		try {
			for (int i = 0; i < edges.length; i++) {
				String line = dataReader.readLine().trim();
				StringTokenizer tokenizer = new StringTokenizer(line, dataDelim);
				for (int j = 0; j < edges[0].length; j++)
					edges[i][j] = Double.parseDouble(tokenizer.nextToken()
							.trim());
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

}
