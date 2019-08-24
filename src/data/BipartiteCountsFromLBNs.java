package data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.TreeSet;
import java.util.Vector;

public class BipartiteCountsFromLBNs {

	public static int lineThreshold = 100000;

	public static double ONE_DATE = 24 * 60 * 60 * 1000;

	public static DateFormat dateFormatter = new SimpleDateFormat(
			"yyyy-MM-dd'T'HH:mm:ss'Z'");

	public static HashMap<String, Integer> locationIDMap = new HashMap<String, Integer>();
	public static int currentLocationIDIndex = 0;

	public static HashMap<Integer, Coordinate> locationCoordinateMap = new HashMap<Integer, Coordinate>();

	public static Vector<Integer> locationEventCounter = new Vector<Integer>();

	public static int numOfUsers = 58228;

	public static int[] userEventCounter = new int[numOfUsers];

	public static Vector<HashMap<Integer, Integer>> usersLocations = new Vector<HashMap<Integer, Integer>>();

	public static void main(String[] args) throws Exception {

		String checkInEvents = "/media/OS/workspace/data/location-based/Brightkite_totalCheckins.txt";

		readDataFile(checkInEvents);

		outputData(50, 50);

		outputData(100, 100);

	}

	protected static void outputData(int userThreshold, int locationThreshold)
			throws IOException {
		String outputFile = "data/poisson/brightkite/Y_"
				+ Integer.toString(userThreshold) + "_"
				+ Integer.toString(locationThreshold) + ".txt";

		PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(
				new File(outputFile))));

		TreeSet<Integer> selectedLocations = new TreeSet<Integer>();
		HashMap<Integer, Integer> newLocationIDMap = new HashMap<Integer, Integer>();
		int currentNewLocationIndex = 0;
		for (int j = 0; j < currentLocationIDIndex; j++)
			if (locationEventCounter.get(j) >= locationThreshold) {
				selectedLocations.add(j);
				newLocationIDMap.put(j, currentNewLocationIndex++);
			}

		int currentNewUserIndex = 0;
		for (int i = 0; i < numOfUsers; i++) {
			if (userEventCounter[i] >= userThreshold) {
				HashMap<Integer, Integer> myMap = usersLocations.get(i);
				for (int location : selectedLocations) {
					if (myMap.containsKey(location))
						writer.write(currentNewUserIndex + "\t"
								+ newLocationIDMap.get(location) + "\t"
								+ myMap.get(location) + "\n");
				}
				currentNewUserIndex++;
			}
		}

		writer.close();

		System.out.println(currentNewUserIndex);
		System.out.println(selectedLocations.size());

	}

	protected static void readDataFile(String checkInEvents) {

		for (int i = 0; i < numOfUsers; i++)
			usersLocations.add(new HashMap<Integer, Integer>());

		try {
			System.out.println("Reading check-in events!!!");
			BufferedReader userEnterEventsReader = new BufferedReader(
					new FileReader(new File(checkInEvents)));

			// Read users and enter times
			String line = null;
			int lineCounter = 0;
			while ((line = userEnterEventsReader.readLine()) != null) {

				lineCounter++;
				if (lineCounter % lineThreshold == 0)
					System.out.println("Reading up to " + lineCounter);

				StringTokenizer myTokenizer = new StringTokenizer(line, "\t");

				int userID = Integer.parseInt(myTokenizer.nextToken().trim());
				userEventCounter[userID]++;

				double eventTime = dateFormatter.parse(
						myTokenizer.nextToken().trim()).getTime()
						/ ONE_DATE;

				double latitude = Double.parseDouble(myTokenizer.nextToken()
						.trim());

				double longitude = Double.parseDouble(myTokenizer.nextToken()
						.trim());

				String locationID = myTokenizer.nextToken().trim();
				int mappedLocationID;
				if (locationIDMap.containsKey(locationID)) {
					mappedLocationID = locationIDMap.get(locationID);
					locationEventCounter.set(mappedLocationID,
							locationEventCounter.get(mappedLocationID) + 1);

				} else {
					mappedLocationID = currentLocationIDIndex++;
					locationIDMap.put(locationID, mappedLocationID);
					locationCoordinateMap.put(mappedLocationID, new Coordinate(
							latitude, longitude));
					locationEventCounter.add(1);
				}

				HashMap<Integer, Integer> myMap = usersLocations.get(userID);
				if (myMap.containsKey(mappedLocationID))
					myMap.put(mappedLocationID, myMap.get(mappedLocationID) + 1);
				else
					myMap.put(mappedLocationID, 1);

				System.out.println(userID + "\t" + eventTime);
			}
			userEnterEventsReader.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Errors in reading user enter events: "
					+ checkInEvents);
		}
	}
}
