import java.util.*;

public class Predecessor {


    /**
     * Compute list of predecessors of a currentCoordinate
     *
     * @param currentCoordinate current coordinates of time series
     * @param window window size
     * @return list of predecessor coordinates
     */
    public static ArrayList<int[]> getPredecessors(int[] currentCoordinate, int window) {

        // init list with just some random capacity
        ArrayList<int[]> listPredecessors = new ArrayList<>(7);

        // Start from 1, because an empty M Index set is not allowed
        outer:
        for (int i = 1; i < (1 << currentCoordinate.length); i++) {
            int m = 1;

            int[] coordsPredecessor = currentCoordinate.clone();

            for (int j = 0; j < currentCoordinate.length; j++) {

                if ((i & m) > 0) {
                    // increase coordinate of m index, all others are split coordinates, no increase
                    // necessary
                    if (currentCoordinate[j] == 0) {
                        continue outer;
                    }
                    coordsPredecessor[j] -= 1;
                }
                m = m << 1;
            }

            // discard as predecessor if coordinates do not fit into window size
            if (window >= 0) {
                int minCoord = Integer.MAX_VALUE;
                int maxCoord = 0;
                for (int c : coordsPredecessor) {
                    if (c < minCoord)
                        minCoord = c;
                    if (c > maxCoord)
                        maxCoord = c;
                }
                if (maxCoord - minCoord > window) {
                    continue outer;
                }
            }

            listPredecessors.add(coordsPredecessor);
        }
        return listPredecessors;

    }

}
