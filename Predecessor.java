import java.util.*;

public class Predecessor {

    // Compute list of predecessors of a currentCoordinate
    public static ArrayList<int[]> getPredecessors(int[] currentCoordinate, int diag) {

        ArrayList<int[]> listPredecessors = new ArrayList<>(7);

        // Start from 1, because an empty M Index set is not allowed

        outer: for (int i = 1; i < (1 << currentCoordinate.length); i++) {
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

            if (diag >= 0) {
                int minCoord = Integer.MAX_VALUE;
                int maxCoord = 0;
                for (int c : coordsPredecessor) {
                    if (c < minCoord)
                        minCoord = c;
                    if (c > maxCoord)
                        maxCoord = c;
                }
                if (maxCoord - minCoord > diag) {
                    continue outer;
                }
            }

            listPredecessors.add(coordsPredecessor);
        }
        return listPredecessors;

    }

    public static void main(String[] args) {
        int[] currentCoords = new int[] { 1, 2, 3 };
        ArrayList<int[]> pred = Predecessor.getPredecessors(currentCoords, -1);
        if (pred.isEmpty())
            System.out.println("empty");
        int c = 0;
        for (int[] p : pred) {
            c++;
            System.out.println(Arrays.toString(p) + " " + c);
        }
    }

}
