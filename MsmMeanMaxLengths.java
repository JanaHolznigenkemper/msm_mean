import java.util.ArrayList;
import java.util.Arrays;

// class extending MSMMean for computing a mean that is as most as long as the longest time series.
public class MsmMeanMaxLengths extends MsmMean {

    public MsmMeanMaxLengths(double[][] timeseries, double c) {
        super(timeseries, c);
        this.lLengths = this.maxDimension;
        this.loopMultiplier = 1;
    }

    public MsmMeanMaxLengths(double[][] timeseries, double c, int diag) {
        this(timeseries, c);
        this.diag = diag;
    }

    public void calculateTableEntry() {

        setOriginCoordinates();

        // iterate through all sums of coordinates
        for (int i = 1; i <= this.sumLengthTimeseries; i++) {

            ArrayList<int[]> sumICoords = this.sumCoords.get(i);
            for (int[] currentCoord : sumICoords) {

                int flattenCurrentCoord = table.flattenCoordinate(currentCoord);
                double[][] firstDimTableCurrentCoord = table.getFirstDim(flattenCurrentCoord);

                // get min index of coordinate
                int minCoord = Integer.MAX_VALUE;
                int maxCoord = 0;
                for (int c : currentCoord) {
                    if (c < minCoord)
                        minCoord = c;
                    if (c > maxCoord)
                        maxCoord = c;
                }

                ArrayList<int[]> predecessorList = Predecessor.getPredecessors(currentCoord, this.diag);

                for (int l = Math.max(minCoord - 1, 0); l < lLengths + 1; l++) {

                    double[] secondDimTableCurrentCoord = table.getSecondDim(firstDimTableCurrentCoord, l);

                    outer: for (int s = 0; s < this.distinctTsValues.length; s++) {
                        boolean isValid = false;
                        for (int q = 0; q < k; q++) {
                            if (appDist[s][q] <= currentCoord[q]) {
                                isValid = true;
                                break;
                            }
                        }

                        if (!isValid) {
                            continue outer;
                        }

                        // median coordinate = 0. All other timeseries are only able to merge, the move
                        // and split case is not applied.
                        if (l == 0) {
                            double cost = costBorderCoords(currentCoord, l, s);
                            if (cost == Double.POSITIVE_INFINITY)
                                continue outer;
                            table.set(secondDimTableCurrentCoord, s, cost);
                            continue outer;
                        }

                        // normal case:

                        double cost = Double.POSITIVE_INFINITY;

                        for (int[] predecessor : predecessorList) {
                            int flattenCoordsPredecessor = table.flattenCoordinate(predecessor);
                            double[][] firstDimTablePredecessor = table.getFirstDim(flattenCoordsPredecessor);
                            int[] predecessorMOSP = predecessor.clone();
                            double costMOSP = calcMOSPCost(currentCoord, l, s, predecessorMOSP,
                                    firstDimTablePredecessor);
                            int[] predecessorME = predecessor.clone();
                            double costME = calcMECost(currentCoord, l, s, predecessorME, flattenCoordsPredecessor,
                                    firstDimTablePredecessor);

                            cost = Math.min(Math.min(cost, costME), costMOSP);
                        }

                        if (cost == Double.POSITIVE_INFINITY)
                            continue outer;
                        table.set(secondDimTableCurrentCoord, s, cost);

                    }

                }
            }
        }

    }

    public static void main(String[] args) {
        // small example for k=3, c=0.1
        int n = 14;

        double[][] exp1 = new double[3][n];

        for (int i = 0; i < 3; i++) {

            for (int j = i; j < (n + i); j++) {
                exp1[i][j - i] = j;
            }

            System.out.println(Arrays.toString(exp1[i]));
        }
        System.out.println("Start Computing: k=3 and n=" + n);

        MsmMeanMaxLengths msmMean = new MsmMeanMaxLengths(exp1, 0.1);

        Pair<double[], Double> mean = msmMean.calcMean();

        System.out.println("Mean");

        System.out.println(Arrays.toString(mean.getFirst()));
        System.out.println(mean.getSecond());

    }

}
