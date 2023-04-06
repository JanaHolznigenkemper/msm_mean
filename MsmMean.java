import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

public class MsmMean {

    // MultiDimArray for storing all Coordinates of a table with their respective
    // cost
    protected final MultiDimArray table;

    protected final double[][] timeseries;

    protected final int k;

    protected int sumLengthTimeseries = 0;

    protected final double[] distinctTsValues;

    int lLengths;

    int sLengths;

    // max coordinates (length-1) of time series
    protected int[] maxCoordsTS;

    protected int minDimension = Integer.MAX_VALUE;

    protected int maxDimension = 0;

    // Storing all first appearances of a value of V(X) for each time series
    protected final int[][] appDist;

    // Array at each entry there is a list of coordinates that sum up to the index
    ArrayList<ArrayList<int[]>> sumCoords;

    int loopMultiplier;

    int diag = -1;

    double c;

    public MsmMean(double[][] timeseries, double c) {
        this.c = c;
        this.timeseries = timeseries;
        this.k = timeseries.length;
        for (double[] ts : this.timeseries) {
            this.sumLengthTimeseries += ts.length - 1;
            this.maxDimension = Math.max(maxDimension, ts.length);
            this.minDimension = Math.min(minDimension, ts.length);
        }

        // create list of all time series data points
        TreeSet<Double> values = new TreeSet<Double>();
        for (var ts : this.timeseries) {
            for (double v : ts) {
                values.add(v);
            }
        }
        this.distinctTsValues = new double[values.size()];

        int i = 0;
        for (double v : values)
            this.distinctTsValues[i++] = v;

        // compute max coords = last coordinate of each time series
        // maxDimCoords has one entry more, this is the size of the axis for all y
        // values

        this.maxCoordsTS = new int[timeseries.length];

        for (int j = 0; j < timeseries.length; j++) {
            maxCoordsTS[j] = timeseries[j].length - 1;
        }

        this.loopMultiplier = k;

        // length of mean
        this.lLengths = k * (maxDimension - 1) + 1;

        this.sLengths = this.distinctTsValues.length;

        this.table = new MultiDimArray(this.maxCoordsTS, this.lLengths, this.sLengths);

        // Create array for distinct values appearance
        this.appDist = new int[this.distinctTsValues.length][k];

        for (int t = 0; t < k; t++) {
            double[] ts = this.timeseries[t];
            for (int d = 0; d < this.distinctTsValues.length; d++) {
                double distinctValue = this.distinctTsValues[d];
                this.appDist[d][t] = Integer.MAX_VALUE;
                for (int z = 0; z < ts.length; z++) {
                    if (ts[z] == distinctValue) {
                        this.appDist[d][t] = z;
                        break;
                    }
                }
            }
        }

        // fill in the ArrayList of sum Coords
        this.sumCoords = new ArrayList<>(sumLengthTimeseries);

        for (int l = 0; l <= sumLengthTimeseries; l++) {
            this.sumCoords.add(new ArrayList<>());
        }

        int[] originCoordinte = new int[timeseries.length];

        for (int o = 0; o < originCoordinte.length - 1; o++) {
            originCoordinte[o] = 0;
        }

        this.sumCoords.get(addCoordinates(originCoordinte)).add(originCoordinte);

        int[] successor = getNextTSCoordinate(originCoordinte);

        while (successor != null) {
            this.sumCoords.get(addCoordinates(successor)).add(successor);
            successor = getNextTSCoordinate(successor);
        }

    }

    // alternative constructor for diagonal calc

    public MsmMean(double[][] timeseries, double c, int diagParam) {
        this(timeseries, c);
        this.diag = diagParam;
    }

    public Pair<double[], Double> calcMean() {

        calculateTableEntry();

        // Start Coord for backtracking is the Table Entry with max coordinates for all
        // time series and minimal cost
        // l = max mean length

        double cost = Double.POSITIVE_INFINITY;
        int startS = 0;
        int startL = 1;

     //   int initL = (Math.max(1, this.minDimension) - 1);
        int initL = 0;


        for (int s = 0; s < this.distinctTsValues.length; s++) {

            // The mean has to be as least as long as the shortest time series-1

            int flattenMaxCoords = table.flattenCoordinate(maxCoordsTS);
            for (int l = initL; l < lLengths; l++) {
                double tmp = table.get(flattenMaxCoords, l, s);

                if ((tmp < cost) && tmp >= 0.) {
                    cost = tmp;
                    startS = s;
                    startL = l;
                }
            }

        }
        return traceback(maxCoordsTS, startL, startS);
    }

    private Pair<double[], Double> traceback(int[] startCoordsTraceback, int startL, int startS) {

        double[] meanArray = new double[startL + 1];
        int idx = startL;
        int flattenStartCoordsTraceback = table.flattenCoordinate(startCoordsTraceback);
        double cost = this.table.get(flattenStartCoordsTraceback, startL, startS);

        CoordLS predecessor = new CoordLS(startCoordsTraceback, startL, startS);

        meanArray[idx--] = this.distinctTsValues[startS];

        while (predecessor != null) {

            int flattenPredecessorCoord = table.flattenCoordinate(predecessor.coord);
            CoordLS predecessorTmp = getPredecessorT(predecessor.coord, predecessor.l, predecessor.s,
                    flattenPredecessorCoord);

            if (predecessorTmp == null)
                break;

            if (predecessor.l != predecessorTmp.l) {

                meanArray[idx--] = this.distinctTsValues[predecessorTmp.s];
            }

            predecessor = predecessorTmp;
        }

        return new Pair<>(meanArray, cost);

    }

    private CoordLS getPredecessorT(int[] currentCoord, int currentL, int currentS, int flattenCurrentCoord) {

        if (currentL == 0) {
            boolean isOrigin = true;
            for (int c = 0; c < currentCoord.length; c++) {
                if (currentCoord[c] > 0) {
                    isOrigin = false;
                    break;
                }
            }
            if (isOrigin)
                return null;

            int[] predecessor = currentCoord.clone();

            for (int idx = 0; idx < timeseries.length; idx++) {
                if (currentCoord[idx] != 0) {
                    predecessor[idx]--;

                }
            }

            int flattenPredecessorCoord = table.flattenCoordinate(predecessor);

            int[] predecessorME = tracebackME(currentCoord, currentL, currentS, predecessor, flattenPredecessorCoord,
                    flattenCurrentCoord);
            if (predecessorME != null)
                return new CoordLS(predecessorME, currentL, currentS);
        }

        // normal case:
        ArrayList<int[]> predecessorList = Predecessor.getPredecessors(currentCoord, this.diag);

        for (int[] predecessor : predecessorList) {
            int flattenPredecessorCoord = table.flattenCoordinate(predecessor);
            int[] predecessorM = predecessor.clone();
            CoordLS predecessorMOSP = tracebackMOSP(currentCoord, currentL, currentS, predecessorM,
                    flattenPredecessorCoord, flattenCurrentCoord);
            if (predecessorMOSP != null)
                return predecessorMOSP;
            int[] predecessorME = tracebackME(currentCoord, currentL, currentS, predecessor, flattenPredecessorCoord,
                    flattenCurrentCoord);
            if (predecessorME != null)
                return new CoordLS(predecessorME, currentL, currentS);

        }

        return null;
    }

    private CoordLS tracebackMOSP(final int[] currentCoord, int currentL, int currentS, int[] coordsPredecessor,
            int flattenPredecessorCoord, int flattenCurrentCoord) {

        int predecessorL = currentL - 1;

        final double y = distinctTsValues[currentS];
        double sumMove = 0;
        for (int i = 0; i < currentCoord.length; i++) {
            // only sum up if in M Index
            sumMove += (currentCoord[i] - coordsPredecessor[i]) * Math.abs(timeseries[i][currentCoord[i]] - y);
        }

        double currentCost = table.get(flattenCurrentCoord, currentL, currentS);

        for (int s = 0; s < distinctTsValues.length; s++) {

            final double costPred = table.get(flattenPredecessorCoord, predecessorL, s);
            if (costPred < 0)
                continue;

            double sumSplit = 0;
            for (int i = 0; i < currentCoord.length; i++) {
                // only sum up if in split index (the index that doesnt change)
                if ((currentCoord[i] - coordsPredecessor[i]) == 0) {
                    sumSplit += C(y, timeseries[i][currentCoord[i]], distinctTsValues[s]);
                }
            }

            double costMOSP = costPred + sumMove + sumSplit;

            if (costMOSP == currentCost) {
                return new CoordLS(coordsPredecessor, predecessorL, s);
            }
        }

        return null;

    }

    private int[] tracebackME(final int[] currentCoord, int currentL, int currentS, int[] coordsPredecessor,
            int flattenPredecessorCoord, int flattenCurrentCoord) {

        double currentCost = table.get(flattenCurrentCoord, currentL, currentS);

        // evtl falsch
        final double costPred = table.get(flattenPredecessorCoord, currentL, currentS);
        if (costPred < 0.0)
            return null;

        double sumMerge = 0;

        final double y = distinctTsValues[currentS];

        for (int i = 0; i < currentCoord.length; i++) {
            final int currentCoordI = currentCoord[i];
            if ((currentCoordI - coordsPredecessor[i]) == 1) {
                sumMerge += C(timeseries[i][currentCoordI], timeseries[i][currentCoordI - 1], y);
            }
        }

        double costME = costPred + sumMerge;

        if (costME == currentCost) {
            return coordsPredecessor;
        }

        return null;

    }

    public void setOriginCoordinates() {
        // Set cost for origin coordinates
        int[] originCoordinte = new int[timeseries.length];

        int flattenOriginCoordinate = table.flattenCoordinate(originCoordinte);
        double[][] firstDimTableOrigin = table.getFirstDim(flattenOriginCoordinate);
        double[] secondDimTableOrigin = table.getSecondDim(firstDimTableOrigin, 0);

        for (int s = 0; s < this.distinctTsValues.length; s++) {

            double sumMove = 0;
            for (int i = 0; i < k; i++) {
                sumMove += Math.abs(timeseries[i][originCoordinte[i]] - this.distinctTsValues[s]);
            }

            table.set(secondDimTableOrigin, s, sumMove);

        }
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

                for (int l = 0; l < Math.min(lLengths,
                        this.loopMultiplier * (maxCoord) + 1); l++) {

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

                        // mean coordinate = 0. All other timeseries are only able to merge, the move
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

    public double costBorderCoords(final int[] currentCoord, final int l, final int s) {
        int[] predecessor = currentCoord.clone();

        for (int idx = 0; idx < timeseries.length; idx++) {
            if (currentCoord[idx] != 0) {
                predecessor[idx]--;

            }
        }
        int flattenCoordsPredecessor = table.flattenCoordinate(predecessor);
        double[][] firstDimTablePredecessor = table.getFirstDim(flattenCoordsPredecessor);

        return calcMECost(currentCoord, l, s, predecessor, flattenCoordsPredecessor, firstDimTablePredecessor);

    }

    public double calcMOSPCost(final int[] currentCoord, final int l, final int s, int[] coordsPredecessor,
            double[][] firstDimTablePredecessor) {
        // If there exists a Merge operation, all other time series are pausing. So we
        // have to consider only Merge operations OR only Split and Move Operations

        // create coords for predeceasing table entry.
        // Case Move and Split

        int lPred = l - 1;
        double[] secondDimTablePredecessor = table.getSecondDim(firstDimTablePredecessor, lPred);

        double cost = Double.POSITIVE_INFINITY;

        final double y = distinctTsValues[s];
        double sumMove = 0;
        for (int i = 0; i < currentCoord.length; i++) {
            // only sum up if in M Index
            sumMove += (currentCoord[i] - coordsPredecessor[i]) * Math.abs(timeseries[i][currentCoord[i]] - y);
        }

        outer: for (int sPred = 0; sPred < distinctTsValues.length; sPred++) {

            final double costPred = table.getThirdDim(secondDimTablePredecessor, sPred);
            if (costPred < 0)
                continue;

            double sumSplit = 0;
            for (int i = 0; i < currentCoord.length; i++) {
                // only sum up if in split index (the index that doesnt change)
                if ((currentCoord[i] - coordsPredecessor[i]) == 0) {
                    sumSplit += C(y, timeseries[i][currentCoord[i]], distinctTsValues[sPred]);
                }
            }

            final double costMOSP = costPred + sumMove + sumSplit;

            if (costMOSP < cost) {
                cost = costMOSP;
            }
        }

        return cost;
    }

    public double calcMECost(final int[] currentCoord, final int l, final int s, int[] coordsPredecessor,
            int flattenCoordsPredecessor, double[][] firstDimTablePredecessor) {

        double[] seconDimTablePredecessor = table.getSecondDim(firstDimTablePredecessor, l);
        final double costPred = table.getThirdDim(seconDimTablePredecessor, s);
        if (costPred < 0.0)
            return Double.POSITIVE_INFINITY;

        double sumMerge = 0;

        final double y = distinctTsValues[s];

        for (int i = 0; i < currentCoord.length; i++) {
            final int currentCoordI = currentCoord[i];
            if ((currentCoordI - coordsPredecessor[i]) == 1) {
                sumMerge += C(timeseries[i][currentCoordI], timeseries[i][currentCoordI - 1], y);
            }
        }
        return costPred + sumMerge;

    }

    public double C(double new_point, double x, double y) {

        // c - cost of Split/Merge operation. Change this value to what is more
        // appropriate for your data.

        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }

    // Function to get the successing coordinate of a current coordinate
    public int[] getNextTSCoordinate(int[] currentCoord) {
        int[] next = currentCoord.clone();

        boolean carry = true;

        int i = 0;

        while (carry) {
            if (i >= currentCoord.length)
                return null;
            if (currentCoord[i] + 1 <= this.maxCoordsTS[i]) {
                next[i] = next[i] + 1;
                return next;
            } else {
                next[i] = 0;
                i++;
            }
        }

        return next;

    }

    private int addCoordinates(int[] currentCoord) {
        int sum = 0;
        for (int i : currentCoord) {
            sum += i;
        }
        return sum;
    }

    public static void main(String[] args) {

      int n = 10;

        double[][] exp1 = new double[3][n];

        for (int i = 0; i < 3; i++) {

            for (int j = i; j < (n + i); j++) {
                exp1[i][j - i] = j;
            }

            System.out.println(Arrays.toString(exp1[i]));
        }
        System.out.println("Start Computing: k=3 and n=" + n);

        MsmMean msmMean = new MsmMean(exp1, 0.1);

        Pair<double[], Double> mean = msmMean.calcMean();

        System.out.println("Mean");
        System.out.println(Arrays.toString(mean.getFirst()));
        System.out.println(mean.getSecond());

        


    }

}
