import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

public class MsmMean {

    // MultiDimArray for storing all coordinates of a table with their respective
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

    // Parameter for the window heuristic, if set to -1, no window
    int window = -1;

    // Cost for Merge and Split Operation
    double c;

    /**
     * MSM Mean computation
     *
     * @param timeseries input time series
     * @param c          constant cost for merge and split
     */
    public MsmMean(double[][] timeseries, double c) {
        this.c = c;
        this.timeseries = timeseries;
        this.k = timeseries.length;
        for (double[] ts : this.timeseries) {
            this.sumLengthTimeseries += ts.length - 1;
            this.maxDimension = Math.max(maxDimension, ts.length);
            this.minDimension = Math.min(minDimension, ts.length);
        }

        // create list of all time series data points (V(X))
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
        this.maxCoordsTS = new int[timeseries.length];

        for (int j = 0; j < timeseries.length; j++) {
            maxCoordsTS[j] = timeseries[j].length - 1;
        }

        this.loopMultiplier = k;

        //max length of mean
        this.lLengths = k * (maxDimension - 1) + 1;

        // number of distinct values in timeseries
        this.sLengths = this.distinctTsValues.length;

        this.table = new MultiDimArray(this.maxCoordsTS, this.lLengths, this.sLengths);

        // Create array for distinct values appearance
        this.appDist = new int[this.distinctTsValues.length][k];
        // Filling array -> first appearances = position of a value of V(X) for each time series
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

    /**
     * alternative constructor for window heuristic
     *
     * @param timeseries  input time series
     * @param c           constant cost for merge and split
     * @param windowParam window size
     */
    public MsmMean(double[][] timeseries, double c, int windowParam) {
        this(timeseries, c);
        this.window = windowParam;
    }

    /**
     * Compute mean
     *
     * @return mean and cost
     */
    public Pair<double[], Double> calcMean() {

        // fill multidimensional table
        calculateTableEntry();

        // Start Traceback:
        // Start entry is the entry where all time series reached their max coordinate
        // and where the cost is minimum for all possible l and s
        double cost = Double.POSITIVE_INFINITY;
        int startS = 0;
        int startL = 1;
        int initL = 0;


        for (int s = 0; s < this.distinctTsValues.length; s++) {

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


    /**
     * set entry for table origin
     */
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

    /**
     * Fill multidimensional table
     */
    public void calculateTableEntry() {

        setOriginCoordinates();

        // iterate through all sums of coordinates
        for (int i = 1; i <= this.sumLengthTimeseries; i++) {

            ArrayList<int[]> sumICoords = this.sumCoords.get(i);
            // Iterate through all time series such that the sum of the current coordinates of the time series is i
            for (int[] currentCoord : sumICoords) {

                int flattenCurrentCoord = table.flattenCoordinate(currentCoord);
                double[][] firstDimTableCurrentCoord = table.getFirstDim(flattenCurrentCoord);

                // get max index of coordinates to limit the intermediate length of the mean during computation
                int maxCoord = 0;
                for (int c : currentCoord) {
                    if (c > maxCoord) maxCoord = c;
                }

                ArrayList<int[]> predecessorList = Predecessor.getPredecessors(currentCoord, this.window);

                for (int l = 0; l < Math.min(lLengths, this.loopMultiplier * (maxCoord) + 1); l++) {

                    double[] secondDimTableCurrentCoord = table.getSecondDim(firstDimTableCurrentCoord, l);

                    outer:
                    for (int s = 0; s < this.distinctTsValues.length; s++) {
                        // Check if the value is valid to be set (if it did not occur in any of the time series until the current position, it is not set
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

                        // if mean coordinate = 0. All other timeseries are only able to merge, the move
                        // and split case is not applied.
                        if (l == 0) {
                            double cost = costBorderCoords(currentCoord, l, s);
                            if (cost == Double.POSITIVE_INFINITY) continue outer;
                            table.set(secondDimTableCurrentCoord, s, cost);
                            continue outer;
                        }

                        // regular case: find minimum cost of move/split or merge
                        // If there exists a merge operation, all other time series are pausing. So we
                        // have to consider only merge operations OR only split and Move Operations
                        double cost = Double.POSITIVE_INFINITY;

                        for (int[] predecessor : predecessorList) {
                            int flattenCoordsPredecessor = table.flattenCoordinate(predecessor);
                            double[][] firstDimTablePredecessor = table.getFirstDim(flattenCoordsPredecessor);
                            int[] predecessorMOSP = predecessor.clone();
                            double costMOSP = calcMOSPCost(currentCoord, l, s, predecessorMOSP, firstDimTablePredecessor);
                            int[] predecessorME = predecessor.clone();
                            double costME = calcMECost(currentCoord, l, s, predecessorME, firstDimTablePredecessor);

                            cost = Math.min(Math.min(cost, costME), costMOSP);
                        }

                        if (cost == Double.POSITIVE_INFINITY) continue outer;
                        table.set(secondDimTableCurrentCoord, s, cost);

                    }

                }
            }
        }

    }

    /**
     * if the mean coordinate is zero, only merges are possible
     *
     * @param currentCoord current coordinates of time series
     * @param l            current position of mean
     * @param s            current value position
     * @return cost
     */
    public double costBorderCoords(final int[] currentCoord, final int l, final int s) {
        int[] predecessor = currentCoord.clone();

        for (int idx = 0; idx < timeseries.length; idx++) {
            if (currentCoord[idx] != 0) {
                predecessor[idx]--;

            }
        }
        int flattenCoordsPredecessor = table.flattenCoordinate(predecessor);
        double[][] firstDimTablePredecessor = table.getFirstDim(flattenCoordsPredecessor);

        return calcMECost(currentCoord, l, s, predecessor, firstDimTablePredecessor);

    }

    /**
     * Compute cost for the case move and split
     *
     * @param currentCoord             current coordinates of time series
     * @param l                        current position of mean
     * @param s                        current value position
     * @param coordsPredecessor        coordinates of all possible predecessors
     * @param firstDimTablePredecessor table entries of all predeceasing entries
     * @return cost
     */
    public double calcMOSPCost(final int[] currentCoord, final int l, final int s, int[] coordsPredecessor, double[][] firstDimTablePredecessor) {

        // for the case move/split the predeceasing entry for mean is decremented
        int lPred = l - 1;

        // get all entries for current time series and mean positions
        double[] secondDimTablePredecessor = table.getSecondDim(firstDimTablePredecessor, lPred);

        double cost = Double.POSITIVE_INFINITY;

        final double y = distinctTsValues[s];
        double sumMove = 0;
        for (int i = 0; i < currentCoord.length; i++) {
            // only sum up if in M Index, that is, the position had changed in comparison to the predeceasing entry
            sumMove += (currentCoord[i] - coordsPredecessor[i]) * Math.abs(timeseries[i][currentCoord[i]] - y);
        }

        // iterate through all possible distinct values for the predeceasing entry in split case
        outer:
        for (int sPred = 0; sPred < distinctTsValues.length; sPred++) {

            final double costPred = table.getThirdDim(secondDimTablePredecessor, sPred);
            if (costPred < 0) continue;

            double sumSplit = 0;
            for (int i = 0; i < currentCoord.length; i++) {
                // only sum up if in split index (the index that does not change)
                if ((currentCoord[i] - coordsPredecessor[i]) == 0) {
                    sumSplit += C(y, timeseries[i][currentCoord[i]], distinctTsValues[sPred]);
                }
            }

            final double costMOSP = costPred + sumMove + sumSplit;

            //find minimum cost
            if (costMOSP < cost) {
                cost = costMOSP;
            }
        }

        return cost;
    }

    /**
     * Compute cost for the case merge
     *
     * @param currentCoord             current coordinates of time series
     * @param l                        current position of mean
     * @param s                        current value position
     * @param coordsPredecessor        coordinates of all possible predecessors
     * @param firstDimTablePredecessor table entries of all predeceasing entries
     * @return cost
     */
    public double calcMECost(final int[] currentCoord, final int l, final int s, int[] coordsPredecessor, double[][] firstDimTablePredecessor) {

        double[] secondDimTablePredecessor = table.getSecondDim(firstDimTablePredecessor, l);
        final double costPred = table.getThirdDim(secondDimTablePredecessor, s);
        if (costPred < 0.0) return Double.POSITIVE_INFINITY;

        double sumMerge = 0;

        final double y = distinctTsValues[s];

        for (int i = 0; i < currentCoord.length; i++) {
            final int currentCoordI = currentCoord[i];
            // only add the cost for those entries that changed
            if ((currentCoordI - coordsPredecessor[i]) == 1) {
                sumMerge += C(timeseries[i][currentCoordI], timeseries[i][currentCoordI - 1], y);
            }
        }
        return costPred + sumMerge;

    }


    /**
     * traceback through the table to find the corresponding mean with lowest cost
     * @param startCoordsTraceback start coordinates time series = max coordinates time series
     * @param startL start coordinate mean
     * @param startS start value
     * @return mean with corresponding cost
     */
    private Pair<double[], Double> traceback(int[] startCoordsTraceback, int startL, int startS) {

        double[] meanArray = new double[startL + 1];
        int idx = startL;
        int flattenStartCoordsTraceback = table.flattenCoordinate(startCoordsTraceback);
        double cost = this.table.get(flattenStartCoordsTraceback, startL, startS);

        CoordLS predecessor = new CoordLS(startCoordsTraceback, startL, startS);

        meanArray[idx--] = this.distinctTsValues[startS];

        while (predecessor != null) {

            int flattenPredecessorCoord = table.flattenCoordinate(predecessor.coord);
            CoordLS predecessorTmp = getPredecessorT(predecessor.coord, predecessor.l, predecessor.s, flattenPredecessorCoord);

            if (predecessorTmp == null) break;

            if (predecessor.l != predecessorTmp.l) {

                meanArray[idx--] = this.distinctTsValues[predecessorTmp.s];
            }

            predecessor = predecessorTmp;
        }

        return new Pair<>(meanArray, cost);

    }

    /**
     * compute predeceasing entry in a computed table - the entry that leads to the minimum cost in the succeeding entry
     * @param currentCoord current coordinate of time series
     * @param currentL current coordinate of mean
     * @param currentS current value position
     * @param flattenCurrentCoord flattened current coordinates
     * @return coordinates for: time series, mean and value
     */
    private CoordLS getPredecessorT(int[] currentCoord, int currentL, int currentS, int flattenCurrentCoord) {

        // border case: only merge
        if (currentL == 0) {
            boolean isOrigin = true;
            for (int c = 0; c < currentCoord.length; c++) {
                if (currentCoord[c] > 0) {
                    isOrigin = false;
                    break;
                }
            }
            if (isOrigin) return null;

            int[] predecessor = currentCoord.clone();

            for (int idx = 0; idx < timeseries.length; idx++) {
                if (currentCoord[idx] != 0) {
                    predecessor[idx]--;

                }
            }

            int flattenPredecessorCoord = table.flattenCoordinate(predecessor);

            int[] predecessorME = tracebackME(currentCoord, currentL, currentS, predecessor, flattenPredecessorCoord, flattenCurrentCoord);
            if (predecessorME != null) return new CoordLS(predecessorME, currentL, currentS);
        }

        // regular case:
        ArrayList<int[]> predecessorList = Predecessor.getPredecessors(currentCoord, this.window);

        for (int[] predecessor : predecessorList) {
            int flattenPredecessorCoord = table.flattenCoordinate(predecessor);
            int[] predecessorM = predecessor.clone();
            CoordLS predecessorMOSP = tracebackMOSP(currentCoord, currentL, currentS, predecessorM, flattenPredecessorCoord, flattenCurrentCoord);
            if (predecessorMOSP != null) return predecessorMOSP;
            int[] predecessorME = tracebackME(currentCoord, currentL, currentS, predecessor, flattenPredecessorCoord, flattenCurrentCoord);
            if (predecessorME != null) return new CoordLS(predecessorME, currentL, currentS);

        }

        return null;
    }

    /**
     * traceback for move and split case
     * @param currentCoord current coordinate of time series
     * @param currentL current coordinate of mean
     * @param currentS current value position
     * @param coordsPredecessor predeceasing coordinates
     * @param flattenPredecessorCoord flattened predeceasing coordinates
     * @param flattenCurrentCoord flattened current coordinates
     * @return coordinates for: time series, mean and value
     */
    private CoordLS tracebackMOSP(final int[] currentCoord, int currentL, int currentS, int[] coordsPredecessor, int flattenPredecessorCoord, int flattenCurrentCoord) {

        int predecessorL = currentL - 1;

        final double y = distinctTsValues[currentS];
        double sumMove = 0;
        for (int i = 0; i < currentCoord.length; i++) {
            // only sum up if in M Index, that is, the position had changed in comparison to the predeceasing entry
            sumMove += (currentCoord[i] - coordsPredecessor[i]) * Math.abs(timeseries[i][currentCoord[i]] - y);
        }

        double currentCost = table.get(flattenCurrentCoord, currentL, currentS);

        for (int s = 0; s < distinctTsValues.length; s++) {

            final double costPred = table.get(flattenPredecessorCoord, predecessorL, s);
            if (costPred < 0) continue;

            double sumSplit = 0;
            for (int i = 0; i < currentCoord.length; i++) {
                // only sum up if in split index (the index that does not change)
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
    /**
     * traceback for merge case
     * @param currentCoord current coordinate of time series
     * @param currentL current coordinate of mean
     * @param currentS current value position
     * @param coordsPredecessor predeceasing coordinates
     * @param flattenPredecessorCoord flattened predeceasing coordinates
     * @param flattenCurrentCoord flattened current coordinates
     * @return coordinates for: time series, mean and value
     */
    private int[] tracebackME(final int[] currentCoord, int currentL, int currentS, int[] coordsPredecessor, int flattenPredecessorCoord, int flattenCurrentCoord) {

        double currentCost = table.get(flattenCurrentCoord, currentL, currentS);

        final double costPred = table.get(flattenPredecessorCoord, currentL, currentS);
        if (costPred < 0.0) return null;

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

    /**
     * cost for merge and split
     * @param new_point
     * @param x
     * @param y
     * @return cost
     */
    public double C(double new_point, double x, double y) {

        // c - cost of Split/Merge operation. Change this value to what is more
        // appropriate for your data.

        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }


    /**
     * get the succeeding coordinate of the current coordinates
     *
     * @param currentCoord current coordinates of the time series
     * @return succeeding coordinates
     */
    public int[] getNextTSCoordinate(int[] currentCoord) {
        int[] next = currentCoord.clone();

        boolean carry = true;

        int i = 0;

        while (carry) {
            // no increase possible, reached end of each time series
            if (i >= currentCoord.length) return null;
            // Increase the next posible coordinate of a time series and return new coordinates
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

    /**
     * Sum up the current coordinates
     *
     * @param currentCoord current coordinates of the time series
     * @return sum
     */
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
