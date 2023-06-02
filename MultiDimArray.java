public class MultiDimArray {

    // private final ArrayList<Double> array;
    double[][][] array;
    final int[] dimensionLengths;
    int sumLengths;
    int[] multiplier;

    /**
     * Set up a multidimensional array of size: Number of time series + number of different values + size of mean
     * @param dimensionLengths Number of time series
     * @param lLengths size of mean
     * @param sLengths number of different values
     */
    public MultiDimArray(int[] dimensionLengths, int lLengths, int sLengths) {

        this.dimensionLengths = dimensionLengths;
        this.multiplier = new int[dimensionLengths.length];
        this.multiplier[0] = 1;
        this.sumLengths = dimensionLengths[0] + 1;

        // the coordinates of dimensionLengths are flattened
        for (int i = 1; i < dimensionLengths.length; i++) {
            multiplier[i] = sumLengths;
            sumLengths *= (dimensionLengths[i] + 1);
        }

        this.array = new double[sumLengths][lLengths][sLengths];

        for (int i = 0; i < sumLengths; i++) {
            for (int j = 0; j < lLengths; j++) {
                for (int k = 0; k < sLengths; k++) {
                    this.array[i][j][k] = -1.0;
                }
            }

        }
    }

    /**
     * Get table entry
     * @param flattenCoords time series coordinates
     * @param l mean position
     * @param s index of value
     * @return
     */
    public double get(final int flattenCoords, final int l, final int s) {
        return array[flattenCoords][l][s];
    }

    /**
     * get entries for current time series coordinates
     * @param flattenCoords flattened time series coordinates
     * @return two-dimensional array for different mean lengths and values in V(X)
     */
    public double[][] getFirstDim(final int flattenCoords) {
        return this.array[flattenCoords];
    }

    /**
     *  get entries for current time series coordinates and mean position
     * @param arrayFirstDim first dimension array
     * @param l mean position
     * @return one-dimensional array of entries for different values in V(X)
     */
    public double[] getSecondDim(final double[][] arrayFirstDim, final int l) {
        return arrayFirstDim[l];
    }

    /**
     * get table entry
     * @param arraySecondDim second dimension array
     * @param s index for value in V(X)
     * @return table entry
     */
    public double getThirdDim(final double[] arraySecondDim, final int s) {
        return arraySecondDim[s];
    }

    /**
     * set table entry
     * @param arraySecondDim second dimension array
     * @param s index for value in V(X)
     * @param element to set
     */
    public void set(final double[] arraySecondDim, final int s, final double element) {
        arraySecondDim[s] = element;
    }

    /**
     * flatten coordinate
     * @param coords coordinates
     * @return int
     */
    public int flattenCoordinate(final int[] coords) {
        int index = 0;
        for (int i = 0; i < coords.length; i++) {
            index += multiplier[i] * coords[i];
        }
        return index;
    }

}
