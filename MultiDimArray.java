
public class MultiDimArray {

    // private final ArrayList<Double> array;
    double[][][] array;
    final int[] dimensionLengths;
    int sumLengths;
    int[] multiplier;

    public MultiDimArray(int[] dimensionLengths, int lLengths, int sLengths) {

        this.dimensionLengths = dimensionLengths;
        this.multiplier = new int[dimensionLengths.length];
        this.multiplier[0] = 1;
        this.sumLengths = dimensionLengths[0] + 1;

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

    public double get(final int flattenCoords, final int l, final int s) {
        return array[flattenCoords][l][s];
    }

    public double[][] getFirstDim(final int flattenCoords) {
        return this.array[flattenCoords];
    }

    public double[] getSecondDim(final double[][] arrayFirstDim, final int l) {
        return arrayFirstDim[l];
    }

    public double getThirdDim(final double[] arraySecondDim, final int s) {
        return arraySecondDim[s];
    }

    public void set(final double[] arraySecondDim, final int s, final double element) {
        arraySecondDim[s] = element;
    }

    public int flattenCoordinate(final int[] coords) {
        int index = 0;
        for (int i = 0; i < coords.length; i++) {
            index += multiplier[i] * coords[i];
        }
        return index;
    }

}
