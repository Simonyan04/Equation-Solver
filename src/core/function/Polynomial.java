package core.function;

import core.Expression;
import core.Point;

import java.util.ArrayList;
import java.util.List;

public class Polynomial extends Function implements Cloneable{
    private static final Polynomial MINUS_ONE = generateFromCoefficients(new double[] {-1.0});
    public static final int HIGHEST_POSSIBLE_DEGREE = 100000;
    private double[] coefficients;
    private Polynomial derivative;
    private final ArrayList<Double> rootIntervals = new ArrayList<>();
    private int numberOfRoots = -1;
    private double[] solutions;

    public Polynomial(Expression expression){
        super(expression);
        coefficients = expression.coefficientsOfPolynomialExpression();
    }
    public static Polynomial generateFromCoefficients(double[] coefficients){
        Polynomial polynomial = new Polynomial(Expression.coefficientsToPolynomialExpression(coefficients));
        polynomial.coefficients = coefficients.clone();
        return polynomial;
    }
    public double[] getCoefficients() {
        double[] array = new double[coefficients.length];
        System.arraycopy(coefficients, 0, array, 0, coefficients.length);
        return array;
    }
    public Polynomial getDerivative(){
        if (derivative == null){
            derivative = derivative();
        }
        return derivative.clone();
    }
    public double[] getSolutions() {
        if (solutions == null){
            solutions = this.solve();
        }
        return solutions;
    }
    private ArrayList<Double> getRootIntervals(){
        if (numberOfRoots == -1){
            ArrayList<Polynomial> sturmSequence = sturmSequenceGenerator();
            double bound = getBound(sturmSequence);
            rootIsolator(-bound, bound+1, sturmSequence);
        }
        return rootIntervals;
    }
    public int getNumberOfRoots(){
        if (numberOfRoots == -1){
            numberOfRoots = numberOfRoots();
        }
        return numberOfRoots;
    }
    public int getHighestDegree(){
        return coefficients.length - 1;
    }
    public Polynomial clone(){
        try {
            Polynomial clone = (Polynomial) super.clone();
            clone.coefficients = getCoefficients().clone();
            if (derivative != null){
                clone.derivative = derivative.clone();
            }
            return clone;
        }
        catch (CloneNotSupportedException e) {
            return null;
        }
    }
    public double valueAt(double x){
        double result = coefficients[0];
        for(int i = 1; i < coefficients.length; i++) {
            result += coefficients[i] * Math.pow(x, i);
        }
        return result;
    }
    protected Polynomial derivative(){
        double[] coefficientsOfDerivative;
        if (coefficients.length == 0){
            coefficientsOfDerivative = new double[0];
            return generateFromCoefficients(coefficientsOfDerivative);
        }
        coefficientsOfDerivative = new double[coefficients.length - 1];
        for (int i = 1; i < coefficients.length; i++) {
            coefficientsOfDerivative[i-1] = coefficients[i] * i;
        }
        return generateFromCoefficients(coefficientsOfDerivative);
    }
    public static Polynomial multiplication(Polynomial firstPolynomial, Polynomial secondPolynomial){
        return generateFromCoefficients(multiplication(firstPolynomial.getCoefficients(),secondPolynomial.getCoefficients()));
    }
    private static double[] multiplication(double[] first, double[] second){
        double[] result = new double[first.length + second.length - 1];
        for (int i = 0; i < first.length; i++) {
            for (int j = 0; j < second.length; j++) {
                result[i + j] += first[i] * second[j];
            }
        }
        return result;
    }
    public static Polynomial scale(double coefficient, Polynomial polynomial) {
        double[] coefficients = new double[polynomial.coefficients.length];
        double[] polynomialCoefficients = polynomial.coefficients;
        for (int i = 0; i < polynomialCoefficients.length; i++) {
            coefficients[i] = coefficient * polynomialCoefficients[i];
        }
        return generateFromCoefficients(coefficients);
    }
    private double newtonMethod(double x0) {
        double x1 = x0 - valueAt(x0) / (getDerivative().valueAt(x0) == 0 ? getDerivative().valueAt(x0-0.0001) : getDerivative().valueAt(x0));
        int iteration = 0;
        while (Math.abs(x1 - x0) > PRECISION && iteration < MAX_ITERATIONS) {
            x0 = x1;
            x1 = x0 - valueAt(x0) / getDerivative().valueAt(x0);
            iteration++;
        }
        return x1;
    }
    private double bisectionMethod(double x1, double x2) {
        double x = x1;
        int iteration = 0;
        while (Math.abs(x2 - x1) > PRECISION && iteration++ < MAX_ITERATIONS) {
            x = (x1 + x2) / 2;
            double fx1 = valueAt(x1);
            double fx = valueAt(x);
            if (fx == 0) return x;
            x1 = fx1 * fx >= 0 ? x: x1;
            x2 = fx1 * fx < 0 ? x: x2;
        }
        return x;
    }
    private static double[] subtractArray(double[] minuend, double[] subtrahend) {
        double epsilon = 1e-15;
        for (int i = 0; i < subtrahend.length; i++) {
            minuend[i] -= subtrahend[i];
        }
        int newLength = 0;
        for (int i = 0; i < minuend.length; i++) {
            if (Math.abs(minuend[i]) > epsilon) {
                newLength = i + 1;
            }
        }
        if (newLength == 0) {
            return new double[]{0.0};
        }
        double[] result = new double[newLength];
        System.arraycopy(minuend, 0, result, 0, newLength);
        return result;
    }
    public static Polynomial remainderOfEuclideanDivisionOfPolynomials(Polynomial dividend, Polynomial divisor){
        int d = divisor.getHighestDegree();
        double[] remainder = dividend.getCoefficients();
        double[] divisorCoefficients = divisor.getCoefficients();
        while (remainder.length > 1 && remainder.length - 1 >= d){
            double[] temp = new double[remainder.length-divisorCoefficients.length+1];
            temp[temp.length-1] = remainder[remainder.length-1]/divisorCoefficients[divisorCoefficients.length-1];
            remainder = subtractArray(remainder, multiplication(divisorCoefficients,temp));
        }
        return generateFromCoefficients(remainder);
    }

    private ArrayList<Polynomial> sturmSequenceGenerator(){
        ArrayList<Polynomial> sturmSequence = new ArrayList<>();
        Polynomial firstTerm = this;
        Polynomial current = firstTerm.getDerivative();
        sturmSequence.add(firstTerm);
        sturmSequence.add(current);
        int i = 1;
        while (current.getHighestDegree()>0){
            current = multiplication(remainderOfEuclideanDivisionOfPolynomials(sturmSequence.get(i-1), sturmSequence.get(i++)), MINUS_ONE);
            sturmSequence.add(current);
        }
        return sturmSequence;
    }
    private int numberOfRootsOn(double start, double end, ArrayList<Polynomial> sturmSequence){
        boolean sign1 = sturmSequence.getFirst().valueAt(start) > 0;
        boolean sign2 = sturmSequence.getFirst().valueAt(end) > 0;
        int counter1 = 0;
        int counter2 = 0;
        for (int i = 1; i < sturmSequence.size(); i++) {
            if ((sturmSequence.get(i).valueAt(start) > 0 && !sign1) || (sturmSequence.get(i).valueAt(start) < 0 && sign1)){
                counter1++;
                sign1 = !sign1;
            }
            if ((sturmSequence.get(i).valueAt(end) > 0 && !sign2) || (sturmSequence.get(i).valueAt(end) < 0 && sign2)){
                counter2++;
                sign2 = !sign2;
            }
        }
        return counter1-counter2;
    }
    public int numberOfRoots(){
        return getRootIntervals().size()-1;
    }
    private int numberOfRoots(ArrayList<Polynomial> sturmSequence){
        if(numberOfRoots != -1){
            return numberOfRoots;
        }
        double[] current = sturmSequence.getFirst().getCoefficients();
        boolean sign1;
        boolean sign2;
        if (current.length % 2 == 1){
            sign1 = current[current.length-1]>0;
            sign2 = sign1;
        }
        else{
            sign1 = current[current.length-1]<0;
            sign2 = !sign1;
        }
        int counter1 = 0;
        int counter2 = 0;
        for (int i = 1; i < sturmSequence.size(); i++) {
            current = sturmSequence.get(i).getCoefficients();
            if (current.length % 2 == 1){
                if ((current[current.length-1]>0 && !sign1) || (current[current.length-1]<0 && sign1)){
                    counter1++;
                    sign1 = !sign1;
                }
            }
            else{
                if ((current[current.length-1]<0 && !sign1) || (current[current.length-1]>0 && sign1)){
                    counter1++;
                    sign1 = !sign1;
                }
            }
            if ((current[current.length-1]>0 && !sign2) || (current[current.length-1]<0 && sign2)){
                counter2++;
                sign2 = !sign2;
            }
        }
        numberOfRoots = counter1-counter2;
        return numberOfRoots;
    }
    private void rootIsolator(double first, double second, ArrayList<Polynomial> sturmSequence){
        if (rootIntervals.size() == numberOfRoots(sturmSequence)){
            rootIntervals.add(first);
            return;
        }
        double lowerBound = first;
        double upperBound = second;
        while (numberOfRootsOn(lowerBound, upperBound, sturmSequence) > 1 || (upperBound - lowerBound) > 10){
            double mid = (lowerBound + upperBound) / 2;
            if (numberOfRootsOn(lowerBound, mid, sturmSequence) != 0) upperBound = mid;
            else lowerBound = mid;
        }
        rootIntervals.add(lowerBound);
        rootIsolator(upperBound, second, sturmSequence);
    }
    private double getBound(ArrayList<Polynomial> sturmSequence){
        int n = numberOfRoots(sturmSequence);
        int bound = 10;
        while (numberOfRootsOn(-bound, bound, sturmSequence) != n) bound *= 2;
        return bound;
    }

    private double[] solve(){
        int n = getNumberOfRoots();
        double[] solutions = new double[n];
        ArrayList<Double> intervals = getRootIntervals();
        for(int i = 0; i < n; i++) {
            double x1 = intervals.get(i);
            double x2 = intervals.get(i+1);
            solutions[i] = valueAt(x1) * valueAt(x2) < 0 ? bisectionMethod(x1, x2) : newtonMethod((x1 + x2)/2);
        }
        return solutions;
    }
    public boolean isSolvable(){
        return getNumberOfRoots() != 0;
    }
    public static Polynomial interpolate(List<Point> points) {
        int n = points.size();
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = points.get(i).getX();
            y[i] = points.get(i).getY();
        }
        double[][] dividedDifferences = new double[n][n];
        System.arraycopy(y, 0, dividedDifferences[0], 0, n);
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < n-i ; j++) {
                dividedDifferences[i][j] = (dividedDifferences[i-1][j+1] - dividedDifferences[i-1][j]) / (x[i+j] - x[j]);
            }
        }
        double[][] coeffsArray = new double[n][n];
        coeffsArray[0][0]=1;
        for (int i=1;i<n;i++){
            double[] newTerm = new double[2];
            newTerm[0]=-x[i-1];
            newTerm[1]=1;
            coeffsArray[i]=multiplication(coeffsArray[i-1],newTerm);
        }
        for (int i=0;i<n;i++) {
            for (int j = 0; j <= i; j++) {
                coeffsArray[i][j] = dividedDifferences[i][0] * coeffsArray[i][j];
            }
        }
        double[] columnSums = new double[n];
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                sum += coeffsArray[i][j];
            }
            columnSums[j] = sum;
        }
        return generateFromCoefficients(columnSums);
    }
    public static double interpolateValue(double xValue, List<Point> points) {
        Polynomial interpolatingPolynomial = interpolate(points);
        return interpolatingPolynomial.valueAt(xValue);
    }
}