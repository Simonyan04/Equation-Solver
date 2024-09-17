package core.system;

import core.Expression;

import java.util.ArrayList;
import java.util.Arrays;

public class Matrix implements Cloneable{
    private double[][] coefficients;
    private Matrix rref;
    public Matrix(double[][] matrix)
    {
        coefficients = matrix.clone();
    }
    public double getCoefficient(int row, int column){
        return coefficients[row][column];
    }
    public double[] getRow(int row){
        return coefficients[row].clone();
    }
    public void setCoefficient(int row, int column, double value){
        coefficients[row][column] = value;
    }
    public int getNumberOfRows(){
        return coefficients.length;
    }
    public int getNumberOfColumns(){
        return coefficients[0].length;
    }
    public Matrix getRREF(){
        if (rref == null){
            rref = reduceToRREF();
        }
        return rref.clone();
    }
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < getNumberOfRows(); i++) {
            stringBuilder.append(Arrays.toString(coefficients[i]));
            if (i == getNumberOfRows() - 1){
                break;
            }
            stringBuilder.append("\n");
        }
        return stringBuilder.toString();
    }
    public Matrix clone(){
        try {
            Matrix clone = (Matrix) super.clone();
            clone.coefficients = this.coefficients;
            return clone;
        }
        catch (CloneNotSupportedException e) {
            return null;
        }
    }
    public void swap(int i, int j){
        if (i == j){
            return;
        }
        double[] temp = getRow(i);
        coefficients[i] = getRow(j);
        coefficients[j] = temp;
    }
    public void scale(int row, double scalar){
        for (int i = 0; i < getNumberOfColumns(); i++) {
            setCoefficient(row, i, getCoefficient(row, i) * scalar);
        }
    }
    private Matrix reduceToRREF() {
        Matrix matrix = this.clone();
        int lead = 0;
        int rowCount = matrix.getNumberOfRows();
        int colCount = matrix.getNumberOfColumns();
        for (int r = 0; r < rowCount; r++) {
            if (colCount <= lead) {
                break;
            }
            int i = r;
            while (matrix.getCoefficient(i, lead) == 0) {
                i++;
                if (i == rowCount) {
                    i = r;
                    lead++;
                    if (colCount == lead) {
                        lead--;
                        break;
                    }
                }
            }
            matrix.swap(i, r);
            if (matrix.getCoefficient(r,lead) != 0) {
                matrix.scale(r, 1.0 / matrix.coefficients[r][lead]);
            }
            for (i = 0; i < rowCount; i++) {
                if (i != r) {
                    double leadFactor = matrix.coefficients[i][lead];
                    for (int j = 0; j < colCount; j++) {
                        matrix.coefficients[i][j] -= leadFactor * matrix.coefficients[r][j];
                    }
                }
            }
            lead++;
        }
        return matrix;
    }
    public boolean isInconsistent(){
        Matrix rref = getRREF();
        int rowCount = rref.getNumberOfRows();
        int colCount = rref.getNumberOfColumns()-1;
        for (int i = 0; i < rowCount; i++) {
            boolean isAllZero = true;
            for (int j = 0; j < colCount; j++) {
                if (rref.getCoefficient(i, j) != 0) {
                    isAllZero = false;
                    break;
                }
            }
            if (isAllZero && rref.getCoefficient(i, colCount) != 0) {
                return true;
            }
        }
        return false;
    }
    private boolean[] freeVariableIndexes(){
        Matrix rref = getRREF();
        int rowCount = rref.getNumberOfRows();
        int colCount = rref.getNumberOfColumns() - 1;
        boolean[] isFreeVariable = new boolean[colCount];
        Arrays.fill(isFreeVariable, true);
        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < colCount; j++) {
                if (rref.getCoefficient(i, j) == 1) {
                    isFreeVariable[j] = false;
                    break;
                }
            }
        }
        return isFreeVariable;
    }
    public boolean isFreeVariable(int column){
        return freeVariableIndexes()[column];
    }
    public boolean isPivotVariable(int column){
        return !isFreeVariable(column);
    }
    public String[] solve(ArrayList<Character> listOfVariables) {
        String[] solutions = new String[listOfVariables.size()];
        if(isInconsistent()){
            return solutions;
        }
        Matrix rref = getRREF();
        for (int j = 0; j < solutions.length; j++) {
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.append(listOfVariables.get(j));
            if (isFreeVariable(j)) {
                stringBuilder.append(" is a free variable");
                solutions[j] = stringBuilder.toString();
                stringBuilder.setLength(0);
            }
            else {
                stringBuilder.append(" = ");
                for (int i = 0; i < rref.getNumberOfRows(); i++) {
                    if(rref.getCoefficient(i, j) == 1){
                        if (rref.getCoefficient(i, solutions.length) != 0){
                            stringBuilder.append(rref.getCoefficient(i, solutions.length));
                        }
                        for (int k = j+1; k < solutions.length; k++) {
                            double currentCoefficient = -rref.getCoefficient(i, k);
                            if (currentCoefficient != 0){
                                stringBuilder.append(-rref.getCoefficient(i, k));
                                stringBuilder.append(listOfVariables.get(k));
                            }
                        }
                        if (stringBuilder.toString().equals(listOfVariables.get(i) + " = ")){
                            stringBuilder.append(0.0);
                        }
                        solutions[j] = stringBuilder.toString();
                        stringBuilder.setLength(0);
                        break;
                    }
                }
            }
        }
        return solutions;
    }
}
