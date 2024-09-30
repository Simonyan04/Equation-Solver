package core;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Pattern;

import static core.function.Polynomial.HIGHEST_POSSIBLE_DEGREE;

public class Expression implements Cloneable{
    private final String expression;
    public Expression(String expression){
        this.expression = expression;
    }
    public String getExpression() {
        return expression;
    }
    public String toString() {
        return expression;
    }
    public Expression clone(){
        try {
            return (Expression) super.clone();
        }
        catch (CloneNotSupportedException e) {
            return null;
        }
    }
    public boolean equals(Object that){
        if (that == null){
            return false;
        }
        if (that.getClass() != getClass()){
            return false;
        }
        Expression givenExpression = (Expression) that;
        return getExpression().equals(givenExpression.getExpression());
    }
    public boolean isPolynomial() {
        String pattern = "^[+-]?(\\d*x(\\^\\d+)?|\\d+)([+-](\\d*x(\\^\\d+)?|\\d+))*$";
        return Pattern.matches(pattern, expression);
    }
    public HashMap<Character, Double> coefficientsOfLinearExpression(){
        HashMap<Character, Double> map = new HashMap<>();
        boolean otherSide = false;
        if (expression.isEmpty()){
            return map;
        }
        int i = 0;
        while (i<expression.length()){
            char currentChar = expression.charAt(i);
            if (currentChar == '='){
                otherSide = true;
            }
            if (Character.isAlphabetic(currentChar)){
                double value = getCoefficientAtIndex(i);
                if (otherSide){
                    value = -value;
                }
                if (map.containsKey(currentChar)) {
                    double currentValue = map.get(currentChar);
                    map.put(currentChar, currentValue + value);
                } else {
                    map.put(currentChar, value);
                }
            }
            else if (Character.isDigit(currentChar)){
                double[] zeroCoefficientAndEndIndex = zeroDegreeCoefficient(i);
                double value = zeroCoefficientAndEndIndex[0];
                if (otherSide){
                    value = -value;
                }
                if (map.containsKey('\0')) {
                    map.put('\0', map.get('\0') + value);
                } else {
                    map.put('\0', value);
                }
                i = (int)zeroCoefficientAndEndIndex[1];
                continue;
            }
            i++;
        }
        if (!map.containsKey('\0')){
            map.put('\0', 0.0);
            return map;
        }
        if (map.get('\0') != 0){
            map.put('\0', -map.get('\0'));
        }
        return map;
    }
    public double[] coefficientsOfPolynomialExpression(){
        if (expression.isEmpty()){
            return new double[0];
        }
        double[] coefficients = new double[HIGHEST_POSSIBLE_DEGREE];
        int highestDegree = 0;
        int i = 0;
        while (i<expression.length()){
            if(Character.isAlphabetic(expression.charAt(i))){
                int[] degreeAndEndIndex = getDegreeAtIndex(i);
                coefficients[degreeAndEndIndex[0]] += getCoefficientAtIndex(i);
                if(degreeAndEndIndex[0]>highestDegree){
                    highestDegree = degreeAndEndIndex[0];
                }
                i = degreeAndEndIndex[1]+1;
                continue;
            }
            else if(Character.isDigit(expression.charAt(i))){
                double[] zeroCoefficientAndEndIndex = zeroDegreeCoefficient(i);
                coefficients[0] += zeroCoefficientAndEndIndex[0];
                i = (int)zeroCoefficientAndEndIndex[1];
                continue;
            }
            i++;
        }
        double[] newArr = new double[highestDegree+1];
        System.arraycopy(coefficients, 0, newArr, 0, highestDegree+1);
        return newArr;
    }

//    public double[] coefficientsOfPolynomialExpression(){
//        if (expression.isEmpty()){
//            return new double[0];
//        }
//        ArrayList<Double> coefficients = new ArrayList<>(5);
//        int highestDegree = 0;
//        int i = 0;
//        while (i<expression.length()){
//            if(Character.isAlphabetic(expression.charAt(i))){
//                int[] degreeAndEndIndex = getDegreeAtIndex(i);
//                coefficients.set(degreeAndEndIndex[0], coefficients.get(degreeAndEndIndex[0]) + getCoefficientAtIndex(i));
//                i = degreeAndEndIndex[1]+1;
//                continue;
//            }
//            else if(Character.isDigit(expression.charAt(i))){
//                double[] zeroCoefficientAndEndIndex = zeroDegreeCoefficient(i);
//                coefficients.set(0, coefficients.get(0) + zeroCoefficientAndEndIndex[0]);
//                i = (int)zeroCoefficientAndEndIndex[1];
//                continue;
//            }
//            i++;
//        }
//        double[] newArr = new double[coefficients.size()];
//        for (int j = 0; j < newArr.length; j++) {
//            newArr[j] = coefficients.get(j);
//        }
//        //System.arraycopy(coefficients, 0, newArr, 0, highestDegree+1);
//        return newArr;
//    }

    public static Expression coefficientsToLinearExpression(HashMap<Character, Double> coefficients){
        StringBuilder stringBuilder = new StringBuilder();
        Set<Character> set = coefficients.keySet();
        boolean first = true;
        for (Character element : set) {
            if (element == '\0'){
                continue;
            }
            double currentCoefficient = coefficients.get(element);
            if (currentCoefficient > 0 && !first){
                stringBuilder.append("+");
            }
            else if (currentCoefficient < 0 && (first || currentCoefficient == -1)){
                stringBuilder.append("-");
            }
            if (Math.abs(currentCoefficient) != 1){
                stringBuilder.append(currentCoefficient);
            }
            stringBuilder.append(element);
            first = false;
        }
        stringBuilder.append("=");
        stringBuilder.append(coefficients.get('\0'));
        return new Expression(stringBuilder.toString());
    }
    public static Expression coefficientsToPolynomialExpression(double[] coefficients) {
        if (coefficients == null || coefficients.length == 0){
            return new Expression("");
        }
        StringBuilder res = new StringBuilder();
        for (int i = coefficients.length - 1; i >= 0; i--) {
            if (coefficients[i] == 0) {
                continue;
            }
            if (!res.isEmpty()) {
                if (coefficients[i] > 0) res.append("+");
                else if (coefficients[i] < 0) res.append("-");
            } else {
                if (coefficients[i] < 0) res.append("-");
            }
            if (Math.abs(coefficients[i]) != 1 || i == 0) {
                res.append(Math.abs(coefficients[i]));
            }
            if (i >= 1) {
                res.append("x");
                if (i != 1) {
                    res.append("^").append(i);
                }
            }
        }
        return new Expression(res.toString());
    }
    private double getCoefficientAtIndex(int index){
        if (index == 0){
            return 1;
        }
        int start = index-1;
        while (start >= 0 && (Character.isDigit(expression.charAt(start)) || expression.charAt(start) == '.')){
            start--;
        }
        if (start == index-1){
            if (expression.charAt(start) == '-'){
                return -1;
            }
            return 1;
        }
        double result = Double.parseDouble(expression.substring(start+1, index));
        if(start >= 0 && expression.charAt(start) == '-'){
            return -result;
        }
        return result;
    }
    private double[] zeroDegreeCoefficient(int index){
        double[] res = new double[2];
        if((index == 0 || expression.charAt(index-1) != '^')){
            int end = index;
            while (end<expression.length() && (Character.isDigit(expression.charAt(end)) || expression.charAt(end) == '.')){
                end++;
            }
            res[1] = end;
            boolean isPositive = index == 0 || expression.charAt(index-1) == '+' || expression.charAt(index-1) == '=';
            if(end == expression.length() || !Character.isAlphabetic(expression.charAt(end))){
                double result = Double.parseDouble(expression.substring(index, end));
                if (isPositive){
                    res[0] = result;
                }
                else{
                    res[0] = -result;
                }
            }
        }
        return res;
    }
    private int[] getDegreeAtIndex(int index){
        int[] result = new int[2];
        if(index+1 >= expression.length() || expression.charAt(index+1) != '^'){
            result[0] = 1;
            result[1] = index+1;
        }
        else{
            index+= 2;
            int end = index + 1;
            while (end < expression.length() && Character.isDigit(expression.charAt(end))){
                end++;
            }
            result[0] = Integer.parseInt(expression.substring(index, end));
            result[1] = end;
        }
        return result;
    }
}