package core.system;
import core.Expression;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class SystemOfEquations {
    private boolean simplified = false;
    private final Expression[] equations;
    private Matrix matrix;
    private ArrayList<Character> listOfVariables;
    private HashMap<Character, Double>[] coefficientsMap;

    public SystemOfEquations(Expression[] equations){
        this.equations = new Expression[equations.length];
        System.arraycopy(equations, 0, this.equations, 0, equations.length);
        setMatrix();
    }
    public String toString(){
        Expression[] equationArray = getEquations();
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < equationArray.length; i++) {
            stringBuilder.append(equationArray[i]);
            if (i != equationArray.length - 1){
                stringBuilder.append('\n');
            }
        }
        return stringBuilder.toString();
    }
    public Expression[] getEquations() {
        Expression[] copyEquations = new Expression[equations.length];
        if (!simplified){
            for (int i = 0; i < equations.length; i++) {
                copyEquations[i] = Expression.coefficientsToLinearExpression(coefficientsMap[i]);
            }
            simplified = true;
        }
        else {
            for (int i = 0; i < equations.length; i++) {
                copyEquations[i] = equations[i].clone();
            }
        }
        return copyEquations;

    }
    private HashMap<Character, Double>[] getCoefficientsMap(){
        if (coefficientsMap == null){
            coefficientsMap = coefficients();
        }
        return coefficientsMap;
    }
    public Matrix getMatrix() {
        if (matrix == null){
            setMatrix();
        }
        return matrix;
    }
    public ArrayList<Character> getListOfVariables() {
        if (listOfVariables == null){
            listOfVariables = listOfVariables();
        }
        return listOfVariables;
    }
    private void setMatrix() {
        ArrayList<Character> list = getListOfVariables();
        double[][] coefficients = new double[equations.length][list.size() + 1];
        HashMap<Character, Double>[] maps = getCoefficientsMap();

        for (int i = 0; i < equations.length; i++) {
            for (int j = 0; j < list.size(); j++) {
                coefficients[i][j] = maps[i].getOrDefault(list.get(j), 0.0);
            }
            coefficients[i][list.size()] = maps[i].getOrDefault('\0', 0.0);  // Handling constant term explicitly
        }
        matrix = new Matrix(coefficients);
    }
    private HashMap<Character, Double>[] coefficients(){
        HashMap<Character, Double>[] mapOfCoefficients = new HashMap[equations.length];
        for (int i = 0; i < mapOfCoefficients.length; i++) {
            mapOfCoefficients[i] = equations[i].coefficientsOfLinearExpression();
        }
        return mapOfCoefficients;
    }
    private ArrayList<Character> listOfVariables() {
        HashMap<Character, Double>[] maps = getCoefficientsMap();
        ArrayList<Character> variablesList = new ArrayList<>();
        if (maps.length == 0) {
            return variablesList;
        }
        Set<Character> set = new HashSet<>(maps[0].keySet());
        for (int i = 1; i < maps.length; i++) {
            set.addAll(maps[i].keySet());
        }
        set.remove('\0');
        variablesList.addAll(set);
        return variablesList;
    }
    public String[] getSolutions(){
        return matrix.solve(getListOfVariables());
    }
}
