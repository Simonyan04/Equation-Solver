package core.function;

import core.Expression;
import exceptions.InvalidExpressionExceptions;

public class Linear extends Polynomial{
    public Linear(Expression equation) throws InvalidExpressionExceptions {
        super(equation);
    }
    public double[] solve(){
        double[] coefficients = getCoefficients();
        double[] solutions = new double[1];
        solutions[0] = -(double)coefficients[0]/coefficients[1];
        return solutions;
    }
    public boolean isSolvable() {
        return true;
    }
}
