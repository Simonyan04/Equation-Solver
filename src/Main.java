import core.Expression;
import core.function.Polynomial;
import core.system.SystemOfEquations;
import java.util.Arrays;

public class Main {
    public static void main(String[] args) {
        Expression expression1 = new Expression("5z+y=1");
        Expression expression2 = new Expression("2x+3y=4y-1");
        Expression expression = new Expression("x^6-x^5+x");
        //Expression expression = new Expression("x^12-x");
        Polynomial polynomial = new Polynomial(expression);
//        System.out.println(Arrays.toString(polynomial.getCoefficients()));
//        System.out.println(polynomial.getDerivative());
        System.out.println(Arrays.toString(polynomial.solve()));

    }
}