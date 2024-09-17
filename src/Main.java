import core.Expression;
import core.function.Polynomial;
import core.system.SystemOfEquations;
import java.util.Arrays;

public class Main {
    public static void main(String[] args) {
        Expression expression1 = new Expression("5z+y=1");
        Expression expression2 = new Expression("2x+3y=4y-1");
        Expression[] expressions = new Expression[2];
        expressions[0] = expression1;
        expressions[1] = expression2;
        SystemOfEquations system = new SystemOfEquations(expressions);
        System.out.println(Arrays.toString(system.getSolutions()));
        //Expression expression = new Expression("x^2-x^4");
        Expression expression = new Expression("x^12");
        //System.out.println(expression.coefficientsOfLinearExpression());
        Polynomial polynomial = new Polynomial(expression);
        System.out.println(Arrays.toString(polynomial.solve()));
        System.out.println(polynomial.numberOfRootsOn(-100,100));
    }
}