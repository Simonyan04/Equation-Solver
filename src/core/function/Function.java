package core.function;

import core.Expression;

public abstract class Function {
    public static final double PRECISION = 1e-16;
    public static final int MAX_ITERATIONS = 1000000;
    private Expression expression;
    public Function(Expression expression){
        setExpression(expression);
    }
    protected void setExpression(Expression expression) {
        this.expression = expression.clone();
    }
    public Expression getExpression() {
        return expression.clone();
    }
    public String toString(){
        return expression.toString();
    }
    public abstract boolean isSolvable();
    protected abstract double[] getSolutions();
    protected abstract Function derivative();
    public abstract double valueAt(double x);
}
