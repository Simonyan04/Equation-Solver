package exceptions;

public class InvalidExpressionExceptions extends Exception{
    public InvalidExpressionExceptions(){
        super("Invalid Expression");
    }
    public InvalidExpressionExceptions(String message){
        super(message);
    }
}
