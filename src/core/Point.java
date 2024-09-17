package core;

public class Point {
    private double x;
    private double y;

    public Point() {
    }
    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }
    public Point(Point that) {
        this.x = that.x;
        this.y = that.y;
    }
    public double getX() {
        return x;
    }
    public double getY() {
        return y;
    }
}
