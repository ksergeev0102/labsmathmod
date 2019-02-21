import java.util.Arrays;

import static java.lang.Math.abs;

public class Polinom {
    private double l;
    private double a, b;
    private int n;
    private double x[];
    private double N[][];

    public Polinom(int n, double a, double b) {
        this.n = n;
        this.a = a;
        this.b = b;
        this.l = ((b - a) / n);
        this.x = new double[n + 1];
        for (int i = 0; i < n + 1; i++) {
            x[i] = a + ((b - a) / n) * i;
        }
        this.N = new double[n][n];
    }


    public double f(double x) {
        return Math.cos(x);
    }

    public double N(int i, double t, int p) {
        double result;
        if (p == 0) {
            if (t >= x[i] && t < x[i + 1]) {
                return 1;
            } else {
                return 0;
            }
        }
        result = ((t - x[i]) / (x[i + p] - x[i])) * N(i, t, p - 1) + ((x[i + p + 1] - t) / (x[i + p + 1] - x[i + 1])) * N(i + 1, t, p - 1);
        return result;
    }

    public double Lagr(double x1, double x2, double x3, double x4, double x) {
        return (f(x1) * (x - x2) * (x - x3) * (x - x4))
                / ((x1 - x2) * (x1 - x3) * (x1 - x4))
                + (f(x2) * (x - x1) * (x - x3) * (x - x4))
                / ((x2 - x1) * (x2 - x3) * (x2 - x4))
                + (f(x3) * (x - x2) * (x - x1) * (x - x4))
                / ((x3 - x2) * (x3 - x1) * (x3 - x4))
                + (f(x4) * (x - x2) * (x - x3) * (x - x1))
                / ((x4 - x2) * (x4 - x3) * (x4 - x1));
    }

    public double maximum() {
        double values[] = new double[100];
        double h = (x[3] - x[0]) / 100;
        for (int i = 0; i < 100; i++) {
            values[i] = abs(Lagr(x[0], x[1], x[2], x[3], x[3] - h * i)
                    - Spline3(x[3] - h * i, CoeffSpline3(-Math.cos(a), -Math.cos(b))));
        }
        Arrays.sort(values);
        return values[values.length - 1];
    }

    public double LagrDiff1(double x1, double x2, double x3, double x4, double x) {
        return (f(x1) * (3 * Math.pow(x, 2) - 2 * x3 * x - 2 * x * x2 + x2 * x3 - 2 * x * x4 + x3 * x4 + x2 * x4))
                / ((x1 - x2) * (x1 - x3) * (x1 - x4))
                + (f(x2) * (3 * Math.pow(x, 2) - 2 * x3 * x - 2 * x * x1 + x1 * x3 - 2 * x * x4 + x3 * x4 + x1 * x4))
                / ((x2 - x1) * (x2 - x3) * (x2 - x4))
                + (f(x3) * (3 * Math.pow(x, 2) - 2 * x2 * x - 2 * x * x1 + x1 * x2 - 2 * x * x4 + x2 * x4 + x1 * x4))
                / ((x3 - x2) * (x3 - x1) * (x3 - x4))
                + (f(x4) * (3 * Math.pow(x, 2) - 2 * x2 * x - 2 * x * x1 + x1 * x2 - 2 * x * x3 + x2 * x3 + x1 * x3))
                / ((x4 - x2) * (x4 - x3) * (x4 - x1));
    }

    public double LagrDiff2(double x1, double x2, double x3, double x4, double x) {
        return (f(x1) * (6 * x - 2 * x3 - 2 * x2 - 2 * x4))
                / ((x1 - x2) * (x1 - x3) * (x1 - x4))
                + (f(x2) * (6 * x - 2 * x3 - 2 * x1 - 2 * x4))
                / ((x2 - x1) * (x2 - x3) * (x2 - x4))
                + (f(x3) * (6 * x - 2 * x1 - 2 * x2 - 2 * x4))
                / ((x3 - x2) * (x3 - x1) * (x3 - x4))
                + (f(x4) * (6 * x - 2 * x3 - 2 * x2 - 2 * x1))
                / ((x4 - x2) * (x4 - x3) * (x4 - x1));
    }

    public double Spline3(double data, double gamma[]) {
        int i = 0;
        for (int m = 0; m < n; m++) {
            if (data >= this.x[m] && data <= this.x[m + 1]) {
                i = m;
            }
        }
        return f(x[i]) * (x[i + 1] - data) / l
                + f(x[i + 1]) * (data - x[i]) / l
                + gamma[i] * (Math.pow(x[i + 1] - data, 3) - Math.pow(l, 2) * (x[i + 1] - data)) / (6 * l)
                + gamma[i + 1] * (Math.pow(data - x[i], 3) - Math.pow(l, 2) * (data - x[i])) / (6 * l);
    }

    public double[] CoeffSpline3(double gamma0, double gammaN) {
        Matrix S = new Matrix(n - 1);
        for (int i = 0; i < n - 1; i++) {
            S.setData(i, i, 4 * l);
        }
        for (int i = 0; i < n - 2; i++) {
            S.setData(i, i + 1, l);
        }
        for (int i = 1; i < n - 1; i++) {
            S.setData(i, i - 1, l);
        }
        Vector func = new Vector(n - 1);
        Vector v = new Vector(n - 1);
        Vector solve0 = new Vector(n - 1);
        Methods methods = new Methods(S, v, func);
        for (int i = 0; i < n - 1; i++) {
            func.setData(i, (f(x[i + 2]) - 2 * f(x[i + 1]) + f(x[i])) / l);
        }
        methods.Running();
        solve0 = methods.getSolve();
        Vector solve = new Vector(n + 1);
        solve.setData(0, gamma0);
        solve.setData(n, gammaN);
        for (int i = 1; i < n; i++) {
            solve.setData(i, solve0.getData(i - 1));
        }
        return solve.getData();
    }
}
