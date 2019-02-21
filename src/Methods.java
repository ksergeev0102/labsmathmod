import java.util.ArrayList;

public class Methods {
    private Matrix matrix_of_system;
    private Vector initial_approximation;
    private Vector function;
    private double t, w;
    private int iteration;
    private Vector solve;
    private Matrix H;
    private Matrix R;

    public Methods(Matrix A, Vector v, Vector f) {
        this.matrix_of_system = A;
        this.initial_approximation = v;
        this.function = f;
    }

    public Vector getSolve(){
        return solve;
    }

    public int getIteration() {
        return iteration;
    }

    public double getT() {
        return t;
    }

    public double getW() {
        return w;
    }

    public Matrix getMatrix_of_system() {
        return matrix_of_system;
    }

    public Vector getFunction() {
        return function;
    }

    public Matrix getH() {
        return H;
    }

    public Matrix getR() {
        return R;
    }

    public Vector getInitial_approximation() {
        return initial_approximation;
    }

    public void PTM() {
        this.t = matrix_of_system.Tau();
        this.w = matrix_of_system.Omega();
        Matrix B = matrix_of_system.forMethod(w);
        Matrix B_inversion = Matrix.inversion(B);
        ArrayList<Vector> solves = new ArrayList<>();
        solves.add(initial_approximation);
        solves.add(Vector.sum(solves.get(solves.size() - 1), Matrix.multByVector(B_inversion.multByNumber(t), function),
                Matrix.multByVector(B_inversion.multByNumber(t), Matrix.multByVector(matrix_of_system, solves.get(solves.size() - 1)))));
        int i = 0;
        while (Vector.vectorDifference(Matrix.multByVector(matrix_of_system, solves.get(solves.size() - 1)), function).module() > 0.0000001) {
            solves.add(Vector.sum(solves.get(solves.size() - 1), Matrix.multByVector(B_inversion.multByNumber(t), function),
                    Matrix.multByVector(B_inversion.multByNumber(t), Matrix.multByVector(matrix_of_system, solves.get(solves.size() - 1)))));
            i++;
        }
        this.solve = solves.get(solves.size()-1);
        this.iteration = i;
    }

    public void Richardson() {
        this.t = 2 / (matrix_of_system.evaluationOwnNumbers()[0]+matrix_of_system.evaluationOwnNumbers()[1]);
        Matrix T = Matrix.diffMatrix(Matrix.getUnitMatrix(matrix_of_system.getSize()), matrix_of_system.multByNumber(t));
        Vector V = function.multByNumber(t);
        ArrayList<Vector> solves = new ArrayList<>();
        solves.add(initial_approximation);
        solves.add(Vector.sum(Matrix.multByVector(T, solves.get(solves.size()-1)), V));
        int i = 0;
        while (Vector.vectorDifference(solves.get(solves.size()-1),solves.get(solves.size()-2)).module() > 0.0000001) {
            solves.add(Vector.sum(Matrix.multByVector(T, solves.get(solves.size() - 1)), V));
            i++;
        }
        this.solve = solves.get(solves.size()-1);
        this.iteration = i;
    }

    public void Householder() {
        Matrix R = Matrix.multByMatrix(matrix_of_system.HausholderN(0), matrix_of_system);
        Matrix H = matrix_of_system.HausholderN(0);
        int k = 1;
        for (int i = 1; i < R.getSize() - 1; i++) {
            H = Matrix.multByMatrix(H, R.HausholderN(i));
            R = Matrix.multByMatrix(R.HausholderN(i), R);
            k++;
        }
        System.out.println(k);
        this.solve = Matrix.multByVector(Matrix.inversion(R), Matrix.multByVector(Matrix.inversion(H), function));
        this.H = Matrix.inversion(H);
        this.R = R;
    }

    public void Running() {
        double ksi[] = new double[matrix_of_system.getSize() + 1];
        double n[] = new double[matrix_of_system.getSize() + 1];
        double c[] = new double[matrix_of_system.getSize()];
        ksi[0] = 0;
        n[0] = 0;
        ksi[ksi.length - 1] = 0;
        c[c.length - 1] = 0;
        for (int i = 0; i < c.length - 1; i++) {
            c[i] = matrix_of_system.getData(i, i + 1);
        }
        double f[] = new double[matrix_of_system.getSize()];
        for (int i = 0; i < c.length; i++) {
            f[i] = function.getData(i);
        }
        double a[] = new double[matrix_of_system.getSize()];
        a[0] = 0;
        for (int i = 1; i < a.length; i++) {
            a[i] = matrix_of_system.getData(i, i - 1);
        }
        double b[] = new double[matrix_of_system.getSize()];
        for (int i = 0; i < b.length; i++) {
            b[i] = matrix_of_system.getData(i, i);
        }
        for (int i = 0; i < ksi.length - 1; i++) {
            ksi[i + 1] = -c[i] / (a[i] * ksi[i] + b[i]);
        }
        for (int i = 0; i < n.length - 1; i++) {
            n[i + 1] = (f[i] - a[i] * n[i]) / (a[i] * ksi[i] + b[i]);
        }
        double x[] = new double[matrix_of_system.getSize() + 1];
        x[x.length - 1] = 0;
        for (int i = x.length - 2; i >= 0; i--) {
            x[i] = ksi[i + 1] * x[i + 1] + n[i + 1];
        }
        Vector solve = new Vector(matrix_of_system.getSize());
        for (int i = 0; i < solve.getSize(); i++) {
            solve.setData(i, x[i]);
        }
        this.solve =  solve;
    }
}
