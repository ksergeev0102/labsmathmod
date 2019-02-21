import javax.swing.plaf.basic.BasicInternalFrameTitlePane;
import java.util.Arrays;
import java.util.Scanner;

public class Matrix {
    private int size;
    private double[][] data;

    public Matrix(int size) {
        this.size = size;
        data = new double[size][size];
    }

    public int getSize() {
        return size;
    }

    public void setData(int x, int y, double value) {
        data[x][y] = value;
    }

    public double getData(int x, int y) {
        return data[x][y];
    }

    public Matrix additiveDecompositionL() {
        //нижнетреугольная матрица аддитивного разложения для самосопряженной матрицы
        Matrix B = new Matrix(getSize());
        for (int i = 1; i <= B.getSize() - 1; i++) {
            for (int j = 0; j < i; j++) {
                B.setData(i, j, getData(i, j));
            }
        }
        for (int k = 0; k < B.getSize(); k++) {
            B.setData(k, k, getData(k, k) / 2);
        }
        return B;
    }

    public Matrix additiveDecompositionR() {
        //верхнетреугольная матрица аддитивного разложения для самосопряженной матрицы
        Matrix B = new Matrix(getSize());
        for (int i = B.getSize() - 2; i >= 0; i--) {
            for (int j = B.getSize() - 1; j > i; j--) {
                B.setData(i, j, getData(i, j));
            }
        }
        for (int k = 0; k < B.getSize(); k++) {
            B.setData(k, k, getData(k, k) / 2);
        }
        return B;
    }

    public Matrix multByNumber(double value) {
        // умножение матрицы на скаляр
        Matrix A = new Matrix(getSize());
        for (int i = 0; i < getSize(); i++) {
            for (int j = 0; j < getSize(); j++) {
                A.setData(i, j, getData(i, j) * value);
            }
        }
        return A;
    }

    public static Matrix getUnitMatrix(int size) {
        // создать единичную матрицу размера size
        Matrix unitMatrix = new Matrix(size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    unitMatrix.setData(i, j, 1);
                } else {
                    unitMatrix.setData(i, j, 0);
                }
            }
        }
        return unitMatrix;
    }

    public static Matrix multByMatrix(Matrix A, Matrix B) {
        double[][] value = new double[A.getSize()][A.getSize()];
        Matrix C = new Matrix(A.getSize());
        for (int i = 0; i < C.getSize(); ++i) {
            for (int j = 0; j < C.getSize(); ++j) {
                for (int k = 0; k < C.getSize(); ++k) {
                    value[i][j] += A.getData(i, k) * B.getData(k, j);
                    C.setData(i, j, value[i][j]);
                }
            }
        }
        return C;
    }

    public void editMatrix() {
        Scanner scanner = new Scanner(System.in);
        for (int i = 0; i < getSize(); i++) {
            for (int j = 0; j < getSize(); j++) {
                setData(i, j, scanner.nextDouble());
            }
        }
    }

    public void showMatrix() {
        for (int i = 0; i < getSize(); i++) {
            for (int j = 0; j < getSize(); j++) {
                System.out.printf("  %f  ", getData(i, j));
            }
            System.out.println();
        }
    }

    public static Vector multByVector(Matrix A, Vector v) {
        Vector C = new Vector(A.getSize());
        double[] value = new double[A.getSize()];
        for (int i = 0; i < A.getSize(); i++) {
            for (int j = 0; j < A.getSize(); j++) {
                value[i] += v.getData(j) * A.getData(i, j);
            }
            C.setData(i, value[i]);
        }
        return C;
    }

    public static Matrix sumMatrix(Matrix A, Matrix B) {
        Matrix C = new Matrix(A.getSize());
        for (int i = 0; i < C.getSize(); ++i) {
            for (int j = 0; j < C.getSize(); ++j) {
                C.setData(i, j, A.getData(i, j) + B.getData(i, j));
            }
        }
        return C;
    } //сумма матриц

    public static Matrix diffMatrix(Matrix A, Matrix B) {
        Matrix C = new Matrix(A.getSize());
        for (int i = 0; i < C.getSize(); ++i) {
            for (int j = 0; j < C.getSize(); ++j) {
                C.setData(i, j, A.getData(i, j) - B.getData(i, j));
            }
        }
        return C;
    } //разность матриц

    public Matrix forMethod(double w) {
        // создаем матрицу перехода А
        Matrix A = multByMatrix(sumMatrix(getUnitMatrix(getSize()), additiveDecompositionL().multByNumber(w)),
                sumMatrix(getUnitMatrix(getSize()), additiveDecompositionR().multByNumber(w)));
        return A;
    }

    public static Matrix inversion(Matrix A) {
        double temp;
        Matrix E = getUnitMatrix(A.getSize()); //единичная матрица
        for (int k = 0; k < A.getSize(); k++) {
            temp = A.getData(k, k);
            for (int j = 0; j < A.getSize(); j++) {
                A.setData(k, j, A.getData(k, j) / temp);
                E.setData(k, j, E.getData(k, j) / temp);
            }
            for (int i = k + 1; i < A.getSize(); i++) {
                temp = A.getData(i, k);
                for (int j = 0; j < A.getSize(); j++) {
                    A.setData(i, j, A.getData(i, j) - A.getData(k, j) * temp);
                    E.setData(i, j, E.getData(i, j) - E.getData(k, j) * temp);
                }
            }
        }
        for (int k = A.getSize() - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                temp = A.getData(i, k);

                for (int j = 0; j < A.getSize(); j++) {
                    A.setData(i, j, A.getData(i, j) - A.getData(k, j) * temp);
                    E.setData(i, j, E.getData(i, j) - E.getData(k, j) * temp);
                }
            }
        }

        for (int i = 0; i < A.getSize(); i++)
            for (int j = 0; j < A.getSize(); j++)
                A.setData(i, j, E.getData(i, j));
        return A;
    }

    public double[] evaluationOwnNumbers() {
        double numbers1[] = new double[getSize()];
        double numbers2[] = new double[getSize()];
        double own[] = new double[2];
        for (int i = 0; i < getSize(); i++) {
            double temp = 0;
            for (int j = 0; j < getSize(); j++) {
                temp += Math.abs(getData(i, j));
            }
            numbers1[i] = getData(i, i) - temp + Math.abs(getData(i, i));
            numbers2[i] = getData(i, i) + temp - Math.abs(getData(i, i));
        }
        Arrays.sort(numbers1);
        Arrays.sort(numbers2);
        own[0] = min(numbers1[0], numbers2[0]);
        own[1] = max(numbers1[numbers1.length - 1], numbers2[numbers2.length - 1]);
        return own;
    } //оценка собственных значений сверху и снизу с помощью кругов гершгорина


    public double Omega() {
        if (evaluationOwnNumbers()[0] <= 0) {
            return 2 / (Math.sqrt(0.1 * evaluationOwnNumbers()[1]));
        } else return 2 / (Math.sqrt(evaluationOwnNumbers()[0] * evaluationOwnNumbers()[1]));
    } //параметр для ПТМ

    public double Tau() {
        if (evaluationOwnNumbers()[0] <= 0) {
            return 2 / ((0.25 / (1 + Math.sqrt(0.1 / evaluationOwnNumbers()[1])))
                    + Math.sqrt(0.1 * evaluationOwnNumbers()[1]) / 4);
        } else
            return 2 / ((0.5 * evaluationOwnNumbers()[0] / (1 + Math.sqrt(evaluationOwnNumbers()[0] / evaluationOwnNumbers()[1])))
                    + Math.sqrt(evaluationOwnNumbers()[0] * evaluationOwnNumbers()[1]) / 4);
    } //параметр для ПТМ

    public static double max(double a, double b) {
        if (a > b) {
            return a;
        } else return b;
    }

    public static double min(double a, double b) {
        if (a > b) {
            return b;
        } else return a;
    }

    public Matrix HausholderN(int N) {
        Vector x = Vector.matrixToVector(this, N); //вектор под диагональю
        Matrix E = Matrix.getUnitMatrix(getSize());
        Matrix E1 = Matrix.getUnitMatrix(getSize() - N);
        Matrix B = Matrix.diffMatrix(E1, Vector.vectorGenerating(Vector.getUtilVector(x.getSize(), 0), x).multByNumber(2));
        for (int i = getSize() - B.getSize(); i < getSize(); i++) {
            for (int j = getSize() - B.getSize(); j < getSize(); j++) {
                E.setData(i, j, B.getData(i + B.getSize() - getSize(), j + B.getSize() - getSize()));
            }
        }
        return E;
    }
}