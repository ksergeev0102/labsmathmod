import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class Main {
    public static void main(String[] args)  throws IOException {
        Matrix a = new Matrix(3);
        Scanner sc = new Scanner(new File("/Users/macbook/IdeaProjects/labsmathmod/src/data"));
        for (int i = 0; i < a.getSize(); i++)
        {
            for (int j = 0; j < a.getSize(); j++)
            {
                a.setData(i,j,sc.nextDouble());
            }
        }
        Vector v = new Vector(3);
        for(int i = 0;i<v.getSize();i++){
            v.setData(i,sc.nextDouble());
        }
        Vector w = new Vector(3);
        for(int i = 0;i<v.getSize();i++) {
            w.setData(i, sc.nextDouble());
        }
        Methods method = new Methods(a,v,w);
//        method.Householder();
//        method.Running();
//        method.PTM();
//        method.Richardson();
//        FileWriter fileWriter = new FileWriter("results.txt");
//        fileWriter.write("Матрица системы:\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            for(int j = 0;j<method.getInitial_approximation().getSize();j++){
//                fileWriter.write(Double.toString(method.getMatrix_of_system().getData(i,j)));
//                fileWriter.write("  ");
//            }
//            fileWriter.write("\n");
//        }
//        fileWriter.write("\nНачальное приближение\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            fileWriter.write(Double.toString(method.getInitial_approximation().getData(i)));
//            fileWriter.write("  ");
//        }
//        fileWriter.write("\nПравая часть\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            fileWriter.write(Double.toString(method.getFunction().getData(i)));
//            fileWriter.write("  ");
//
//        }
//        fileWriter.write("\nРешение СЛАУ\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            fileWriter.write(Double.toString(method.getSolve().getData(i)));
//            fileWriter.write("  ");
//
//        }
//        fileWriter.write("\nКоличество итераций\n");
//        fileWriter.write(Integer.toString(method.getIteration()));
//        fileWriter.write("\nМатрица Хаусхолдера:\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            for(int j = 0;j<method.getInitial_approximation().getSize();j++){
//                fileWriter.write(Double.toString(method.getH().getData(i,j)));
//                fileWriter.write("  ");
//            }
//            fileWriter.write("\n");
//        }
//        fileWriter.write("\nВерхнетреугольная матрица:\n");
//        for(int i = 0;i<method.getInitial_approximation().getSize();i++){
//            for(int j = 0;j<method.getInitial_approximation().getSize();j++){
//                fileWriter.write(Double.toString(method.getR().getData(i,j)));
//                fileWriter.write("  ");
//            }
//            fileWriter.write("\n");
//        }
//        fileWriter.close();
        Polinom polinom = new Polinom(10,0,10);
        System.out.println(polinom.maximum());
        System.out.println(polinom.H(0,3,4));
    }
}
