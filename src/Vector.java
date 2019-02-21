import java.util.Scanner;

public class Vector {
    private int size;
    private double[] data;

    public Vector(int size){
        this.size = size;
        data = new double[size];
    }

    public double[] getData() {
        return data;
    }

    public int getSize(){
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public double getData(int x) {
        return data[x];
    }

    public void setData(int x,double data) {
        this.data[x] = data;
    }

    public Vector multByNumber(double value){
        for(int i =0;i<size;i++){
            setData(i,getData(i)*value);
        }
        return this;
    }

    public void editVector(){
        Scanner scanner = new Scanner(System.in);
        for (int i = 0; i < getSize(); i++) {
            setData(i,scanner.nextDouble());
        }
    }

    public void showVector(){
        for(int i = 0;i<getSize();i++){
            System.out.printf("  %f  ",getData(i));
        }
    }

    public static Vector vectorDifference(Vector v,Vector w){
        Vector u  = new Vector(v.getSize());
        for(int i = 0;i<v.getSize();i++){
            u.setData(i,v.getData(i)- w.getData(i));
        }
        return u;
    }

    public double module(){
        double mod = 0;
        for(int i = 0;i<getSize();i++){
            mod+=Math.pow(getData(i),2);
        }
        return Math.sqrt(mod);
    }

    public static Vector sum(Vector v,Vector w,Vector u){
        Vector z = new Vector(v.getSize());
        for(int i =0;i<v.getSize();i++){
            z.setData(i,v.getData(i)+w.getData(i)-u.getData(i));
        }
        return z;
    }

    public  static  Vector sum(Vector v,Vector w){
        Vector z = new Vector(v.getSize());
        for(int i =0;i<v.getSize();i++){
            z.setData(i,v.getData(i)+w.getData(i));
        }
        return z;
    }

    public Vector Normalization(){
        double x = module();
        for(int i = 0;i<getSize();i++){
            setData(i,getData(i)/x);
        }
        return this;
    }

    public static Vector matrixToVector(Matrix a,int q){
        Vector v = new Vector(a.getSize()-q);
        for(int i = 0;i<v.getSize();i++){
            v.setData(i,a.getData(a.getSize()-v.getSize()+i,a.getSize()-v.getSize()));
        }
        return v;
    }

    public static Vector getUtilVector(int size,int one){
        Vector e = new Vector(size);
        for(int i = 0;i<size;i++){
            e.setData(i,0);
        }
        e.setData(one,1);
        return e;
    }

    public static Matrix vectorGenerating(Vector e,Vector u){ //u вектор под диагональю который надо обнулить
        Vector w = Vector.vectorDifference(u,e.multByNumber(u.module())).Normalization(); // w вектор хаусхолдера
        Matrix uuT = new Matrix(u.getSize());
        for(int i = 0;i<u.getSize();i++){
            for(int j = 0;j<u.getSize();j++){
                uuT.setData(i,j,w.getData(i)*w.getData(j));
            }
        }
        return uuT;
    }
}
