import java.util.Scanner;

public class HMM0 {
    // Next Emission Distribution
    public static double[][] str2Mat(String sr){
        String[] splitSr = sr.split(" ");
        int nRows = Integer.parseInt(splitSr[0]);
        int nCols = Integer.parseInt(splitSr[1]);

        double[][] matrix = new double[nRows][nCols];

        int l = 0;
        for(int i = 0; i < nRows ; i++){
            for (int j = 0; j < nCols;j++){
                matrix[i][j] =  Double.parseDouble(splitSr[l+2]);
                l++;
            }
        }
        return matrix;
    }

    public static String mat2str(double[][] c){
        String[] d = new String[c.length * c[0].length + 2];
        d[0] = String.valueOf(c.length);
        d[1] = String.valueOf(c[0].length);
        int k =0;
        for (int i = 0; i < c.length; i++){
            for (int j = 0; j < c[0].length; j++){
                d[k+2] = String.valueOf((c[i][j]));
                k++;
            }
        }
        return String.join(" ", d);
    }

    public static double[][] matMul(double[][] a, double[][] b){
        //creating another matrix to store the multiplication of two matrices
        double [][] c = new double[a.length][b[0].length];

        if (a[0].length == b.length){
            for (int i = 0; i < a.length; i++){
                for (int j = 0; j < b[0].length; j++){
                    for (int k = 0; k < a[0].length; k++){
                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        } else {
            System.out.println("Incorrect dimensions for matrix multiplication");
            System.exit(0);
        }
        return c;
    }

    public static void main( String args[]) {
        Scanner scanString = new Scanner(System.in);

        String lineA = scanString.nextLine();
        String lineB = scanString.nextLine();
        String linePi = scanString.nextLine();

        double[][] A = str2Mat(lineA);
        double[][] B = str2Mat(lineB);
        double[][] pi = str2Mat(linePi);

        double[][] stateProb = matMul(pi, A);
        double[][] emissionProb = matMul(stateProb, B);

        System.out.println(mat2str(emissionProb));
    }
}

/* Sample input
4 4 0.2 0.5 0.3 0.0 0.1 0.4 0.4 0.1 0.2 0.0 0.4 0.4 0.2 0.3 0.0 0.5
4 3 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.2 0.6 0.2
1 4 0.0 0.0 0.0 1.0

Sample output
1 3 0.3 0.6 0.1
 */