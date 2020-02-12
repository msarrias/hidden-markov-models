import java.util.Scanner;


public class HMM1 {

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

    public static int[] str2list(String sr){
        String[] splitSr = sr.split(" ");
        int [] emSeq = new int[splitSr.length-1];
        for (int i = 0; i < splitSr.length - 1; i++){
            emSeq[i] = Integer.parseInt(splitSr[i+1]);
        }
        return emSeq;
    }

    public static double alphaPass(double[][] A, double[][] B, double[][] pi, int[] O){
        int hiddenStates = A[0].length; //N of hidden states
        int nEmissions = O.length;
        double alphaSum = 0.0;

        double[][] alpha = new double[nEmissions][hiddenStates];

        for (int t = 0; t < nEmissions; t++){
            if (t == 0){
                for (int i = 0; i < hiddenStates; i++){
                    alpha[0][i] = pi[0][i] * B[i][O[0]];
                }
            } else {
                for (int i = 0; i < hiddenStates; i++){
                    double result = 0;
                    for (int j = 0; j < hiddenStates; j++){
                        result += alpha[t-1][j] * A[j][i];
                    }
                    alpha[t][i] = result * B[i][O[t]];
                }
            }
        }
        for (int i = 0; i < hiddenStates; i++){
            alphaSum += alpha[nEmissions-1][i];
        }
        return alphaSum;
    }

    public static void main( String args[]) {

        Scanner scanString = new Scanner(System.in);

        String lineA = scanString.nextLine();
        String lineB = scanString.nextLine();
        String linePi = scanString.nextLine();
        String lineObs = scanString.nextLine();

        double[][] A = str2Mat(lineA);
        double[][] B = str2Mat(lineB);
        double[][] pi = str2Mat(linePi);
        int[] O = str2list(lineObs);

        double alpha = alphaPass(A, B, pi, O);
        System.out.println(alpha);
    }

}
