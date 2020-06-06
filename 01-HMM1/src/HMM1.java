import java.util.Scanner;


public class HMM1 {
    // Probability of Emission Sequence
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

    public static int[] str2Array(String sr){
        String[] splitSr = sr.split(" ");
        int [] emSeq = new int[splitSr.length-1];
        for (int i = 0; i < splitSr.length - 1; i++){
            emSeq[i] = Integer.parseInt(splitSr[i+1]);
        }
        return emSeq;
    }

    public static double forwardAlgorithm(double[][] A, double[][] B, double[][] pi, int[] O){
        int N = A[0].length; // number of states in the model
        int T = O.length; // length of the observation sequence
        double alphaSum = 0.0;

        double[][] alpha = new double[T][N];

        for (int t = 0; t < T; t++){
            if (t == 0){
                for (int i = 0; i < N; i++){
                    alpha[0][i] = pi[0][i] * B[i][O[0]];
                }
            } else {
                for (int i = 0; i < N; i++){
                    double result = 0;
                    for (int j = 0; j < N; j++){
                        result += alpha[t-1][j] * A[j][i];
                    }
                    /*this is the probability of the partial
                    observation sequence up to time "t" where
                    the underlying Markov process is in state
                    O[t] at time "t".
                     */
                    alpha[t][i] = result * B[i][O[t]];
                }
            }
        }
        for (int i = 0; i < N; i++){
            alphaSum += alpha[T-1][i];
        }
        return alphaSum;
    }

    public static void main( String args[]) {

        Scanner scanString = new Scanner(System.in);

        String lineA = scanString.nextLine();
        String lineB = scanString.nextLine();
        String linePi = scanString.nextLine();
        String lineObs = scanString.nextLine();

        double[][] A = str2Mat(lineA); // state transition probabilities
        double[][] B = str2Mat(lineB); // observation probability matrix
        double[][] pi = str2Mat(linePi); // initial state distribution
        int[] O = str2Array(lineObs); // observation sequence:

        double alpha = forwardAlgorithm(A, B, pi, O); // alpha-pass
        System.out.println(alpha);
    }
}

/* Sample input:
4 4 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.8 0.1 0.1 0.0
4 4 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.1 0.0 0.0 0.9
1 4 1.0 0.0 0.0 0.0
8 0 1 2 3 0 1 2 3

Sample output:
0.090276
 */