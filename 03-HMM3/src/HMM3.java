import java.util.Scanner;
import java.util.ArrayList;
import java.util.List;

public class HMM3 {

    public static Double[] scale;

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

    public static int[] str2Array(String sr){
        String[] splitSr = sr.split(" ");
        int [] emSeq = new int[splitSr.length-1];
        for (int i = 0; i < splitSr.length - 1; i++){
            emSeq[i] = Integer.parseInt(splitSr[i+1]);
        }
        return emSeq;
    }

    public static double[][] forwardAlgorithm(double[][] A, double[][] B, double[][] pi, int[] O){
        int N = A.length; // number of states in the model
        int T = O.length; // length of the observation sequence

        scale = new Double[T];

        //alpha measures the relevant probability up to time t
        // we will scale the numbers to avoid underflow problems.
        double[][] alpha = new double[T][N];

        for (int t = 0; t < T; t++){
            scale[t] = 0.0;
            if (t == 0){
                for (int j = 0; j < N; j++) {
                    alpha[0][j] = pi[0][j] * B[j][O[0]];
                    scale[0] += alpha[0][j];
                }
                scale[0] = 1.0 / scale[0];
                for (int i = 0; i < N; i++) {
                    alpha[0][i] = alpha[0][i] * scale[0];
                }
            } else {
                // predicted distribution
                for (int i = 0; i < N; i++){
                    double result = 0;
                    for (int j = 0; j < N; j++){
                        result += alpha[t-1][j] * A[j][i];
                    }
                    // improved distribution
                    alpha[t][i] = result * B[i][O[t]];
                    scale[t] += alpha[t][i];
                }
                scale[t] = 1.0 / scale[t];
                for (int i = 0; i < N; i++){
                    alpha[t][i] = alpha[t][i] * scale[t];
                }
            }
        }
        return alpha;
    }

    public static double[][] backwardAlgorithm(double[][] A, double[][] B, int[] O){
        int N = A.length; // number of states in the model
        int T = O.length; // length of the observation sequence

        // beta measures the relevant probability after time t.
        double[][] beta = new double[T][N];

        for (int t = T-1 ; t >= 0; t--){
            if (t == T-1){
                for (int i = 0; i < N; i++){
                    beta[T-1][i] = scale[T-1];
                }
            } else {
                for (int i = 0; i < N; i++){
                    double result = 0.0;
                    for (int j = 0; j < N; j++){
                        result += beta[t+1][j] * A[i][j] * B[j][O[t+1]];
                    }
                    beta[t][i] = result * scale[t];
                }
            }
        }
        return beta;
    }

    public static List<double[][]> BaumWelchAlgorithm(double[][] A, double[][] B, double[][] pi, int[] O) {
        // Reference: https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
        int N = A.length; // number of states
        int M = B[0].length; // number of observation outcomes
        int T = O.length; // number of emissions
        int maxIterations = 100; // max number of re-estimation iterations
        int iteration = 0; // current iteration
        double oldLogProb = Integer.MIN_VALUE;
        double[][] gamma = new double[T][N];
        double[][][] diGamma = new double[T][N][N];
        List<double[][]> results = new ArrayList<double[][]>();

        while (true) {
            // 1. initialize lambda = (A,B,pi) the first guess is given
            // so we skip this step.
            // 2.1 compute alpha (Forward algorithm/ alpha-pass)
            double[][] alpha = forwardAlgorithm(A, B, pi, O);
            // 2.2 compute beta (Backward algorithm / beta-pass)
            double[][] beta = backwardAlgorithm(A, B, O);
            // 2.3 compute gama(i) and diGamma(i,j)
            for (int t = 0; t < T; t++) {
                double evidence = 0.0;
                if (t == T - 1) {
                    for (int i = 0; i < N; i++) {
                        evidence += alpha[T - 1][i];
                    }
                    for (int i = 0; i < N; i++) {
                        gamma[T - 1][i] = alpha[T - 1][i] / evidence;
                    }
                } else {
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            evidence += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
                        }
                    }
                    for (int i = 0; i < N; i++) {
                        gamma[t][i] = 0.0;
                        for (int j = 0; j < N; j++) {
                            diGamma[t][i][j] = (alpha[t][i] * A[i][j] * beta[t + 1][j] * B[j][O[t + 1]]) / evidence;
                            gamma[t][i] += diGamma[t][i][j];
                        }
                    }
                }
            }
            // 3. re-estimate the model lambda = (A,B,pi)
            // 3.1 update pi
            for (int i = 0; i < N; i++) {
                pi[0][i] = gamma[0][i];
            }
            // 3.2 update A
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    double tempGamma = 0.0;
                    double tempDiGamma = 0.0;
                    for (int t = 0; t < T - 1; t++) {
                    /* tempDiGamma gives the expected number of
                     transitions from state ki to state kj.
                     tempGamma gives the expected number of transitions
                     from state ki to any state.
                     */
                        tempDiGamma += diGamma[t][i][j];
                        tempGamma += gamma[t][i];
                    }
                /* this ratio gives the probability of transitioning
                from state ki to state kj: aij.
                 */
                    A[i][j] = tempDiGamma / tempGamma;
                }
            }
            // 3.2 update B
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    double numerator = 0.0;
                    double denominator = 0.0;
                    for (int t = 0; t < T; t++) {
                        if (O[t] == j) {
                            numerator += gamma[t][i];
                        }
                        denominator += gamma[t][i];
                    }
                    B[i][j] = numerator / denominator;
                }
            }
            // 4. Repeat 2 until convergence
            double logProb = 0.0;
            for (int t = 0; t < T; t++) {
                logProb += Math.log(scale[t]);
            }
            logProb = -logProb;

            iteration++;
            if (iteration < maxIterations && logProb > oldLogProb) {
                oldLogProb = logProb;
            } else {
                break;
            }
        }
        results.add(A);
        results.add(B);
        return results;
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
        int[] O = str2Array(lineObs);

        List<double[][]> estimateModelParam = BaumWelchAlgorithm(A, B, pi, O);

        System.out.println(mat2str(estimateModelParam.get(0))); // A
        System.out.println(mat2str(estimateModelParam.get(1))); // B
    }
}


/* Sample input 1:
4 4 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4
4 4 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4
1 4 0.241896 0.266086 0.249153 0.242864
1000 0 1 2 3 3 0 0 1 1 1 2 2 2 3 0 0 0 1 1 1 2 3 3 0 0 0 1 1 1 2 3 3 0 1 2 3 0 1 1 1 2 3 3 0 1 2 2 3 0 0 0 1 1 2 2 3 0 1 1 2 3 0 1 2 2 2 2 3 0 0 1 2 3 0 1 1 2 3 3 3 0 0 1 1 1 1 2 2 3 3 3 0 1 2 3 3 3 3 0 1 1 2 2 3 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 3 3 3 3 0 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 0 1 2 3 0 1 1 1 2 3 0 1 1 2 2 2 2 2 3 0 1 1 1 2 2 2 2 3 0 0 0 0 0 1 1 1 1 2 2 3 3 0 1 2 3 3 0 0 0 0 0 0 1 1 2 2 3 0 0 1 1 1 1 1 1 2 3 3 0 0 1 1 1 2 3 0 0 1 2 3 0 1 1 2 3 3 0 0 0 1 2 3 3 3 0 1 1 1 1 2 3 3 3 3 3 3 0 1 2 2 2 2 2 2 3 0 1 1 1 2 2 3 3 3 3 0 1 2 3 0 0 0 1 1 2 2 3 0 0 0 0 0 0 0 1 2 2 2 3 3 3 3 0 0 1 2 2 2 3 3 3 0 0 1 2 2 3 0 0 0 0 1 1 1 2 3 3 3 3 3 3 3 3 0 1 2 3 0 0 1 2 3 3 3 0 0 0 0 0 1 1 1 1 2 3 0 0 0 1 2 2 3 3 0 0 0 1 1 1 1 1 2 3 3 3 3 0 1 1 1 2 2 3 0 1 2 3 3 3 3 0 0 0 0 1 2 3 3 0 1 2 2 3 3 0 0 1 1 2 3 3 0 1 2 2 3 3 3 0 0 1 1 2 3 3 3 3 0 0 1 1 2 3 3 0 1 2 3 0 1 1 2 2 3 0 1 2 3 3 0 1 1 1 2 2 2 3 3 0 0 1 1 1 1 1 2 3 3 3 0 1 1 2 2 2 2 3 3 0 0 1 2 3 0 1 1 2 2 2 2 3 0 0 1 2 2 3 0 0 0 0 0 1 1 1 2 3 0 0 1 2 3 3 0 0 0 1 2 2 2 3 3 0 0 0 1 2 2 2 2 2 3 0 1 1 2 3 0 0 1 1 1 2 2 3 0 0 0 0 1 1 1 2 2 3 0 1 1 1 2 2 2 3 3 0 0 1 2 2 3 3 3 0 1 1 2 3 0 0 0 0 0 1 2 2 2 3 3 3 0 0 0 1 2 3 0 1 1 2 3 3 3 0 1 2 2 2 3 0 0 1 1 1 1 2 3 3 0 0 0 0 1 2 3 3 3 0 0 0 1 1 2 3 0 1 1 1 1 2 2 2 2 2 2 3 0 0 0 0 1 2 2 2 2 3 0 1 2 2 3 0 1 2 3 0 1 2 3 0 0 0 1 1 2 2 3 3 0 1 1 1 1 2 2 3 3 0 1 1 1 2 2 2 3 3 3 0 1 1 2 3 3 0 1 2 3 0 0 0 0 1 2 3 0 0 0 0 0 0 1 2 2 3 3 0 0 1 2 3 0 1 2 2 3 0 0 0 1 1 2 2 2 2 2 3 3 3 3 3 0 1 2 2 3 3 3 3 3 0 0 1 1 2 2 3 0 0 1 2 2 3 3 3 0 0 0 1 2 2 2 2 3 3 0 1 2 3 0 0 1 1 1 2 2 3 0 0 1 1 2 2 2 3 3 0 0 1 1 1 1 1 2 3 3 3 0 1 2 2 2 2 3 3 3 3 3 3 0 0 0 0 0 0 1 2 3 0 0 1 1 1 2 3 0 0 1 1 2 2 2 2 3 3 3 0 1 1 2 2 2 3 3 0 0 0 0 0 0 1 2 2 3 3 0 0 0 0 0 0 1 2 3 3 3 0 1 1 1 2 2 2 2 2 3 3 3 0 1 2 2 2 3 3 3 3 0 0 0 0 1 2 3 3 3 3 3 3 0 0 1 1 1 1 2 3 0 1 2 3 0 1 1 2 3 3 3 0 0 0 0 1 1 2 3 3 3 3 0 0 1 1 1 2 2 2 2 2 2 3 3 0 0 0 1 2 3 0 0 1 1 2 2 3 3 3 3 3 0 0 1 2 2 2 2 3 0 0 1 1 1 1 1 2 3 3 0 0 1 1 1 2 3 3 3 0 0

Sample output 1:
4 4 0.545455 0.454545 0.0 0.0 0.0 0.506173 0.493827 0.0 0.0 0.0 0.504132 0.495868 0.478088 0.0 0.0 0.521912
4 4 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0

Sample input 2:
4 4 0.2 0.4 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.4 0.2
4 4 0.7 0.1 0.1 0.1 0.1 0.7 0.1 0.1 0.1 0.1 0.7 0.1 0.1 0.1 0.1 0.7
1 4 1.0 0.0 0.0 0.0
1000 1 0 1 0 2 3 0 1 0 1 2 3 2 0 1 0 3 2 3 2 3 2 1 0 1 0 1 2 3 2 0 2 1 0 1 3 2 3 2 3 2 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 2 3 2 0 1 0 1 2 3 2 3 0 3 2 1 0 3 2 1 2 3 0 1 0 1 2 0 3 0 1 0 3 0 1 2 3 1 2 3 2 1 0 1 0 3 2 3 2 1 2 1 0 2 1 3 2 3 2 3 1 0 1 0 3 2 3 2 0 1 0 2 1 0 1 0 3 1 0 1 0 2 1 0 3 2 3 2 3 2 3 0 1 2 1 0 1 0 1 0 3 2 3 2 0 1 3 2 3 2 3 0 1 0 3 2 3 2 3 0 1 3 2 3 1 0 1 0 3 1 0 1 2 3 0 1 0 3 2 3 0 1 0 1 2 1 0 3 0 1 0 1 3 2 3 1 0 1 2 1 0 3 2 1 0 2 3 1 2 3 2 3 2 3 2 3 0 3 0 1 2 3 2 1 3 2 0 3 2 3 2 3 2 1 0 1 3 0 1 3 2 3 2 0 3 2 3 2 3 2 3 0 1 0 1 2 3 2 3 2 1 0 1 0 3 0 1 0 2 3 0 1 0 1 0 1 2 1 0 1 0 1 3 2 3 2 3 2 3 0 2 1 2 3 2 3 2 0 1 0 1 0 1 0 3 2 3 2 3 0 1 0 1 0 2 3 2 3 2 3 2 3 2 1 0 2 3 2 3 0 1 2 3 2 1 0 1 2 0 3 1 0 1 0 3 0 2 3 1 0 1 2 1 2 3 0 3 0 1 0 3 2 1 0 1 0 3 2 0 1 0 1 2 3 0 3 2 1 0 3 2 0 1 0 1 0 1 0 1 0 1 2 0 1 3 2 3 2 3 0 2 1 0 3 2 1 0 1 0 1 0 2 3 2 1 0 1 0 1 0 1 2 3 0 1 2 0 1 0 1 0 1 0 1 0 3 2 1 0 1 0 1 0 1 2 3 2 3 0 1 0 1 0 3 2 3 2 1 2 1 0 1 0 1 2 3 2 3 2 3 2 3 0 3 0 1 0 3 2 3 2 3 2 1 2 3 2 3 2 3 0 1 3 0 2 1 0 1 0 1 3 0 1 0 1 0 1 0 3 2 3 2 1 0 1 3 1 0 1 0 3 0 1 0 2 3 2 3 2 3 2 3 2 1 2 3 2 3 2 1 0 1 0 1 2 3 2 3 2 1 0 1 2 1 0 1 3 0 3 2 3 2 3 0 1 0 1 2 3 2 1 0 1 0 1 3 2 1 0 3 2 1 0 3 2 3 2 0 1 0 1 2 3 0 2 3 2 3 2 1 0 1 2 3 2 1 0 1 3 2 0 1 2 3 0 2 3 2 3 0 1 0 1 2 3 2 0 1 0 3 2 0 1 0 1 0 1 0 1 0 1 0 1 3 2 3 2 1 2 3 2 1 0 1 0 1 0 1 0 3 0 1 0 1 0 3 0 1 0 1 2 3 0 1 0 3 0 1 0 1 0 1 2 3 0 1 0 3 0 1 0 1 2 3 0 1 0 1 2 1 3 0 1 0 1 0 1 0 3 2 3 2 3 2 1 0 1 0 1 3 0 3 0 1 2 3 1 3 2 3 2 1 3 1 2 1 0 1 0 1 0 2 1 0 1 0 1 2 3 2 1 0 2 3 2 3 2 3 2 3 0 3 2 3 0 1 2 3 2 3 2 3 2 3 1 0 1 0 1 0 1 2 0 3 0 3 2 0 3 2 1 0 1 0 1 0 1 2 3 2 3 2 1 2 3 0 1 0 3 2 0 3 2 0 1 0 1 0 3 0 1 0 1 2 3 2 3 2 1 0 1 0 1 2 3 2 3 2 3 2 1 0 1 0 1 0 3 0 1 2 1 0 3 2 3 2 3 2 0 3 2 3 2 3 0 1 2 3 1 0 1 0 3 2 3 2 1 2 3 0 3 2 3 2 3 2 1 2 3 2 3 2 1 2 3 2 3 2 3 2 1 0 1 0 1 0 3 2 3 0 3 2 3 2 3 0 1 2 1 0 3 0 3 2 3 2 3 2 0 1 2 3 0 1 0 1 2 0 1 0 1 0 1 0 1 3 2 1 0 1 0 1 0 1 0 3 2 3 2 1 0 1 0 1 0 3 2 1 0 1 0 1 0 1 3 2 3 0 3 0 1

Sample output 2:
4 4 0.0 0.694045 0.070896 0.235060 0.684412 0.0 0.228137 0.087451 0.105932 0.271192 0.0 0.622876 0.266083 0.060087 0.673830 0.0
4 4 1.0 0.0 0.0 0.0 0.0 0.999980 0.0 0.000020 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0

Sample input 3:
3 3 0.800000 0.100000 0.100000 0.100000 0.800000 0.100000 0.100000 0.100000 0.800000
3 2 0.600000 0.400000 0.400000 0.600000 0.400000 0.600000
1 3 1.000000 0.000000 0.000000
1300 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 0 1 1 1 0 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 0 1 1 0 0 1 1 1 1 1 0 1 1 1 0 1 0 0 0 1 0 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0 1 0 1 1 1 0 1 1 1 0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 1 0 1 1 0 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 1 0 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 0 0 0 1 0 0 0 0 1 0 0 1 1 0 1 0 1 0 0 0 0 1 1 1 0 1 0 0 1 1 0 1 1 1 1 1 1 1 0 0 0 1 1 1 0 1 0 0 1 1 0 1 1 0 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1 0 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 0 0 0 1 0 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 0 1 0 0 1 1 0 1 1 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1 1 0 1 1 1 1 0 0 1 0 1 0 1 1 0 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 0 0 1 0 1 0 1 1 1 1 1 0 1 1 0 1 0 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 0 0 0 0 1 0 1 0 1 0 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 1 0 0 0 1 0 1 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 0 0 1 0 1 1 0 1 1 0 0 0 0 0 0 0 0 1 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 0 1 0 1 1 1 0 0 0 0 0 0 0 1 1 1 0 1 1 1 1 0 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 1 0 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 1 0 0 0 1 1 0 1 1 0 1 1 1 1 1 0 1 0 1 0 0 0 1 1 1 0 0 1 0 0 0 0 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 0 0 0 1 0 0 1 1 0 0 0 1 0 1 1 0 0 0 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0 1 0 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 1 0 1 0 1 0 0 0 1 0 1 0 1 0 1 0 1 1 0 1 0 0 0 1 1 0 1 1 1 1 0 1 0 1

Sample output 3:
3 3 0.845052 0.0774739 0.0774739 0.0841673 0.814073 0.101759 0.0841673 0.101759 0.814073
3 2 0.684137 0.315863 0.126646 0.873354 0.126646 0.873354
 */

