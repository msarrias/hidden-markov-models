import java.util.Scanner;

public class HMM2 {
    //Estimate Sequence of States
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

    public static String array2Str(int [] path){
        String[] sr = new String[path.length];
        for (int i = 0; i < path.length; i++){
            sr[i] = String.valueOf(path[i]);
        }
        String bestPathString = String.join(" ", sr);
        return bestPathString;
    }


    public static int[] ViterbiAlgorithm(double[][] A, double[][] B, double[][] pi, int[] O){
        int nHiddenStates = A[0].length; //N of hidden states
        int nEmissions = O.length; //N of emissions
        double max;
        int maxIdx;
        double[][] delta = new double[nEmissions][nHiddenStates];
        int [][] deltaIndex = new int[nEmissions][nHiddenStates];
        int [] bestPath = new int[O.length];

        for (int t = 0; t < nEmissions; t++){
            if (t == 0){
                for (int i = 0; i < nHiddenStates; i++){
                    delta[0][i] = Math.log(pi[0][i] * B[i][O[0]]);
                }
            } else {
                for (int i = 0; i < nHiddenStates; i++){
                    maxIdx = 0;
                    max = -16;
                    for (int j = 0; j < nHiddenStates; j++){
                        double temp = delta[t-1][j] + Math.log(A[j][i]) + Math.log(B[i][O[t]]);
                        if (temp > max){
                            maxIdx = j;
                            max = temp;
                        }
                    }
                    delta[t][i] = max;
                    deltaIndex[t][i] = maxIdx;
                }
            }
        }
        for (int t = nEmissions - 1; t >= 0; t--){
            if (t == nEmissions - 1){
                max = -16;
                for (int i = 0; i < nHiddenStates; i++){
                    if (delta[t][i] > max){
                        max = delta[t][i];
                        bestPath[t] = i;
                    }
                }
            } else bestPath[t] = deltaIndex[t + 1][bestPath[t + 1]];
        }
        return bestPath;
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

        int[] bestPath = ViterbiAlgorithm(A, B, pi, O);
        System.out.println(array2Str(bestPath));

    }
}

/* Sample input:
4 4 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.8 0.1 0.1 0.0
4 4 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.1 0.0 0.0 0.9
1 4 1.0 0.0 0.0 0.0
4 1 1 2 2

Sample Output:
0 1 2 1
 */