/*    
 * Math Homework - Least Squares - Cory Brzycki
 */
public class LeastSquares {
	private static double[][] dataSet;
	public static void main(String[] args) {
		//adds the data set to the Array
		dataSet = new double[2][5];
		dataSet[0][0] = 0;
		dataSet[1][0] = 1.5;
		dataSet[0][1] = 1; 
		dataSet[1][1] = 2.5;
		dataSet[0][2] = 2;
		dataSet[1][2] = 3.5;
		dataSet[0][3] = 3;
		dataSet[1][3] = 5;
		dataSet[0][4] = 4;
		dataSet[1][4] = 7.5;
		//Method call to get fit to line
	    System.out.println("fit to line "+fitToLine(dataSet));
	    //method call to get fit to paraboloa 
	    System.out.println("fit to parabola "+fitToParabola(dataSet));
	}
	//method that uses a dataset to fit to a ax^2+bx+c parabola
	public static String fitToParabola(double[][] dataSet) {
		double[] y = new double[dataSet[0].length];
		double[][] A = new double[dataSet[0].length][dataSet.length+1];
		//creates the matrix in the needed form to begin least squares (col x^2s, xs then 1s)
		for (int x = 0 ; x < dataSet[0].length; x++) {
			y[x] = dataSet[1][x];
			A[x][0] = dataSet[0][x] * dataSet[0][x];
			A[x][1] = dataSet[0][x];
			A[x][2] = 1;
		}
		//gets AtY and stores in B
		double[] B = getAtA(A, y);
		//gets AtA and stores back in A
		A = getAtA(A);
		//solves AtA for AtY
		A = solve(A,B);
		//last section formats the result nicely
		String line = "";
		for (int x = 0; x < A.length; x++) {
			line = line + A[x][A[0].length-1] + " " ;
		}
		String line1 = line.substring(0, line.indexOf(" "));
		line = line.substring(line.indexOf(" ")+1);
		String line2 = line.substring(0, line.indexOf(" "));
		line = line.substring(line.indexOf(" ")+1);
		String line3 = line.substring(0, line.indexOf(" "));
		return line1+"X^2 "+line2+"X "+line3;
	}
	//method that uses a dataset to fit to a ax+b line
	public static String fitToLine(double[][] dataSet) {
		double[] y = new double[dataSet[0].length];
		double[][] A = new double[dataSet[0].length][dataSet.length];
		//creates the matrix in the needed form to begin least squares (col x then 1s)
		for (int x = 0 ; x < dataSet[0].length; x++) {
			y[x] = dataSet[1][x];
			A[x][0] = dataSet[0][x];
			A[x][1] = 1;
		}
		//gets AtY and stores in B
		double[] B = getAtA(A, y);
		//gets AtA and stores back in A
		A = getAtA(A);
		//solves AtA for AtY
		A = solve(A,B);
		String line = "";
		for (int x = 0; x < A.length; x++) {
			line = line + A[x][A[0].length-1] + " " ;
		}
		//last section formats the result nicely
		line = line.substring(0, line.indexOf(" "))+"X"+line.substring(line.indexOf(" "));
		return line;
	}
	//solves a matris using gaussian elimination
	public static double[][] solve (double[][] A, double[] Y) {
		double factor = 0;
		//adds y to the end of A to make augmented matrix
		double[][] B = new double[A.length][A[0].length+1];
		for (int x = 0; x < B.length; x++) {
			B[x][B[0].length - 1] = Y[x]; 
			for (int y = 0; y <B.length; y++) {
				B[x][y] = A[x][y];
			}
		}
		//reorders B to ensure better results
		B = reOrder(B);
		for (int x = 0; x <B.length; x++) {
				double divisor = B[x][x];
			for (int z = 0; z <B.length; z++) {
				if (z != x) {
					factor = B[z][x]/divisor;
					B = subtract(B, z,x, factor);
				}
				
			}
		}
		//reirders B to find final result
		B = reOrder(B);
		//stores reults in the correct collumn of B, rounding to 11 sigFigs to account for inpercision
		for (int x = 0; x < B.length; x++) {
			B[x][B[0].length- 1] = B[x][B[0].length-1] / B[x][x];
			B[x][B[0].length- 1] = (double)Math.round(B[x][B[0].length- 1] * 1000000000) / 1000000000;
			B[x][x] = 1;
		}
		return B;
	}
	//swaps row row1 and row row2 of A
	public static double[][] swap(double[][] A, int row1, int row2){
	    	double[] temp = new double[A[0].length];
	    	for (int a = 0; a < A[0].length; a++) {
	    		temp[a] = A[row1][a];
	    	}
	    	for (int a = 0; a < A[0].length; a++) {
	    		A[row1][a] = A[row2][a];
	    		A[row2][a] = temp[a];
	    	}
	    		return A;
	    }
	//method that reorders to max diagonal
	public static double[][] reOrder(double[][] A){
	    	double max = 0;
			int index = 0;
			for (int y = 0; y < A[0].length; y++) {
				index = 0;
				max = 0;
				for (int x = y; x < A.length; x++) {
					if (Math.abs(A[x][y]) > max) {
						max = A[x][y];
						index = x;
						swap(A, x, y);
					}
				}
			}
			return A;
	    }
	//multiplies AtY and returns resulting matrix
	public static double[] getAtA (double[][] A, double[] Y) {
		double[] result = new double[A[0].length];
		for (int x = 0; x < result.length; x++) {
			for (int y = 0; y < result.length - 1; y++) {	
				double total = 0;
				for (int z = 0; z < A.length; z++) {
					total = total + A[z][x] * Y[z]; 
				}
				result[x] = total;
			}
		}
		return result;
	}
    //multiplies AtY and returns resulting matrix
	public static double[][] getAtA(double[][] A){
		double[][] result = new double[A[0].length][A[0].length];
		for (int x = 0; x < result.length; x++) {
			for (int y = 0; y < result.length; y++) {	
				double total = 0;
				for (int z = 0; z < A.length; z++) {
					total = total + A[z][y] * A[z][x]; 
				}
				result[x][y] = total;
			}
		}
		return result;
	}
    //method that subtracts row2*factor from row1
	public static double[][] subtract(double[][] A, int row1, int row2, double factor){
			for (int y = 0; y < A[0].length; y++) {
				A[row1][y] = A[row1][y] - factor*A[row2][y];
			}
			return A;
		}
}