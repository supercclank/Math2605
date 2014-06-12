/*
 * Math Homework - Cory Brzycki
 */
import java.util.Arrays;
import java.util.ArrayList;
public class LUFactorization{
	public static void main(String[] args) {
		//gets n, the first entry in the command line arguments
		int n = Integer.parseInt(args[0]);
		//creates an array that is n by n
		double[][] A = new double[n][n];
		//starting at the 2nd entry of the command line argument array,
		//loops at adds the arguments to array left to right, top to bottom
		int i = 1;
		for (int x = 0; x < n; x++) {
			for (int y = 0; y < n; y++) {
				A[x][y] = Double.parseDouble(args[i]);
				i++;
			}
		}
		//checks to see if matrix is invertible
		if (determinant(A) == 0) {
			System.out.println("Singular matrix (non-invertible), program terminating");
			return;
		}
		//prints out original matrix
		System.out.println("The original matrix...");
		print(A);
		//pivots A if an extra argument is passed by command line (follows after all ints to be added to matrix)
		if (args.length == n*n+2){
			A = pivot(A);
			System.out.println("The pivoted matrix...");
		    print(A);
		}
	        A = luFactorize(A);
	        System.out.println("The LUFactorized matrix...");
		    print(A);
	}
	
		//prints out matrix
		public static void print(double[][] A){
			for (int x = 0; x < A.length; x++) {
				for (int y = 0; y < A.length; y++) {
					System.out.print(A[x][y]+"    ");
				}
				    System.out.println();
			}
			System.out.println("--------");
		}

		//method that preforms LUfactorization on input Matrix A
		public static double[][] luFactorize(double[][] A){
			int i = 0;
				for (int x = 1; x < A.length; x++) {
					i = 0;
					for (int y = 0; y < x; y++) {
						i++;
						A[x][y]= A[x][y]/A[y][y];
						A = subtract(A, x, y, i);
					}
		    }
		    return A;
		}

		//method that subtracts row2 from row1 starting at y of index of matrix A
		public static double[][] subtract(double[][] A, int row1, int row2, int index){
			for (int y = index; y < A.length; y++) {
				A[row1][y] = A[row1][y] - A[row2][y]*A[row1][index - 1];
			}
			return A;
		}
		
		//method that preforms the pivoting of A
	    public static double[][] pivot(double[][] A){
	    	A = eliminateFractions(A);
	    	A = reOrder(A);
	    	return A;
	    }

	    //method that eliminates any fractions present in matrix A so that
	    //it can be properly pivoted
	    public static double[][] eliminateFractions(double[][] A){
	    	for (int x=0; x <A.length; x++){
	    		for (int y=0; y<A.length; y++){
	    			double fraction = 1/A[x][y];
	    			if (Math.abs(fraction)>1&&!Double.isInfinite(fraction)&&!Double.isNaN(fraction)){
	    				A = multiplyRow(A, x, fraction);
	    			}
	    		}
	    	}
	    	return A;
	    }

	    //method that reOrders rows base on the max in each row
	    public static double[][] reOrder(double[][] A){
	    	double max = 0;
			int index = 0;
			for (int y = 0; y < A.length; y++) {
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

	    //mulitplies row X of A by factor
	    public static double[][] multiplyRow(double[][] A, int x, double factor){
	    	for (int y=0; y <A.length; y++){
	    		A[x][y] = A[x][y]*factor;
	    	}
	    	return A;
	    }

	    //swaps row row1 and row row2 of A
	    public static double[][] swap(double[][] A, int row1, int row2){
	    	double[] temp = new double[A.length];
	    	for (int a = 0; a < A.length; a++) {
	    		temp[a] = A[row1][a];
	    	}
	    	for (int a = 0; a < A.length; a++) {
	    		A[row1][a] = A[row2][a];
	    		A[row2][a] = temp[a];
	    	}
	    		return A;
	    }

	    //overloaded recursive method to see if a matrix is singular, calls itself with extra parameters needed
	    public static double determinant(double[][]A){
	    	return determinant(A, new ArrayList<double[][]>(), new ArrayList<Double>());
	    }

	    //recursive method to see if a matrix is singular
	    public static double determinant(double[][]A, ArrayList<double[][]> matrices, ArrayList<Double> scalars){
	    	if (A.length > 2) {
	    		for (int x = 0; x < A.length; x++ ) {
	    			scalars.add(A[0][x]);
	    		}
	    		for (int z = 0; z< A.length; z++){
	    			double[][] temp = new double[A.length-1][A.length-1];
		    		for (int x = 1; x< A.length; x++){
		    		    int position=0;
		    			for (int y=0; y< A.length; y++){
		    				if (y!=z){
		    					temp[x-1][position] = A[x][y];
		    					position++;
		    				}
		    			}
		    		}
		    		matrices.add(temp);
		    	}
	    	}
	    	else {
	    		double[][] B = A;
	    		return scalars.remove(scalars.size()-1) * (B[0][0]*B[1][1]-B[0][1]*B[1][0]);
	    	}
	    	int total = 0 ;
	    	int z =0;
	    	    while (matrices.size()>0){
	    	       total += determinant(matrices.remove(matrices.size()-1), matrices, scalars);
	    	       z++;
	    	    }
	    	    return total;
	    }
}