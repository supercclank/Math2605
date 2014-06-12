/*
 * Math Homework - Iteration - Cory Brzycki
 */
import java.util.ArrayList;
public class Iteration {
	public static void main(String[] args) {
		//sets of test case 1 and 2
		double[][] matrix = new double[3][3];
		double[] solution = new double[3];
		//test 1
		matrix[0][0] = 9; matrix[0][1] = 1; matrix[0][2] = 1;
		matrix[1][0] = 2; matrix[1][1] = 10; matrix[1][2] = 3;
		matrix[2][0] = 3; matrix[2][1] = 4; matrix[2][2] = 11;
		solution[0] = 1; solution[1] = 1; solution[2] = 1;
		System.out.println("Test case 1:");
		print(matrix, solution);
		System.out.println("jacobi "+jacobi(matrix, solution));
		System.out.println("gauss seidel "+gauss(matrix, solution));
		//test 2
		System.out.println("Test case 2:(tridiagonal)");
		matrix[0][0] = 5; matrix[0][1] = 0; matrix[0][2] = 0;
		matrix[1][0] = 0; matrix[1][1] = 9; matrix[1][2] = 0;
		matrix[2][0] = 0; matrix[2][1] = 0; matrix[2][2] = -7;
		solution[0] = -1; solution[1] = 2; solution[2] = 3;
		print(matrix, solution);
		System.out.println("jacobi "+jacobi(matrix, solution));
		System.out.println("gauss seidel "+gauss(matrix, solution));
	}
	//method that preforms gauss seidel by setting up needed U,D,L and calling recursive helper method
	//format the String
	public static String gauss(double[][] matrix, double[] solution) {
		double[][] D = new double[matrix.length][matrix.length];
		double[][] L = new double[matrix.length][matrix.length];
		double[][] U = new double[matrix.length][matrix.length];
		double[] X = new double[solution.length];
		for (int x = 0; x <X.length; x++) {
			X[x] = 1;
		}
		for (int x = 0; x < matrix.length; x++) {
			for (int y = 0; y < matrix.length; y++) {
				if (y==x) {
			        D[x][y] = matrix[x][y];
			    }
			    if (y<x) {
			    	L[x][y] = matrix[x][y];
			    }
			    if (y>x) {
			    	U[x][y] = matrix[x][y];
			    }
		    }
	    }
	    double[][] T = add(D,L);
	    T = invert(T);
	    double[] C = multiply(T,solution);
	    for (int x = 0; x < T.length; x++) {
	    	for (int y = 0; y < T.length; y++) {
	    		T[x][y] = -T[x][y];
	    	}
	    }
	    T = multiply(T,U);
	    X = gaussRecurse(T,C,X,matrix,solution);
	    return "a= "+X[0]+" b= "+X[1]+" c= "+X[2];
    }
    //method that inverts a matrix using the cofactor method
    public static double[][] invert(double[][] A) {
    	//get the cofactor matrix
    	double[][] cofactors = getCofactorMatrix(A);
    	//gets the determinant of c
    	double c = determinant(A, new ArrayList<double[][]>(), new ArrayList<Double>(), 1);
    	double[][] B = new double[A.length][A.length];
    	//transpose the matrix formed by cofactors
    	for (int x = 0; x < A.length; x++) {
    	    for (int y = 0; y < A.length; y++) {
    	    B[y][x] = cofactors[x][y]/c;
    	}	
    	}
    	return B;
    }
    //returns a cofactor matrix formed by determinants of partitions
    public static double[][] getCofactorMatrix(double[][] A){
    	double[][] cofactors = new double[A.length][A.length];
    	int sign = 1;
    	//fix this loop D:
    	for (int z1 = 0; z1 < A.length; z1++) {
    		double[][] B = new double[2][2];
    		for (int z2 = 0; z2 < A.length; z2++) {
    		int xx = 0;
    		int yy = 0;
	    	for (int x = 0; x < A.length; x ++) {
	    		for (int y = 0; y < A.length; y++) {
	    			if (x != z1 && y != z2 && xx != 2) {
	    			    B[xx][yy] = A[x][y];
	    			    yy++;
	    			    if (yy == 2) {
	    			    	xx++;
	    			    	yy = 0;
	    			    }
	    			}
	    		}
	    	}
    	cofactors[z1][z2] = (sign  * ((B[0][0]*B[1][1]-B[0][1]*B[1][0])));
    	sign = sign * -1;
    	}
    	}
    	return cofactors;
    }
    //recursive helper method for gauss seidel that calls itself until ||Ax-b|| < 0.00000000000001
    public static double[] gaussRecurse(double[][] T, double[] C,double[] X, double[][] A, double[] B) {
    	X = multiply(T,X);
    	X = add(X,C);
	    if (subtractWithAbs(multiply(A,X),B) < 0.00000000000001) {
	    	return X;
	    }
    	return gaussRecurse(T,C,X,A,B);
    }
	//method that preforms jacobi iteration by setting up needed U,D,L and calling recursive helper method
	//format the String
	public static String jacobi(double[][] matrix, double[] solution) {
		double[][] D = new double[matrix.length][matrix.length];
		double[][] L = new double[matrix.length][matrix.length];
		double[][] U = new double[matrix.length][matrix.length];
		double[] X = new double[solution.length];
		for (int x = 0; x <X.length; x++) {
			X[x] = 1;
		}
		for (int x = 0; x < matrix.length; x++) {
			for (int y = 0; y < matrix.length; y++) {
				if (y==x) {
			        D[x][y] = -1/matrix[x][y];
			    }
			    if (y<x) {
			    	L[x][y] = matrix[x][y];
			    }
			    if (y>x) {
			    	U[x][y] = matrix[x][y];
			    }
		    }
	    }
	    double[][] T = add(L,U);
	    T = multiply(D,T);
	    double[] C = multiplyNeg(D,solution);
	    X = jacobiRecurse(T,C,X, matrix, solution);
	    return "a= "+X[0]+" b= "+X[1]+" c= "+X[2];
    }
    //recursive helper method for jacobi that calls itself until ||Ax-b|| < 0.00000000000001
    public static double[] jacobiRecurse(double[][] T,double[] C,double[] X, double[][] A, double[] B) {
    	X = multiply(T,X);
    	X = add(X,C);
	    if (subtractWithAbs(multiply(A,X),B) < 0.00000000000001) {
	    	return X;
	    }
    	return jacobiRecurse(T,C,X,A,B);
    }
    //multiplies A by B, returns AB
    public static double[] multiply(double[][] A, double[] B){
		double[] result = new double[A.length];
		for (int x = 0; x < result.length; x++) {
			for (int y = 0; y < result.length; y++) {	
				double total = 0;
				for (int z = 0; z < result.length; z++) {
					total = total + A[x][z] * B[z]; 
				}
				result[x] = total;
			}
		}
		return result;
	}
	//multiplies A by vector B, returns -AB
    public static double[] multiplyNeg(double[][] A, double[] B){
		double[] result = new double[B.length];
		for (int x = 0; x < result.length; x++) {
			for (int y = 0; y < result.length; y++) {	
				double total = 0;
				for (int z = 0; z < result.length; z++) {
					total = total + -1 * A[x][z] * B[z]; 
				}
				result[x] = total;
			}
		}
		return result;
	}
    //multiplies A by B, returns AB
    public static double[][] multiply(double[][] A, double[][] B){
		double[][] result = new double[A.length][A.length];
		for (int x = 0; x < result.length; x++) {
			for (int y = 0; y < result.length; y++) {	
				double total = 0;
				for (int z = 0; z < result.length; z++) {
					total = total + A[x][z] * B[z][y]; 
				}
				result[x][y] = total;
			}
		}
		return result;
	}
	//adds matrices A and B
    public static double[][] add(double[][] A, double[][] B){
    	for (int x = 0; x < A[0].length; x++){
			for (int y = 0; y < A[0].length; y++) {
				A[x][y] = A[x][y] + B[x][y];
			}
		}
			return A;
		}
	//adds vectors A and B
	public static double[] add(double[] A, double[] B){
    	for (int x = 0; x < A.length; x++){
				A[x] = A[x] + B[x];
		}
			return A;
	}
	//returns ||A-B||
	public static double subtractWithAbs(double[] A, double[] B){
		double total = 0;
    	for (int x = 0; x < A.length; x++){
				total += Math.pow(A[x] - B[x], 2);
		}
			return Math.sqrt(total);
	}	
	//recursive method that returns a matrix's determinant
	public static double determinant(double[][]A, ArrayList<double[][]> matrices, ArrayList<Double> scalars, int sign){
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
	    		double hold = scalars.remove(scalars.size()-1);
	    		return  sign * hold * ((B[0][0]*B[1][1]-B[0][1]*B[1][0]));
	    	}
	    	int total = 0 ;
	    	int z =0;
	    	    while (matrices.size()>0){
	    	       total += determinant(matrices.remove(matrices.size()-1), matrices, scalars, sign);
	    	       sign = sign * -1;
	    	       z++;
	    	    }
	    	    return total;
	    }
	    //print matrix
	    public static void print(double[][] A, double[] B){
	    	String[] b = new String[]{"a","b","c"};
	    	for(int x=0; x<A.length; x++){
	    		for(int y=0; y<A.length; y++) {
	    			System.out.print(A[x][y]+" ");
	    		}
	    		System.out.println(b[x]+" "+B[x]);
	    	}
	    }
}