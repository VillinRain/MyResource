package algorithm.ml;
import Jama.*; 

import java.lang.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



public class MCA {
	/*
	  Multiple Correspondence Analysis
	 
	  Computations based on the Burt matrix(a square matrix)
	  Adjust Inertias 
	 */
    private double[] adjExplainRatio_;         //the length of adjExplainRatio is nd, each element represent its explained inertia (in %)
    private double[][] Coord_;            //principal coordinates in nd dimensions, every row represent a point
    
    public MCA(double BurtMatrix[][], int nvars) {  //Burt matrix,  nvars is the the number of variables
    	mainProcess(BurtMatrix, nvars, 2);
    }
    
    public MCA(double BurtMatrix[][], int nvars, int nd) {   //nd is the number of dimensions you want to reduce to,  2 dimensions by default if you omit it
    	mainProcess(BurtMatrix, nvars, nd);
    }

    private void mainProcess(double Burt[][], int nvars, int nd) {  
    	double sum =  0.;
    	int nrow = Burt.length;
    	int ncol = Burt[0].length;
        for (int i = 0; i < nrow; ++ i) {
            for (int j = 0; j < ncol; ++ j) {
                sum += Burt[i][j];
            }
        }        
        
    	Matrix X = new Matrix(Burt);
    	Matrix P = X.times(1./sum);   // correspondence matrix
    	Matrix one = new Matrix(ncol, 1, 1);
        Matrix c = P.times(one);     // column profile
        Matrix Z = diag(c, -.5).times(P.minus(c.times(c.transpose()))).times(diag(c, -.5));   //standardized matrix
        
        EigenvalueDecomposition Z_eign = new EigenvalueDecomposition(Z);
        Matrix EignValue = Z_eign.getD();
        Matrix EignVector_sub = Z_eign.getV().getMatrix(0, ncol-1, 0, nd-1);  // ncol x nd subMatrix of EignVector
        
        int rankZ = ncol - nvars;   //rank(Z) <= ncol -nvars
        double totalInertia = 0.;
        List svAdj = new ArrayList();
        for (int i = 0; i < rankZ; i++){    // the real number of singular values is rank(Z), others are because of the calculation accuracy
        	double thisEignValue = EignValue.get(i, i);
        	totalInertia += Math.pow(thisEignValue, 2.);
        	if (thisEignValue >= 1./nvars){
        		double adjValue = Math.pow((double) nvars / (nvars-1) * (thisEignValue - 1./nvars), 2.);
        		svAdj.add(adjValue);
        	}
        } 
        double totalInertiaAdj = ((double) nvars / (nvars-1)) * (totalInertia - (double) (ncol-nvars) / Math.pow(nvars, 2.0));
        double[] adjExplainRatio = new double[nd];
        for (int i = 0; i < nd; i++){
        	adjExplainRatio[i] = (double) svAdj.get(i) / totalInertiaAdj * 100;
        }
        
        Matrix EignValue_sub = EignValue.getMatrix(0, nd-1, 0, nd-1);     //   nd x nd subMatirx of EignValue 
        Matrix Coord = diag(c, -0.5).times(EignVector_sub).times(EignValue_sub);
        
        adjExplainRatio_ = adjExplainRatio; 
        Coord_ = Coord.getArray();         
    }
    
    public double[] getAdjExplainRatio(){
    	return adjExplainRatio_;
    }
    
    public double[][] getCoord(){
    	return Coord_;
    }
    
    public static  Matrix diag(Matrix r, double exp) {
        int nrow = r.getRowDimension();
        int ncol = r.getColumnDimension();
        int n = nrow + ncol - 1;
        Matrix dr = new Matrix(n, n);
        for (int i = 0; i < nrow; ++ i) {
            for (int j = 0; j < ncol; ++ j) {
                dr.set(i+j, i+j, Math.pow(r.get(i, j), exp));
            }
        }
        return dr;
    }       
              
	public static void printArray(double[][] args){              // to be deleted
		int nrow = args.length;
		int ncol = args[0].length;
		for(int i=0; i<nrow; i++){
			for(int j=0; j<ncol; j++){
				System.out.print(args[i][j] + "\t");
			}
			System.out.print("\n");
		}
		System.out.print("\n");
	}
	public static void printArray(double[] args){                // to be deleted
		for (int i = 0; i < args.length; i++ ){
			System.out.print(args[i] + "\t");
		}
		System.out.print("\n");
	}
            
    public static void main(String args[]) {
//        double testArray[][] = {
//        		{ 119,0,0,0,0,27,28,30,22,12,49,40,18,7,5,15,25,17,34,28 },
//        		{ 0,322,0,0,0,38,74,84,96,30,67,142,60,41,12,22,102,76,68,54 },
//        		{ 0,0,204,0,0,3,48,63,73,17,18,75,70,34,7,10,44,68,58,24 },
//        		{ 0,0,0,178,0,3,21,23,79,52,16,50,40,56,16,9,52,28,54,35 },
//        		{ 0,0,0,0,48,0,3,5,11,29,2,9,9,16,12,4,9,13,12,10 },
//        		{ 27,38,3,3,0,71,0,0,0,0,43,19,4,3,2,9,17,10,10,25 },
//        		{ 28,74,48,21,3,0,174,0,0,0,36,88,34,15,1,16,51,42,45,20 },
//        		{ 30,84,63,23,5,0,0,205,0,0,37,90,57,19,2,10,53,63,51,28 },
//        		{ 22,96,73,79,11,0,0,0,281,0,27,88,75,74,17,6,66,70,92,47 },
//        		{ 12,30,17,52,29,0,0,0,0,140,9,31,27,43,30,19,45,17,28,31 },
//        		{ 49,67,18,16,2,43,36,37,27,9,152,0,0,0,0,25,24,15,38,50 },
//        		{ 40,142,75,50,9,19,88,90,88,31,0,316,0,0,0,15,97,67,89,48 },
//        		{ 18,60,70,40,9,4,34,57,75,27,0,0,197,0,0,5,51,83,41,17 },
//        		{ 7,41,34,56,16,3,15,19,74,43,0,0,0,154,0,6,44,30,51,23 },
//        		{ 5,12,7,16,12,2,1,2,17,30,0,0,0,0,52,9,16,7,7,13 },
//        		{ 15,22,10,9,4,9,16,10,6,19,25,15,5,6,9,60,0,0,0,0 },
//        		{ 25,102,44,52,9,17,51,53,66,45,24,97,51,44,16,0,232,0,0,0 },
//        		{ 17,76,68,28,13,10,42,63,70,17,15,67,83,30,7,0,0,202,0,0 },
//        		{ 34,68,58,54,12,10,45,51,92,28,38,89,41,51,7,0,0,0,226,0 },
//        		{ 28,54,24,35,10,25,20,28,47,31,50,48,17,23,13,0,0,0,0,151 }};
//        double testArray[][] = {
//        						{8,0,5,3,1,7},
//        						{0,2,1,1,1,1},
//        						{5,1,6,0,2,4},
//        						{3,1,0,4,0,4},
//        						{1,1,2,0,2,0},
//        						{7,1,4,4,0,8}};
//        double testArray[][] = {
//				{31287,2908,3275,9071,12693,5389},
//				{2908,141361,11569,116328,119142,83321},
//				{3275,11569,64986,53379,27441,46394},
//				{9071,116328,53379,165513,111378,119267},
//				{12693,119142,27441,111378,145423,74673},
//				{5389,83321,46374,119267,74673,123634}};
        double testArray[][] = {
				{31287,0,0,9071,22216,12693,18594,5389,25898},
				{0,141361,0,116328,25033,119142,22219,83321,58040},
				{0,0,64986,53379,11607,27441,37545,46394,18592},
				{9071,116328,53379,165513,0,111378,54135,119267,46246},
				{22216,25033,11607,0,72121,34045,38076,4367,67754},
				{12693,119142,27441,111378,34045,145423,0,74673,48961},
				{18594,22219,37545,54135,38076,0,92211,48961,65039},
				{5389,83321,46374,119267,4367,74673,48961,123634,0},
				{25898,58040,18592,46246,67754,48961,65039,0,114000}};
        MCA test = new MCA(testArray, 4);    // try  MCA(testArray, 4, 3)
        System.out.println("======adjExplainRatio(in %)===========");
        printArray(test.getAdjExplainRatio());
        System.out.println("======coordinates===========");
        printArray(test.getCoord());

    }
}
