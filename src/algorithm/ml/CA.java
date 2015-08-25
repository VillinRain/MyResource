package algorithm.ml;
import Jama.*; 

public class CA{
    /*
     Correspondence Analysis
     */
    private double[] explainRatio_;  //the length of explainRatio is nd, each element represent its explained inertia (in %)
    private double[][] RowCoord_;     //principal coordinates of every row profile in nd dimensions
    private double[][] ColCoord_;    //principal coordinates of every column profile in nd dimensions
   
    
    public CA(double [][] dataTable) { //dataTable is the cross table
    	mainProcess(dataTable, 2);
    }
    
    public CA(double[][] dataTable, int nd) {  //nd is the number of dimensions you want to reduce to,  2 dimensions by default if you omit it
    	mainProcess(dataTable, nd);
    }

    private void mainProcess(double[][] dataTable, int nd) {  
    	double sum =  0.;
    	int nrow = dataTable.length;
    	int ncol = dataTable[0].length;
        for (int i = 0; i < nrow; ++ i) {
            for (int j = 0; j < ncol; ++ j) {
                sum += dataTable[i][j];
            }
        }        
        
        Matrix X = new Matrix(dataTable);  
        Matrix P = X.times(1./sum);     // correspondence matrix
        Matrix one = new Matrix(ncol, 1, 1);
        Matrix r = P.times(one);
        one = new Matrix(1, nrow, 1);    // row vector
        Matrix c = one.times(P);        // row vector
        Matrix Z = diag(r, -.5).times(P.minus(r.times(c))).times(diag(c, -.5));   //standardized matrix
        
        SingularValueDecomposition Z_svd = new SingularValueDecomposition(Z);                    
        Matrix U_sub = Z_svd.getU().getMatrix(0, nrow-1, 0, nd-1);   //nrow x nd subMatrix of U
        Matrix V_sub = Z_svd.getV().getMatrix(0, ncol-1, 0, nd-1);  // ncol x nd subMatrix of V
        // ! call getS() will Throw ArrayIndexOutOfBoundsException.
        // Matrix SV = Z_svd.getS().getMatrix(0, 1, 0, 1);
        double sv[] = Z_svd.getSingularValues();
        
        int rankZ = Math.min(nrow-1, ncol-1);
        double totalInertia = 0.;
        for (int i = 0; i < rankZ; i++){    // the real number of singular values is rank(Z), others are because of the calculation accuracy
        	totalInertia += Math.pow(sv[i], 2.);
        }   
        double[] explainRatio = new double[nd];
        for (int i = 0; i < nd; i++){
        	explainRatio[i] = Math.pow(sv[i], 2.) / totalInertia * 100;
        }
        
        
        Matrix SV = new Matrix(nd, nd);      // Construct an nd x nd  matrix of zeros.
        for (int i = 0; i < nd; i++){  
        	SV.set(i, i, sv[i]);
        }
        Matrix RowCoord = diag(r, -0.5).times(U_sub).times(SV);
        Matrix ColCoord = diag(c, -0.5).times(V_sub).times(SV);
        
        explainRatio_ = explainRatio;
        RowCoord_ = RowCoord.getArray();
        ColCoord_ = ColCoord.getArray();
    }
    
    public double[] getExplainRatio(){
    	return explainRatio_;
    }
   
    public double[][] getRowCoord(){
    	return RowCoord_;
    }
    
    public double[][] getColCoord(){
    	return ColCoord_;
    }


    private Matrix diag(Matrix r, double exp) {
        int row = r.getRowDimension();
        int col = r.getColumnDimension();
        int n = row + col - 1;
        Matrix dr = new Matrix(n, n);
        for (int i = 0; i < row; ++ i) {
            for (int j = 0; j < col; ++ j) {
                dr.set(i+j, i+j, Math.pow(r.get(i, j), exp));
            }
        }
        return dr;
    }

    public Matrix diag(double r[], double exp) {
        int n = r.length;
        Matrix dr = new Matrix(n, n);
        for (int i = 0; i < n; ++ i) {
            dr.set(i, i, Math.pow(r[i], exp));
        }
        return dr;
    }

	public static void printArray(double[][] args){            // to be deleted
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
	public static void printArray(double[] args){               // to be deleted
		for (int i = 0; i < args.length; i++ ){
			System.out.print(args[i] + "\t");
		}
		System.out.print("\n");
	}
          
    public static void main(String args[]) {                    
//        double testArray[][] = {
//            {2685,9252,9769,11465,8618,8291,4978},
//            {926,10061,19884,19747,10050,5930,1485},
//            {49,809,4620,14195,11847,11075,4769},
//            {4,0,5,27,41,56,150}
//        };
    	double[][] testArray = {
                {382985,463310,306785,320389,356162,245672,131874,89525,72592,156362},
                {9071,12693,5389,4900,4351,2847,2527,979,1140,12677},
                {116328,119142,83321,103780,62077,54991,34965,11022,15390,31688},
                {53378,27440,46373,19717,17413,11814,7886,2061,2808,14890},
                {25086,30679,16818,22176,16472,13558,10496,3242,4462,17744},
                {207876,90237,134773,71726,43337,39396,88561,102356,105545,21042},
                {308,825,146,106,223,102,59,16,21,981},
                {85690,18069,76459,9917,10551,5253,57723,26721,24120,6071}
            };
        CA test = new CA(testArray);  // try  CA(testArray, 3)
        System.out.println("======explainRatio(in %)===========");
        printArray(test.getExplainRatio());
        System.out.println("======RowCoord===========");
        printArray(test.getRowCoord());
        System.out.println("======ColCoord===========");
        printArray(test.getColCoord());
        
    }
}
