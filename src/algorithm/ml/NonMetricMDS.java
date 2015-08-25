package algorithm.ml;
import Jama.*; 

import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NonMetricMDS {
	/*
	  Non-metric multidimensional scaling
	 * */
	
	private double[][] configuration;   ////configuration has nd column, every row represent the coordinates of a point 
	private double stress_;     //stress value, less the value is, better the effect is
	
	public NonMetricMDS(double[][] DissimilarityMatrix){  //Dissimilarity is a square matrix, the diagonal elements of which are all 0
		init(DissimilarityMatrix, 2, false);
	}
	
	public NonMetricMDS(double[][] DissimilarityMatrix, int nd){ //nd is the number of dimensions you want to reduce to,  2 dimensions by default if you omit it
		init(DissimilarityMatrix, nd, false);
	}
	
	public NonMetricMDS(double[][] DissimilarityMatrix, int nd, boolean showDetailsOption){   //by showDetailsOption(true or false), you can decide whether show the details or not
		init(DissimilarityMatrix, nd, showDetailsOption);
			
	}
	
	public void init(double[][] DissimilarityMatrix, int nd, boolean showDetailsOption){
		CMDScale CMDClass = new CMDScale(DissimilarityMatrix, nd);
		mainProcess(CMDClass.getConfiguration(), UpTriangular(DissimilarityMatrix), showDetailsOption);
		
	}
	
	public void mainProcess(double[][] Y, double[] dissimilarities, boolean showDetailsOption){
		double TolX = 1E-6;
		double TolFun = 1E-6;
		double MaxIter = 100;
		int nrow = Y.length;
		int ncol = Y[0].length;
		double stress;
		Matrix YMatrix = new Matrix(Y);
		
		Matrix OneRowMatrix = new Matrix (1, nrow, 1.);
		Matrix OneColMatrix = new Matrix (ncol, 1, 1.);
		
		int iter = 0;
		boolean resetCG = true;
		double oldStress = Double.MAX_VALUE;
		double stepLen = Double.MAX_VALUE;
		Matrix OldGradMatrix = new Matrix(nrow, ncol);
		Matrix StepDirMatrix = new Matrix(nrow, ncol);
		double oldNormGrad = Double.MAX_VALUE;                 //init
		while (true){
			// Center YMatrix
			Matrix colAvg = (new Matrix(1, nrow, 1./nrow)).times(YMatrix);
			for (int i = 0; i < nrow; i++){
				Matrix YMatrixSub = YMatrix.getMatrix(i, i, 0, ncol-1).minus(colAvg);
				YMatrix.setMatrix(i, i, 0, ncol-1, YMatrixSub);
			}
	
			double[] distances = pdist(YMatrix.getArray());
			// scale 
			double scale = max(dissimilarities) / max(distances);
			YMatrix.timesEquals(scale);
		
			int distLen = distances.length;
			for (int i = 0; i < distLen; i++){
				distances[i] = distances[i] * scale;
			}
			double[] disparities = lsqisotonic(dissimilarities, distances);
			Map m =  stressCrit(YMatrix.getArray(), disparities);
			stress =  (double) m.get("s");
			Matrix GradMatrix = new Matrix((double[][]) m.get("grad"));
			double normGrad = GradMatrix.normF();
			
			if (showDetailsOption){
				if (iter == 0){
					System.out.printf("iter: %6d      stress: %12g      Norm of Gradient: %12g\n", iter, stress, normGrad);
				} else {
					System.out.printf("iter: %6d      stress: %12g      Norm of Gradient: %12g      stepLen: %12g\n", iter, stress, normGrad, stepLen);
				}
			}

			if (stress < TolFun) {
				if (showDetailsOption){
					System.out.println("Iterations terminated:  stress less than the TolFun.");
				}
				break;				
			} else if (normGrad < TolFun * stress) {
				if (showDetailsOption){
					System.out.println("Iterations terminated:  Relative norm of gradient less than TolFun.");
				}
				break;
			} else if ((oldStress - stress) < TolFun * stress) {
				if (showDetailsOption){
					System.out.println("Iterations terminated:  Relative change in criterion less than TolFun.");
				}
				break;
			} else if (stepLen < TolX * YMatrix.normF()) {
				if (showDetailsOption){
					System.out.println("Iterations terminated:  Norm of change in configuration less than TolX.");
				}
				break;				
			} else if (iter == MaxIter) {
				if (showDetailsOption){
					System.out.println("Iteration limit exceeded");
				}
				break;
			}
			
			if (resetCG) {
				resetCG = false;
				StepDirMatrix = GradMatrix.times(-1.);
			} else {				
				double beta = OneRowMatrix.times(GradMatrix.minus(OldGradMatrix).arrayTimes(GradMatrix)).times(OneColMatrix).get(0, 0) / Math.pow(oldNormGrad, 2.);
				beta = Math.max(beta, 0.);
				StepDirMatrix = GradMatrix.times(-1.).plus(StepDirMatrix.times(beta));
			}
			oldStress = stress;
			OldGradMatrix = GradMatrix;
			oldNormGrad = normGrad;
			
			double maxStepLen = 2.;
			while (true) {
				stress = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(maxStepLen)).getArray(), disparities).get("s");
				if (stress > oldStress){
					break;
				}
				maxStepLen = 2 * maxStepLen;					
			}
			
			m = goldSearch(0, maxStepLen, YMatrix.getArray(), StepDirMatrix.getArray(), disparities);
			double alpha = (double) m.get("x");
			stress = (double) m.get("fx");

			boolean jumpOut = false;
			if (stress > oldStress) {
				while (true) {
					alpha = alpha / 2;
					if (alpha <= 1E-12){
						jumpOut = true;
						break;
					}
					stress = (double)stressCrit(YMatrix.plus(StepDirMatrix.times(alpha)).getArray(), disparities).get("s");
					if (stress < oldStress){
						break;
					}
				}
				resetCG = true;				
			}
			
			if (jumpOut) {
				stress = oldStress;
				break;
			}
			
			Matrix Ystep = StepDirMatrix.times(alpha);
			stepLen = StepDirMatrix.norm2() * alpha;
			YMatrix = YMatrix.plus(Ystep);
			

			iter += 1;

		}
		configuration = YMatrix.getArray();
		stress_ = stress;
	}
	
	public static Map stressCrit(double[][] Y, double[] disparities){
		Map map = new HashMap();
		double[] distances = pdist(Y);
		Matrix dist = new Matrix(distances, distances.length);
		Matrix disp = new Matrix(disparities, disparities.length);
		
		Matrix YMatrix = new Matrix(Y);
		Matrix diffs = dist.minus(disp);

		double sumDiffSq = diffs.arrayTimes(diffs).norm1();
		double sumDistSq = dist.arrayTimes(dist).norm1();
		double S = Math.sqrt(sumDiffSq / sumDistSq);
		
		int nrow = Y.length;
		int ncol = Y[0].length;
		Matrix grad = new Matrix(nrow, ncol);
		
		if (sumDiffSq > 0) {
			if (allGreatThanZero(distances) == 1){
				Matrix dS = (diffs.times(1./sumDiffSq).minus(dist.times(1./sumDistSq))).arrayRightDivide(dist);
				dS = new Matrix(squareForm(dS.getColumnPackedCopy()));
				
				int[] repcols = new int[nrow];
				for (int i = 0; i < nrow; i++){
					repcols[i] = -1;
				}
				for (int j = 0; j < ncol; j++){
					for (int k = 0; k < nrow; k++){
						repcols[k] += 1;
					}
					Matrix dY = YMatrix.getMatrix(0, nrow-1, repcols);
					dY = dY.minus(dY.transpose());
					Matrix onecol = new Matrix(nrow, 1, 1.0);
					grad.setMatrix(0, nrow-1, j, j, dS.arrayTimes(dY).times(onecol).times(S));
				} 
				
			}else{
				System.out.println("Points in the configuration have co-located, use a different method");  //this is a error because of your data
			}
		}
	
		map.put("s", S);
		map.put("grad", grad.getArray());
		return map;
	}
	
	public static double[] lsqisotonic(double[] x, double[] y){   // isotonic regression		
		Map m = ascendingSort(x, y);
		double[] yhat = (double[]) m.get("ySorted");
		int[] ord = (int[]) m.get("indexSorted");
		List<Double> weightList = new ArrayList<Double>();
		List<Double> yhatUniqueList = new ArrayList<Double>();
		List<Integer> cumNumList = new ArrayList<Integer>();

		int n = x.length;    
		double[] w = new double[n];
		for (int i = 0; i < n; i++){
			w[i] = 1.;
		}
		
		yhatUniqueList.add(0.);           // fill the 0th position, is not in use
		yhatUniqueList.add(yhat[0]);
		weightList.add(0.);               // fill the 0th position, is not in use
		weightList.add(w[0]); 
		cumNumList.add(0);
		cumNumList.add(1);
		
		int j = 1;
		for (int i = 1; i < n; i++){
			j += 1;
			if (yhatUniqueList.size() == j){
				yhatUniqueList.add(yhat[i]);
				weightList.add(w[i]);
			} else {
				yhatUniqueList.set(j, yhat[i]);
				weightList.set(j, w[i]);
			}
			double yhatUniqueList_j = yhatUniqueList.get(j);
			double yhatUniqueList_j_1 = yhatUniqueList.get(j-1);
			while (yhatUniqueList_j < yhatUniqueList_j_1 && j > 1){
				double weightList_j = weightList.get(j);
				double weightList_j_1 = weightList.get(j-1);
				double temp = weightList_j + weightList_j_1;			
				weightList.set(j-1, temp);
				temp = (weightList_j * yhatUniqueList_j + weightList_j_1 * yhatUniqueList_j_1) / temp;
				yhatUniqueList.set(j-1, temp);
				j -= 1;
				yhatUniqueList_j = yhatUniqueList.get(j);
				yhatUniqueList_j_1 = yhatUniqueList.get(j-1);
			}
			if (cumNumList.size() == j){
				cumNumList.add(i+1);
			} else {
				cumNumList.set(j, i+1);
			}
		}
		
		for (int k = 1; k <= j; k++) {
			for (int i = cumNumList.get(k-1); i < cumNumList.get(k); i++){
				yhat[i] = yhatUniqueList.get(k);
			}
		}
		
		double[] newYhat = new double[n];
		for (int i = 0; i < n; i++){
			newYhat[ord[i]] = yhat[i];
		}
	
		return newYhat;				
	}
	
	public static Map goldSearch(double a, double b, double[][] Y, double[][] stepDir, double[] disparities){   // golden section search
		Map map = new HashMap();
		Matrix YMatrix = new Matrix(Y);
		Matrix StepDirMatrix = new Matrix(stepDir);
		
		double Tolx = 1E-6;
		int MaxIter = 500;
		double lambda = (Math.sqrt(5) - 1) / 2;
		double a1 = b - lambda * (b - a);
		double y1 = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(a1)).getArray(), disparities).get("s");
		double a2 = a + lambda * (b - a);
		double y2 = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(a2)).getArray(), disparities).get("s");
		
		int iter = 0;
		while (true){
			if (y1 >= y2){
				a = a1;
				a1 = a2;
				y1 = y2;
				a2 = a + lambda * (b - a);
				y2 = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(a2)).getArray(), disparities).get("s");
			} else{
				b = a2;
				a2 = a1;
				y2 = y1;
				a1 = b - lambda * (b - a);
				y1 = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(a1)).getArray(), disparities).get("s");
			}
			
			iter += 1;
			if (iter == MaxIter || (b-a) < Tolx){             // || Math.abs(y2 - y1)< Tolf
				break;
			} 
		}
		
		double x = (a + b) / 2;
		double fx = (double) stressCrit(YMatrix.plus(StepDirMatrix.times(x)).getArray(), disparities).get("s");
		map.put("x", x);
		map.put("fx", fx);
		return map;
	}
	
	public static Map ascendingSort(double[] xx, double[] yy){

		double[] x = xx.clone();
		double[] y = yy.clone();
		Map map = new HashMap();
		
		int n = xx.length;
		int[] index = new int[n];
		for (int i = 0; i < n; i++){
			index[i] = i;
		}
		
		for (int i = 0; i < n; i++){
			int k = i;                       //the index of the minimum of x, if equal, compare y
			for (int j = i+1; j < n; j++){
				if (x[k] > x[j]){
					k = j;
				} else if(Math.abs(x[k] - x[j]) < 1E-10){           // x[k] == x[j]
					if (y[k] > y[j]){
						k = j;
					}
				}
			}
			if (k != i){
				double temp = x[i];
				x[i] = x[k];
				x[k] = temp;
				
				temp = y[i];
				y[i] = y[k];
				y[k] = temp;
				
				int temp_index = index[i];
				index[i] = index[k];
				index[k] = temp_index;
			}
		}
		
		map.put("xSorted", x);
		map.put("ySorted", y);
		map.put("indexSorted", index);
		return map;
	}
	
	public static double[] pdist(double[][] Y){
		int nrow = Y.length;
		int ncol = Y[0].length;
		double[] distances = new double[nrow * (nrow - 1) / 2];
		int index = 0;
		for (int i = 0; i < nrow; i++){
			for (int j = i + 1; j < nrow; j ++){
				double distSq = 0.;
				for (int k = 0; k < ncol; k++){
					distSq  += Math.pow(Y[i][k] - Y[j][k], 2.);									
				}
				distances[index] = Math.sqrt(distSq);
				index += 1;
			}
		}
		return distances;	
	}
	
	public static double[] UpTriangular(double[][] D){  // D is a square form, zero at the diagonal
		int nrow = D.length;   // nrow = ncol
		double[] upTriangular = new double[nrow * (nrow - 1) / 2];
		int index = 0;
		for (int i = 0; i < nrow ; i++){
			for (int j = i + 1; j < nrow; j++){
				upTriangular[index] = D[i][j];
				index += 1;
			}
		}
		return upTriangular;
		
	} 
	
	public static double max(double[] arg){
		int n = arg.length;
		double maxx = arg[0];
		for (int i = 1; i < n; i++){
			double temp = arg[i];
			if (temp > maxx){
				maxx = temp;
			}
		}
		return maxx;
	}
	
	public static double[][] squareForm(double[] arg){        // transpose arg to be a squared symmetric matrix
		int n = arg.length;
		int nrow = (int) Math.ceil(Math.sqrt(2.*n));  // nrow = ncol
		double[][] square = new double[nrow][nrow];
		
		int index = 0;
		for (int i = 0; i < nrow; i++){
			for (int j = 0; j < nrow; j++){
				if (j == i){
					square[i][j] = 0.;					
				} else if (j > i){
					square[i][j] = arg[index];
					index += 1;
				} else {
					square[i][j] = square[j][i];
				}			
			}
		}
		return square;
	}
	
	public static int allGreatThanZero(double [] args){   //return 0 or 1
		int n = args.length;
		int flag = 1;
		for (int i = 0; i < n; i++){
			if (args[i] <= 0){
				flag = 0;
				return flag;
			}
		}
		return flag;
	}
	
	public double getStress(){
		return stress_;
	}
	
	public double[][] getConfiguration(){
		return configuration;
	}
	
	public static void printArray(double[][] args){
		int nrow = args.length;
		int ncol = args[0].length;
		for(int i=0; i<nrow; i++){
			for(int j=0; j<ncol; j++){
				System.out.print(args[i][j] + "\t");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	
	public static void main(String[] args){		
		double[][] input1={{0,3313,2963,3175,3339,2762,3276,2610,4485,2977,3030,4532,2753,3949,2865,2282,2179,3000,817,3927,1991},
		{3313,0,1318,1326,1294,1498,2218,803,1172,2018,1490,1305,645,636,521,1014,1365,1033,1460,2868,1802},
		{2963,1318,0,204,583,206,966,677,2256,597,172,2084,690,1558,1011,925,747,285,1511,1616,1175},
		{3175,1326,204,0,460,409,1136,747,2224,714,330,2052,739,1550,1059,1077,977,280,1662,1786,1381},
		{3339,1294,583,460,0,785,1545,853,2047,1115,731,1827,789,1347,1101,1209,1160,340,1794,2196,1588},
		{2762,1498,206,409,785,0,760,1662,2436,460,269,2290,714,1764,1035,911,583,465,1497,1403,937},
		{3276,2218,966,1136,1545,760,0,1418,3196,460,269,2971,1458,2498,1778,1537,1104,1176,2050,650,1455},
		{2610,803,677,747,853,1662,1418,0,1975,1118,895,1936,158,1439,425,328,591,513,995,2068,1019},
		{4485,1172,2256,2224,2047,2436,3196,1975,0,2897,2428,676,1817,698,1693,2185,2565,1971,2631,3886,2974},
		{2977,2018,597,714,1115,460,460,1118,2897,0,550,2671,1159,2198,1479,1238,805,877,1751,949,1155},
		{3030,1490,172,330,731,269,269,895,2428,550,0,2280,863,1730,1183,1098,851,457,1683,1500,1205},
		{4532,1305,2084,2052,1827,2290,2971,1936,676,2671,2280,0,1178,668,1762,2250,2507,1799,2700,3231,2937},
		{2753,645,690,739,789,714,1458,158,1817,1159,863,1178,0,1281,320,328,724,471,1048,2108,1157},
		{3949,636,1558,1550,1347,1764,2498,1439,698,2198,1730,668,1281,0,1157,1724,2010,1273,2097,3188,2409},
		{2865,521,1011,1059,1101,1035,1778,425,1693,1479,1183,1762,320,1157,0,618,1109,792,1011,2428,1363},
		{2282,1014,925,1077,1209,911,1537,328,2185,1238,1098,2250,328,1724,618,0,331,856,586,2187,898},
		{2179,1365,747,977,1160,583,1104,591,2565,805,851,2507,724,2010,1109,331,0,821,946,1754,428},
		{3000,1033,285,280,340,465,1176,513,1971,877,457,1799,471,1273,792,856,821,0,1476,1827,1249},
		{817,1460,1511,1662,1794,1497,2050,995,2631,1751,1683,2700,1048,2097,1011,586,946,1476,0,2707,1209},
		{3927,2868,1616,1786,2196,1403,650,2068,3886,949,1500,3231,2108,3188,2428,2187,1754,1827,2707,0,2105},
		{1991,1802,1175,1381,1588,937,1455,1019,2974,1155,1205,2937,1157,2409,1363,898,428,1249,1209,2105,0}};
		
		double[][] input2 = {{0,31800,32949,41280,47989,48147,49051,49162,49250,49470,49491,49571,49696,49704,49723,49724,49896,49913,49936,49938,49948,49959,49961,49974,49978,49978,49979,49981,49982,49986,49988,49989,49989,49989},
				{31800,0,1149,9480,16189,16347,17251,17362,17450,17670,17691,17771,17896,17904,17923,17924,18096,18113,18136,18138,18148,18159,18161,18174,18178,18178,18179,18181,18182,18186,18188,18189,18189,18189},
				{32949,1149,0,8331,15040,15198,16102,16213,16301,16521,16542,16622,16747,16755,16774,16775,16947,16964,16987,16989,16999,17010,17012,17025,17029,17029,17030,17032,17033,17037,17039,17040,17040,17040},
				{41280,9480,8331,0,6709,6867,7771,7882,7970,8190,8211,8291,8416,8424,8443,8444,8616,8633,8656,8658,8668,8679,8681,8694,8698,8698,8699,8701,8702,8706,8708,8709,8709,8709},
				{47989,16189,15040,6709,0,158,1062,1173,1261,1481,1502,1582,1707,1715,1734,1735,1907,1924,1947,1949,1959,1970,1972,1985,1989,1989,1990,1992,1993,1997,1999,2000,2000,2000},
				{48147,16347,15198,6867,158,0,904,1015,1103,1323,1344,1424,1549,1557,1576,1577,1749,1766,1789,1791,1801,1812,1814,1827,1831,1831,1832,1834,1835,1839,1841,1842,1842,1842},
				{49051,17251,16102,7771,1062,904,0,111,199,419,440,520,645,653,672,673,845,862,885,887,897,908,910,923,927,927,928,930,931,935,937,938,938,938},
				{49162,17362,16213,7882,1173,1015,111,0,88,308,329,409,534,542,561,562,734,751,774,776,786,797,799,812,816,816,817,819,820,824,826,827,827,827},
				{49250,17450,16301,7970,1261,1103,199,88,0,220,241,321,446,454,473,474,646,663,686,688,698,709,711,724,728,728,729,731,732,736,738,739,739,739},
				{49470,17670,16521,8190,1481,1323,419,308,220,0,21,101,226,234,253,254,426,443,466,468,478,489,491,504,508,508,509,511,512,516,518,519,519,519},
				{49491,17691,16542,8211,1502,1344,440,329,241,21,0,80,205,213,232,233,405,422,445,447,457,468,470,483,487,487,488,490,491,495,497,498,498,498},
				{49571,17771,16622,8291,1582,1424,520,409,321,101,80,0,125,133,152,153,325,342,365,367,377,388,390,403,407,407,408,410,411,415,417,418,418,418},
				{49696,17896,16747,8416,1707,1549,645,534,446,226,205,125,0,8,27,28,200,217,240,242,252,263,265,278,282,282,283,285,286,290,292,293,293,293},
				{49704,17904,16755,8424,1715,1557,653,542,454,234,213,133,8,0,19,20,192,209,232,234,244,255,257,270,274,274,275,277,278,282,284,285,285,285},
				{49723,17923,16774,8443,1734,1576,672,561,473,253,232,152,27,19,0,1,173,190,213,215,225,236,238,251,255,255,256,258,259,263,265,266,266,266},
				{49724,17924,16775,8444,1735,1577,673,562,474,254,233,153,28,20,1,0,172,189,212,214,224,235,237,250,254,254,255,257,258,262,264,265,265,265},
				{49896,18096,16947,8616,1907,1749,845,734,646,426,405,325,200,192,173,172,0,17,40,42,52,63,65,78,82,82,83,85,86,90,92,93,93,93},
				{49913,18113,16964,8633,1924,1766,862,751,663,443,422,342,217,209,190,189,17,0,23,25,35,46,48,61,65,65,66,68,69,73,75,76,76,76},
				{49936,18136,16987,8656,1947,1789,885,774,686,466,445,365,240,232,213,212,40,23,0,2,12,23,25,38,42,42,43,45,46,50,52,53,53,53},
				{49938,18138,16989,8658,1949,1791,887,776,688,468,447,367,242,234,215,214,42,25,2,0,10,21,23,36,40,40,41,43,44,48,50,51,51,51},
				{49948,18148,16999,8668,1959,1801,897,786,698,478,457,377,252,244,225,224,52,35,12,10,0,11,13,26,30,30,31,33,34,38,40,41,41,41},
				{49959,18159,17010,8679,1970,1812,908,797,709,489,468,388,263,255,236,235,63,46,23,21,11,0,2,15,19,19,20,22,23,27,29,30,30,30},
				{49961,18161,17012,8681,1972,1814,910,799,711,491,470,390,265,257,238,237,65,48,25,23,13,2,0,13,17,17,18,20,21,25,27,28,28,28},
				{49974,18174,17025,8694,1985,1827,923,812,724,504,483,403,278,270,251,250,78,61,38,36,26,15,13,0,4,4,5,7,8,12,14,15,15,15},
				{49978,18178,17029,8698,1989,1831,927,816,728,508,487,407,282,274,255,254,82,65,42,40,30,19,17,4,0,0,1,3,4,8,10,11,11,11},
				{49978,18178,17029,8698,1989,1831,927,816,728,508,487,407,282,274,255,254,82,65,42,40,30,19,17,4,0,0,1,3,4,8,10,11,11,11},
				{49979,18179,17030,8699,1990,1832,928,817,729,509,488,408,283,275,256,255,83,66,43,41,31,20,18,5,1,1,0,2,3,7,9,10,10,10},
				{49981,18181,17032,8701,1992,1834,930,819,731,511,490,410,285,277,258,257,85,68,45,43,33,22,20,7,3,3,2,0,1,5,7,8,8,8},
				{49982,18182,17033,8702,1993,1835,931,820,732,512,491,411,286,278,259,258,86,69,46,44,34,23,21,8,4,4,3,1,0,4,6,7,7,7},
				{49986,18186,17037,8706,1997,1839,935,824,736,516,495,415,290,282,263,262,90,73,50,48,38,27,25,12,8,8,7,5,4,0,2,3,3,3},
				{49988,18188,17039,8708,1999,1841,937,826,738,518,497,417,292,284,265,264,92,75,52,50,40,29,27,14,10,10,9,7,6,2,0,1,1,1},
				{49989,18189,17040,8709,2000,1842,938,827,739,519,498,418,293,285,266,265,93,76,53,51,41,30,28,15,11,11,10,8,7,3,1,0,0,0},
				{49989,18189,17040,8709,2000,1842,938,827,739,519,498,418,293,285,266,265,93,76,53,51,41,30,28,15,11,11,10,8,7,3,1,0,0,0},
				{49989,18189,17040,8709,2000,1842,938,827,739,519,498,418,293,285,266,265,93,76,53,51,41,30,28,15,11,11,10,8,7,3,1,0,0,0}};
//		NonMetricMDS test = new NonMetricMDS(input2,2, true);  // try NonMetricMDS(input2, 2, false) or NonMetricMDS(input2)
		NonMetricMDS test = new NonMetricMDS(input2,2, true);
		System.out.println("======stress===========");
		System.out.println(test.getStress());
		System.out.println("======configuration===========");
		printArray(test.getConfiguration());
	}

}
