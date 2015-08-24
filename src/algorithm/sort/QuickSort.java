package algorithm.sort;

public class QuickSort {

	public static void main(String[] args) {
		int[] s={36,1,2,7,3,1,23,534,65};
		quickSort(s, 0, s.length-1);
		
	}
	public static void quickSort(int[] s,int l,int r){
		 if (l < r)  {  
		        int i = l, j = r, x = s[l];  
		        while (i < j)  {  
		            while(i < j && s[j] >= x) // 从右向左找第一个小于x的数  
		                j--;    
		            if(i < j)   
		                s[i++] = s[j];  
		            while(i < j && s[i] < x) // 从左向右找第一个大于等于x的数  
		                i++;    
		            if(i < j)   
		                s[j--] = s[i];  
		        }  
		        s[i] = x;  
		        quickSort(s, l, i - 1); // 递归调用   
		        quickSort(s, i + 1, r);  
		    } 
	  }
}
