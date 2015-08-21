package algorithm.sort;

public class BubbleSort {

	public static void main(String[] args) {
		int[] array={2,3,1,7,4,9,6,3,1};
		bubbleSort(array);
		
	}
	/**
	 * 冒泡降序
	 * @param array
	 */
	public static void bubbleSort(int[] array){
		for(int i=0;i<array.length-1;i++){
			for(int j=0;j<array.length-1-i;j++){
				if(array[j]<array[j+1]){
					int temp=array[j];
					array[j]=array[j+1];
					array[j+1]=temp;
				}
			}
		}
	}
}
