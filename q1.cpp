#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

using namespace std;

class IntegerSorter{
private:

	int* Keys;
	int TotalKeys;


	void BubbleSort(){
		// Your implementation of Bubble Sort such that it also counts total comparisons, total swaps and total time taken
		// to sort the 100,000 numbers
		double t1,t2;
		t1=clock();
		int temp;
		long long int cmp=0;
		long long int swp=0;
		for(long i = 0; i < TotalKeys; i++)
		{
			for(long j = 0; j <TotalKeys-i-1; j++)
			{
				cmp++;
				if (Keys[j] > Keys[j+1])
				{
					temp        = Keys[j];
					Keys[j]     = Keys[j+1];
					Keys[j+1]   = temp;
					swp++;
					//cmp++;
				}
			}
		}
		t2 = clock();
		TotalOperationsAndTime[0]=swp;
		TotalOperationsAndTime[1]=cmp;
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


	}
	void SelectionSort(){
		// Your implementation of Selection Sort such that it also counts total comparisons, total //swaps and total time taken
		// to sort the 100,000 numbers
		double t1,t2;
		t1=clock();

		rapselectionSort(Keys,TotalKeys);
		t2 = clock();
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


	}
	void rapselectionSort(int arr[], int n)  
	{  
		long i, j, min_idx;
		long long int cmp,swp;  

		// One by one move boundary of unsorted subarray  
		for (i = 0; i < n-1; i++)  
		{  
			// Find the minimum element in unsorted array  
			min_idx = i;  
			for (j = i+1; j < n; j++){
				cmp++;
				if (arr[j] < arr[min_idx]){  
					min_idx = j;
					swap(&arr[min_idx], &arr[i]);  
					swp++;
				}
			}
		}

		TotalOperationsAndTime[0]=swp;
		TotalOperationsAndTime[1]=cmp;

	}
	void swap(int *xp, int *yp)  
	{  
		int temp = *xp;  
		*xp = *yp;  
		*yp = temp;  
	}  



	void InsertionSort(){
		// Your implementation of Insertion Sort such that it also counts total comparisons,
		// total swaps and total time taken
		// to sort the 100,000 numbers
		double t1,t2;
		t1= clock();
		int temp;

		long long int cmp=0;
		long long int swp=0;
		for(long i = 1; i < TotalKeys; i++)
		{
			temp = Keys[i];
			long j;
			for(j = i-1; j >= 0 && Keys[j] > temp; j--)
			{
				Keys[j+1] = Keys[j];
				cmp++;
			}
			Keys[j+1] = temp;
			swp++;
		}
		t2 = clock();
		TotalOperationsAndTime[0]=swp;
		TotalOperationsAndTime[1]=cmp;
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


	}


	void heapify(int arr[], int n, int i) {
		int temp;
		long long int swp,cmp;
		int largest = i;
		int l = 2 * i + 1;
		int r = 2 * i + 2;
		if (l < n && arr[l] > arr[largest])
			largest = l,cmp++;
		if (r < n && arr[r] > arr[largest])
			largest = r,cmp++;
		if (largest != i) {
			temp = arr[i];
			arr[i] = arr[largest];
			arr[largest] = temp;
			swp++;


			heapify(arr, n, largest);
		}
	}
	void heapSortrap(int arr[], int n) {
		int temp;
		for (int i = n / 2 - 1; i >= 0; i--)
			heapify(arr, n, i);
		for (int i = n - 1; i >= 0; i--) {
			temp = arr[0];
			arr[0] = arr[i];
			arr[i] = temp;
			heapify(arr, i, 0);
		}
	}

	void HeapSort(){
		// Your implementation of Heap Sort such that it also counts total comparisons,
		// total swaps and total time taken
		// to sort the 100,000 numbers
		double t1,t2;
		t1 = clock();

		heapSortrap(Keys,TotalKeys);
		t2 = clock();
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


	}
	// main function to do heap sort


	void Merge(int l, int m, int r){
		// To merge the two parts of keys from left - Mid and Mid+1 to Right
		long long int cmp=0,swp=0;

		int i, j, k;
		int n1 = m - l + 1;
		int n2 =  r - m;

		/* create temp arrays */
		int *L = new int [n1];
		int *R = new int [n2];

		/* Copy data to temp arrays L[] and R[] */
		for (i = 0; i < n1; i++)
			L[i] = Keys[l + i],cmp++;
		for (j = 0; j < n2; j++)
			R[j] = Keys[m + 1+ j],cmp++;

		/* Merge the temp arrays back into arr[l..r]*/
		i = 0; // Initial index of first subarray
		j = 0; // Initial index of second subarray
		k = l; // Initial index of merged subarray
		while (i < n1 && j < n2)
		{
			if (L[i] <= R[j])
			{
				Keys[k] = L[i];
				i++;
				cmp++;
			}
			else
			{
				Keys[k] = R[j];
				j++;
				cmp++;
			}
			k++;
		}

		/* Copy the remaining elements of L[], if there
		are any */
		while (i < n1)
		{
			Keys[k] = L[i];
			i++;
			k++;
			cmp++;
		}

		/* Copy the remaining elements of R[], if there
		are any */
		while (j < n2)
		{
			Keys[k] = R[j];
			j++;
			k++;
			cmp++;
		}
	}

	void TopDownMergeSort(){
		double t1,t2;
		t1 = clock();
		rapTopDownMergeSort(Keys,0,100000);
		t2 = clock();
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;

	}



	void rapTopDownMergeSort(int *keys,int l,int r){

		// Your implementation of top-down (i.e. recursive) Merge Sort
		// such that it also counts total comparisons, total swaps and total time taken
		// to sort the 100,000 numbers

		if (l < r)

		{
			// Same as (l+r)/2, but avoids overflow for
			// large l and h
			int m = l+(r-l)/2;

			// Sort first and second halves
			rapTopDownMergeSort(Keys, l, m);
			rapTopDownMergeSort(Keys, m+1, r);

			Merge(l, m, r);
		}

	}

	void BottomUpMergeSort(){
		// Your implementation of bottom-up (i.e. iterative) Merge Sort
		// such that it also counts total comparisons, total swaps and total time taken
		// to sort the 100,000 numbers
		double t1,t2;
		t1=clock();

		int curr_size; 
		int left_start;
		for (curr_size=1; curr_size<=TotalKeys-1; curr_size = 2*curr_size) 
		{ 
			for (left_start=0; left_start<TotalKeys-1; left_start += 2*curr_size) 
			{ 
				int mid = min(left_start + curr_size - 1, TotalKeys-1); 

				int right_end = min(left_start + 2*curr_size - 1, TotalKeys-1); 

				Merge(left_start, mid, right_end); 
			} 
		
		}
		t2 = clock();
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;
	}

	void QuickSort(){

		double t1,t2;
		t1 = clock();
		qicksortrap(0,TotalKeys-1);
		t2 = clock();

		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


		// Your implementation of quick sort
		// Also count the number of comparisons and swaps operations
		// along with the time used for sorting

	}

	void qicksortrap(int left,int right){
		if (left < right)
		{
			long pivot = partition(left, right);
			qicksortrap(left, pivot-1);
			qicksortrap(pivot+1, right);
		}
	}



	long partition(long left, long right)
	{
		int pivot_element = Keys[left];
		int lb = left, ub = right;
		int temp;
		int cmp=0,swp=0;
		while (left < right)
		{
			while(Keys[left] <= pivot_element){
				left++;
				cmp++;
			}
			while(Keys[right] > pivot_element)
			{ right--;
			cmp++;
			}
			if (left < right)
			{
				temp        = Keys[left];
				Keys[left]  = Keys[right];
				Keys[right] = temp;
				swp++;
			}
			cmp++;
		}
		Keys[lb] = Keys[right];
		Keys[right] = pivot_element;
		return right;
	}



	void CountingSort(){
		// Your implementation of counting sort for counting integers in the range 0-MaxValue
		// Also compute number of seconds used for sorting
		double t1,t2;
		t1 = clock();
		rapcount(Keys,TotalKeys);
		t2 = clock();
		TotalOperationsAndTime[0]=0;
		TotalOperationsAndTime[1]=0;
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;


	}
	void rapcount(int *array, int size) {
		int output[100000+1];
		int max = getMax(array, size);
		int count[100000+1];     //create count array (max+1 number of elements)
		for(int i = 0; i<=max; i++)
			count[i] = 0;     //initialize count array to all zero
		for(int i = 1; i <=size; i++)
			count[array[i]]++;     //increase number count in count array.
		for(int i = 1; i<=max; i++)
			count[i] += count[i-1];     //find cumulative frequency
		for(int i = size; i>=1; i--) {
			output[count[array[i]]] = array[i];
			count[array[i]] -= 1; //decrease count for same numbers
		}
		for(int i = 1; i<=size; i++) {
			array[i] = output[i]; //store output array to main array
		}
	}

	int getMax(int array[], int size) {
		int max = array[1];
		for(int i = 2; i<=size; i++) {
			if(array[i] > max)
				max = array[i];
		}
		return max; //the max element from the arr
	}


	void RadixSort(){
		// Your implementation of counting sort for counting integers in the range 0-MaxValue
		// Also compute number of seconds used for sorting
		double t1,t2;
		t1 = clock();
		rapradixsort(Keys,TotalKeys);
		t2 = clock();

		TotalOperationsAndTime[0]=0;
		TotalOperationsAndTime[1]=0;
		TotalOperationsAndTime[2]= (t2-t1)/CLK_TCK ;



	}
	void rapradixsort(int arr[], int n)
	{
		int exp, m;
		m = getMax(arr, n);

		// Calling countSort() for digit at (exp)th place in every input.
		for (exp = 1; m/exp > 0; exp *= 100000)
			countSort(arr, n, exp);
	}

	// Count sort of arr[].
	void countSort(int arr[], int n, int exp)
	{
		// Count[i] array will be counting the number of array values having that 'i' digit at their (exp)th place.  
		int output[100000], i, count[10] = {0};

		// Count the number of times each digit occurred at (exp)th place in every input.
		for (i = 0; i < n; i++)
			count[(arr[i] / exp) % 10]++;

		// Calculating their cumulative count.
		for (i = 1; i < 10; i++)
			count[i] += count[i-1];

		// Inserting values according to the digit '(arr[i] / exp) % 10' fetched into count[(arr[i] / exp) % 10].
		for (i = n - 1; i >= 0; i--)
		{
			output[count[(arr[i] / exp) % 10] - 1] = arr[i];
			count[(arr[i] / exp) % 10]--;
		}

		// Assigning the result to the arr pointer of main().
		for (i = 0; i < n; i++)
			arr[i] = output[i];
	}



public:
	long long int TotalOperationsAndTime[3];

	IntegerSorter(){
		Keys = 0;
		TotalKeys = 0;
	}

	~IntegerSorter(){
		if(Keys)
			delete[] Keys;
	}

	bool LoadData(char FileName[], int TotalKeys){
		if (Keys)
			delete[] Keys;

		this->TotalKeys = TotalKeys;
		Keys = new  int [TotalKeys];

		if(!Keys)
			return false;

		ifstream fin(FileName);


		int  i=0;
		//int a;
		while(!fin.eof())
		{
			fin>>Keys[i];
			//cout<<Keys[i]<<" "<<i<<endl;

			i++;
		}
		// Open File and Load Data into the Array


		return true;
	}

	bool Sort(char Method){
		if(!Keys)
			return false;

		// Start Timer Here
		switch (Method){
		case 'B':
		case 'b':
			BubbleSort();
			break;
			case 'S':
			case 's':
			SelectionSort();
			break;
			case 'I':
			case 'i':
			InsertionSort();
			break;
			case 'm':
			TopDownMergeSort();
			break;
			case 'M':
			BottomUpMergeSort();
			break;

			case 'H':
			case 'h':
			HeapSort();
			break;
			case 'Q':
			case 'q':
			QuickSort();
			break;
			case 'C':
			case 'c':
			CountingSort();
			break;
			case 'R':
			case 'r':
			RadixSort();
			break;
		default:
			cout<<endl<<"Invalid Algorithm Choice"<<endl;
			return false;
			break;
		}

		return true;/*
					End Timer Here*/
	}
};


int main()
{
	IntegerSorter IS;
	//IS.LoadData("IntegerArray.txt",100000);
	bool Continue = true;

	do{
		cout<<"Select Algorithm "<<endl<<
			"B. Bubble"<<endl<<
			"S. Selection"<<endl<<
			"I. Insertion"<<endl<<

			"H. Heap"<<endl<<
			"m. Top-Down Merge"<<endl<<
			"M. Bottom-Up Merge"<<endl<<
			"Q. Quick"<<endl<<
			"C. Counting"<<endl<<
			"R. Radix"<<endl<<
			"E. Exit"<<endl;

		char Choice;
		cin>>Choice;

		if(Choice == 'E' || Choice == 'e')
			break;

		IS.TotalOperationsAndTime[0] = IS.TotalOperationsAndTime[1] = IS.TotalOperationsAndTime[2] = 0;
		char FileName[] = "IntegerArray.txt";


		if(!IS.LoadData(FileName, 100000))
			return -1;

		Continue = IS.Sort(Choice);

		if(Continue){
			cout << "Total Time Used By Sorting Procedure is\t"<< IS.TotalOperationsAndTime[2]<<endl
				<< "Total Operations Used By Sorting Procedure are\t" << IS.TotalOperationsAndTime[0]
			<<"\t"<<IS.TotalOperationsAndTime[1]<<endl;
		}
	}while(1);

	return 0;
}
