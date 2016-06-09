// Algorithms: recursively divide the longer dimension of the screen using the median-cut 
// until the number of tiles equals the number of processors. 
// Hushiqi & QianWang
#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];

//loop indexs
int i = 0;
int j = 0;
int k = 0;
int z = 0;

int numprocs;  /* Number of processors to use */
int myid;
int num_quadrants;

//swap the array elemnt at index i and j
void swap(unsigned int array[], int i, int j);
//find the kth smallest element in array v, implemented using medians of medians whose runtime is O(n)
unsigned int find_kth(unsigned int *v, int n, int k, unsigned int *y);
//global_cost computed at process 0
double global_cost = 0;
//find border of initial quadrant by finding the border of X and Y
void find_quadrants(int num_quadrants, int left, int right, int cur_quadrants,int *coordinates,int *coordinate_pointer,int *xory) {  
    if (num_quadrants <= cur_quadrants) {
        return;
    }
    printf("\nComputing in find_quadrants.%d \n",cur_quadrants);
    int min_x = X_axis[left];
    int max_x = X_axis[left];
    int min_y = Y_axis[left];
    int max_y = Y_axis[left];
    for (i = left; i < right + 1; i++) {
        if (X_axis[i] < min_x) {
            min_x = X_axis[i];
        }
        if (X_axis[i] > max_x) {
            max_x = X_axis[i];
        }
        if (Y_axis[i] < min_y) {
            min_y = Y_axis[i];
        }
        if (Y_axis[i] > max_y) {
            max_y = Y_axis[i];
        }
    }
    int x_dimension;
    int y_dimension;
    x_dimension = max_x - min_x;
    y_dimension = max_y - min_y;
    //check the length of each dimension, choose the longest one to cut
    if(x_dimension >  y_dimension){
        printf("\ncut in x.\n");
        printf("\nleft:  %d \n", left);
        printf("\nright:  %d \n", right);
        //find the median of the corresponding part of X_axis[]
        int x_pivot = find_kth(X_axis+left, right-left + 1,(right - left+1)/2 -1, Y_axis+left);
        printf("\nxpivot:  %d \n", x_pivot);
        k = left;
        j = right;
        //make sure that the numbers smaller than the median are all gathered together on the left half of the array
        while (k <= j && k < left + (right - left)/2) {
            if (X_axis[k] > x_pivot) {
                while (X_axis[j] > x_pivot) {
                    j--;
                }
                if (k < j) {
                    swap(X_axis, k, j);
                    swap(Y_axis, k, j);
                } 
            }
            k++;
        }
        printf("\nk:  %d \n", k);
        printf("\nafter cut x   \n");
        //record the point of seperation for future use 
        int index = *coordinate_pointer;
        printf("\n x_index:  %d \n",*coordinate_pointer); 
        (*coordinate_pointer)++;
        //UPDATE Coordinates[][]
        *(coordinates+index) = x_pivot;
        *(*(coordinates+index)+1) =min_x;
        *(*(coordinates+index)+2) =min_y;
        *(*(coordinates+index)+3) =max_x;
        *(*(coordinates+index)+4) =max_y;

        *(xory+index) = 0;
        if (num_quadrants <= 2 * cur_quadrants)
        {
            printf(" (%d,%d) ", min_x, max_y);
            printf(" (%d,%d) ", x_pivot, max_y);
            printf(" (%d,%d) ", x_pivot,  min_y);
            printf(" (%d,%d) \n", min_x, min_y);

            printf(" (%d,%d) ", x_pivot, max_y);
            printf(" (%d,%d) ", max_x, max_y);
            printf(" (%d,%d) ", max_x,  min_y);
            printf(" (%d,%d) \n",x_pivot,  min_y);
        }
        find_quadrants (num_quadrants,left,left + (right - left)/2-1,cur_quadrants*2,&coordinates[0],&coordinate_pointer[0],&xory[0]);
        find_quadrants (num_quadrants,left + (right - left)/2,right,cur_quadrants*2,&coordinates[0],&coordinate_pointer[0],&xory[0]);
    }
    //if Y dimension is longer, cut on Y dimension, Algorithms SAME as cut on X's dimension
    else {
        printf("\ncut in y.\n");
        int y_pivot = find_kth(Y_axis+left, right-left + 1, (right - left+1)/2 -1, X_axis+left);
        printf("\nypivot:  %d \n", y_pivot);
        k = left;
        j = right;
        while (k <= j && k < left + (right - left)/2) {
            if (Y_axis[k] > y_pivot) {
                while (Y_axis[j] > y_pivot) {
                    j--;
                }
                if (k < j) {
                    swap(X_axis, k, j);
                    swap(Y_axis, k, j);
                } 
            }
            k++;
        }
        printf("\n cut after y.\n");
        int index = *coordinate_pointer;
        printf("\n y_index:  %d \n", *coordinate_pointer); 
        (*coordinate_pointer)++;
        //UPDATE Coordinates[][]
        *(coordinates+index) = y_pivot;
        *(*(coordinates+index)+1) =min_x;
        *(*(coordinates+index)+2) =min_y;
        *(*(coordinates+index)+3) =max_x;
        *(*(coordinates+index)+4) =max_y;

        *(xory+index) = 1;
        if (num_quadrants <= 2 * cur_quadrants)
        {
            printf(" (%d,%d) ", min_x, max_y);
            printf(" (%d,%d) ", max_x, max_y);
            printf(" (%d,%d) ", max_x, y_pivot);
            printf(" (%d,%d) \n", min_x,y_pivot);

            printf(" (%d,%d) ", min_x,y_pivot);
            printf(" (%d,%d) ", max_x, y_pivot);
            printf(" (%d,%d) ", max_x,  min_y);
            printf(" (%d,%d) \n",min_x,  min_y);
        }
        find_quadrants (num_quadrants,left,(left+right)/2-1,cur_quadrants*2,&coordinates[0],&coordinate_pointer[0],&xory[0]);
        find_quadrants (num_quadrants,(left+right)/2,right,cur_quadrants*2,&coordinates[0],&coordinate_pointer[0],&xory[0]);
    }
    return;
}

void find_quad(int num_quadrants,int *coordinates) 
{
    int xory[num_quadrants];
    //int coordinates[num_quadrants][5];
    int coordinate_pointer[1];
    coordinate_pointer[0] = 0;
    find_quadrants(num_quadrants,0,NUM_POINTS,1,&coordinates[0],&coordinate_pointer[0],&xory[0]);
    return;
}
// find the kth smallest key in a, n is the number of elements in array a
// Use medians of medians Algorithms which claims to have linear worst-case runtime for selecting the k th largest
// https://www.quora.com/What-is-the-most-efficient-algorithm-to-find-the-kth-smallest-element-in-an-array-having-n-unordered-elements
unsigned int find_kth(unsigned int *a, int n, int k, unsigned int *y) 
{
    if (n == 1 && k == 0) return a[0];
    int j0 = 0;
    int i1 = 0;
    int j1 = 0;
    //divide the array into n/5 subarrays of 5 elements and find the medians for each subarrays
    int number = (n + 4)/5; // the number of subarrays 
    unsigned int *medians =  (unsigned int *)malloc(number* sizeof(int)); 
    for (i1=0; i1<num; i1++) {
        if (5*i1 + 4 < n) { // sort the subarray 
            unsigned int *w1 = y + 5*i1;
            unsigned int *w = a + 5*i1;
            for (j0=0; j0<3; j0++) {
                int jmin = j0;
                for (j1=j0+1; j1<5; j1++) {
                    if (w[j1] < w[jmin]) jmin = j1;
                }
                swap(w, j0, jmin);
                swap(w1, j0, jmin);
            }
            medians[i1] = w[2];
        } else {
            medians[i1] = a[5*i1];
        }
    }
    //find the median of medians
    int pivot = find_kth(medians, number, number/2, medians);
    free(medians);
    for (i1=0; i1<n; i1++) {
        if (a[i1] == pivot) {
            swap(a, i1, n-1);
            swap(y, i1, n-1);
            break;
        }
    }
 // keep the number of elements that are smaller than the median pivot
    int num = 0;
    for (i1=0; i1<n-1; i1++) {
        if (a[i1] < pivot) {
            swap(a, i1, num);
            swap(y, i1, num);
            num++;
        }
    }
    swap(a, num, n-1);
    swap(y, num, n-1);
    if (num == k) {
       //return when exactly k elemens are smaller than pivot
        return pivot;
    } else if (num > k) {
      //find kth smallest in the left num amount of elements in a
        return find_kth(a, num, k, y);
    } else {
      //find the k-num-1 smallest in the right side of the a staring at index a+store+1
        return find_kth(a+num+1, n-num-1, k-num-1, y+num+1);
    }
}

// swap position x and j of array 
void swap(unsigned int array[], int i, int j) {
    int temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}

int main(argc,argv)
int argc;
char *argv[];
{
    int num_quadrants;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int coordinates[num_quadrants][5];

    /*Time Variables*/
    double startwtime = 0.0, endwtime;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    
    if (argc != 2)
    {
        fprintf (stderr, "Usage: recursive_bisection <#of quadrants>\n");
        MPI_Finalize();
        exit (0);
    }
    
    fprintf (stderr,"Process %d on %s\n", myid, processor_name);
    
    num_quadrants = atoi (argv[1]);
    
    if (myid == 0) {
        fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);
        int i;
        srand (10000);
        for (i = 0; i < NUM_POINTS; i++)
            X_axis[i] = (unsigned int)rand();
        for (i = 0; i < NUM_POINTS; i++)
            Y_axis[i] = (unsigned int)rand();
        //start timer at process 0
        printf("\nComputing Parallely Using MPI.\n");
        startwtime = MPI_Wtime();
    }
    if (myid == 0) find_quad(num_quadrants,&coordinates);
    MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
   //Parallel Commputing-->Compute the Cost and finally reduced to Processor 0 
    double local_cost = 0;
    for (i = myid; i < num_quadrants; i += numprocs)
    {
       int points = NUM_POINTS / num_quadrants;
       for (j = 0; j < points - 1; j++) {
            for (k = j+1; k < points; k++) {
                int x1 = points * i + j;
                int x2 = points * i + k;
                int y1 = points * i + j;
                int y2 = points * i + k;

                double diff_x = abs(X_axis[x1] - X_axis[x2]);
                double diff_y = abs(Y_axis[y1] - Y_axis[y2]);
                // add the cost to the local_cost
                local_cost += sqrt((double)diff_x * diff_x + diff_y * diff_y);
            }
       }     
    }
    //Using Reduce the calculate the total cost at  Processor 0.
    MPI_Reduce(&local_cost, &global_cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("\nTotal cost:  %lf \n", global_cost);

    if (myid == 0) {
       printf("\nFInalTotal cost: %lf \n", global_cost);
      //end timer at process 0
        endwtime = MPI_Wtime();
      //the total execution time
        printf("\nelapsed time = %f\n", endwtime - startwtime);
       //print the total partition cost
    }
    MPI_Finalize();
    return 0;
}


