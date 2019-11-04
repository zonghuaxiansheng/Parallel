#include <omp.h> 
#include <stdio.h>
static long num_steps = 8;
double step; 
#define NUM_THREADS 8 
void main () 
{	  int i; 	  
	  long long sum = 1.0; 
	  omp_set_num_threads(NUM_THREADS); 
	  #pragma omp parallel for reduction(*:sum) private(i) 
	  for (i=1;i<=num_steps; i++){
		  int id = omp_get_thread_num();
		  printf("I am therad %d : i=%d , sum=%lld\n", id, i, sum); 
		  sum = sum - 1 + i; 
	  }  
	  printf("sum = %lld\n",sum);
}
