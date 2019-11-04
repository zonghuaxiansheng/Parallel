#include <omp.h>
#include <stdio.h>
main()
{
	int x;
	x = 0;
	#pragma omp parallel shared(x) 
	    {
	      #pragma omp master
		{
			printf("num of threads : %d\n", omp_get_num_threads());      
		   	x = x + 10;
		}
	      #pragma omp critical
	                x = x + 1;	   	
	}  /* end of parallel section */
		
	printf("out of the parallel region : X = %d\n",x);		
} 


