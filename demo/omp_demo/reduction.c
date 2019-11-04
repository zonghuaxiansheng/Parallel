#include <omp.h>
#include <stdio.h>

#define size 1000000
double a[size], b[size], result;

int main ()  
{
	int   i, n, chunk;
	
	/* Some initializations */
	n = size;
	chunk = 10;
	result = 0.0;
	for (i=0; i < n; i++)  
	{  
		a[i] = i * 1.0; 
		b[i] = i * 2.0; 
 	}
	#pragma omp parallel for default(shared) reduction(+:result) private(i)
	for (i=0; i < n; i++)    
		result = result + (a[i] * b[i]);
	printf("Final result= %f\n",result);
} 

