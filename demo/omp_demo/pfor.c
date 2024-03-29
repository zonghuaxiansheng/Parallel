#include <omp.h>
#define N       100000
#define CHUNKSIZE   1000
double  a[N], b[N], c[N];

int main () {
	int i, chunk;
	
	/* Some initializations */
	for (i=0; i < N; i++)  
		a[i] = b[i] = i * 1.0;
	chunk = CHUNKSIZE;
	#pragma omp parallel for  schedule(static,CHUNKSIZE) private(i)
		for (i=0; i < N; i++)    
			{
			 int id; 

			 c[i] = a[i] + b[i];

			 id = omp_get_thread_num();
			 if (  (i % chunk) == 0 )  
			 	 printf("Iteration #%d in thread #%d\n",i, id);

			}
} 
