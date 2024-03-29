#include <omp.h>
#include <stdio.h>

/*

  for(i=4;i<100;i++){               
S1:       a[i] = b[i-2] + 1;
S2:       c[i] = b[i-1] + f[i];  
S3:       b[i] = a[i-1] + 2;    
S4:       d[i] = d[i+1] + b[i-1];     
  }
  
*/

#define Iter 10000
int a[Iter],b[Iter],c[Iter],d[Iter],f[Iter];
int a1[Iter],b1[Iter],c1[Iter],d1[Iter],f1[Iter];

int main()
{
  int i;
  

  int old_d[Iter]; // duplicating array-d to avoid
                   // anti-dependency  
  
  for(i=0;i<Iter;i++) a[i]=b[i]=c[i]=d[i]=f[i]=old_d[i]=i;

  for(i=0;i<Iter;i++) a1[i]=b1[i]=c1[i]=d1[i]=f1[i]=i;
  

  for(i=4;i<Iter;i++){
     a1[i] = b1[i-2] + 1;
     c1[i] = b1[i-1] + f1[i];
     b1[i] = a1[i-1] + 2;
     d1[i] = d1[i+1] + b1[i-1];
  }

/* The sequencial Loop below containing Statement S1 and S3 */

  for(i=4;i<Iter;i++){ 
     a[i] = b[i-2] + 1; // S1
     b[i] = a[i-1] + 2; // S3
  }

/* The Parallel Loop below only containing both Statement S2 and S4*/ 

#pragma omp parallel for shared(c,b,f,d,old_d) private(i)
  for(i=4;i<Iter;i++)
  {
    c[i] = b[i-1] + f[i] ; // S2
    d[i] = old_d[i+1] + b[i-1] ;  // S4 
  }  
    
/* The Parallel Loop below only containing Statement S4 */ 
    
//#pragma omp parallel for shared(d,b) private(i)
//  for(i=4;i<Iter;i++)
//      d[i] = old_d[i+1] + b[i-1] ; // S4    
  
for(i=4;i<Iter;i++) {

  if ( a[i]!=a1[i]) printf("a[%d] = %d , a1[%d] = %d\n",i,a[i],i,a1[i]);
  if ( b[i]!=b1[i]) printf("b[%d] = %d , b1[%d] = %d\n",i,b[i],i,b1[i]); 
  if ( c[i]!=c1[i]) printf("c[%d] = %d , c1[%d] = %d\n",i,c[i],i,c1[i]);
  if ( d[i]!=d1[i]) printf("d[%d] = %d , d1[%d] = %d\n",i,d[i],i,d1[i]);

}  
  
  
  
 return 0;  
}

