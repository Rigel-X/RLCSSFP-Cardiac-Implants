#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //prhs: Array of pointers: contains all the start adresses of input arguments
    //Printf("hello world"); //prints hello world
    //printf(mxArrayToString(prhs[0])); //prints first argument as a string
    //printf(mxGetString(prhs[0]));
    
    double *N;  //number of points for ram-lak (input)
    double *outptr; //output array of ramlak coeffs
    int count;   //count variable from 0 to N-1
    int index;   //real index from -N/2 to N/2
    int max;    //later uses variable (for normalizing)
    double pi=3.14159265;
    //double Nr;
    //double* points;
    
    ////assign pointers
    ///////////////////
    N = mxGetPr(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(1,*N,mxREAL); //allocate memory for output array
    outptr = mxGetPr(plhs[0]); //now outptr is the output-variable
    //*outptr=*N;
    
    
    //// Compute RamLak coeffs
    //////////////////////////
    count=0;
    index=-*N/2;
    while(count<*N)
    {
        if(index==0)
            *outptr=pi/(2*pow(1,2));  //dt=1
        else if((index%2)==1 || (index%2)==-1)    //if odd
            *outptr=-2/(pi*pow(index*1,2));   //dt=1
        else
            *outptr=0;
        
        outptr++;
        count++;
        index++;
    }
    
    
//     ////normalize
//     //////////////////
//     
//     //find maximum
//     count=0;
//     max=0;
//     outptr = outptr-(int)*N; //set outptr to initial adress
//     while(count<*N)
//     {
//         if(*outptr>max)
//             max=*outptr;
//         count++;
//         outptr++;
//     }
//     *outptr = max;
//     //divide
//     count=0;
//     outptr = outptr-(int)*N; //set outptr to initial adress
//     while(count<*N)
//     {
//         *outptr=*outptr/max;
//         count++;
//         outptr++;
//     }
        
        
        
//     
//     if(mod(*N/2,2)==0)
//         Nr=*N+2;
//     else
//         Nr=*N;
//     
//     for(count=0; count<Nr; count++)
//     {
//         if(count==Nr/2)
//             output[count]=pi/(2*pow(1,2));  //dt=1
//         else if(mod(count,2)==0)
//             output[count]=-2/(pi*pow((count-*N/2)*1,2));   //dt=1
//         else
//             output[count]=0;
//     }
//     
    
}