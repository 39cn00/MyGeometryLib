#include <stdio.h>
#include <time.h>
#include "lageometry.h"

#define PTHREAD

#ifdef PTHREAD
#define NAME "pThread"
#include <pthread.h>
#else
#define NAME "Normal"
#endif

lag::Dmatrixf a;
lag::Dmatrixf b;
lag::Dmatrixf c;

void * threadFunction(void * arg){
	MxChunk * chunk = (MxChunk *)arg;
	lag::multiplyMxPartial(a, b, c *chunk);
}

void runTest(){
	const int noofch = 4;

#ifdef PTHREAD

	pthread_t thread1, thread2;	
#else
	lag::multiplyMx(a, b, c);
#endif
}

int main(){
    printf("%s\n", NAME);

    a.newMatrix(100, 100);
    a.loadTestValues();
    
    //a.setValue(0, 0, 4.0f);
    //a.setValue(1, 0, 5.0f);
    //a.setValue(2, 0, 7.0f);
    //a.setValue(0, 1, 2.0f);
    //a.setValue(1, 1, 1.0f);
    //a.setValue(2, 1, 0.0f);

    b.newMatrix(100, 100);
    b.loadTestValues();
    //b.setValue(0, 0, 2.0f);
    //b.setValue(1, 0, 3.0f);
    //b.setValue(0, 1, 8.0f);
    //b.setValue(1, 1, 9.0f);
    //b.setValue(0, 2, 1.0f);
    //b.setValue(1, 2, 1.0f);

    //c.newMatrix(2, 2);

    clock_t tStart = clock();
    { // Timed section
        //lag::multiplyMx(a, b, c);
    	runTest();
    }
    double tEnd = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    tEnd *= 1000.0f; //ms
    printf("Exe time = %0.3fms\n", tEnd);

    //lag::prontMatrix(a);
    //lag::prontMatrix(b);
    //lag::prontMatrix(c);
    return 0;
}
