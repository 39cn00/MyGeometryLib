#include <stdio.h>
#include <time.h>
#include "lageometry.h"

int main(){
    printf("ELO123\n");

    lag::Dmatrixf a;
    a.newMatrix(100, 100);
    a.loadTestValues();
    
    //a.setValue(0, 0, 4.0f);
    //a.setValue(1, 0, 5.0f);
    //a.setValue(2, 0, 7.0f);
    //a.setValue(0, 1, 2.0f);
    //a.setValue(1, 1, 1.0f);
    //a.setValue(2, 1, 0.0f);

    lag::Dmatrixf b;
    b.newMatrix(100, 100);
    b.loadTestValues();
    //b.setValue(0, 0, 2.0f);
    //b.setValue(1, 0, 3.0f);
    //b.setValue(0, 1, 8.0f);
    //b.setValue(1, 1, 9.0f);
    //b.setValue(0, 2, 1.0f);
    //b.setValue(1, 2, 1.0f);

    lag::Dmatrixf c;
    //c.newMatrix(2, 2);

    clock_t tStart = clock();
    { // Timed section
        lag::multiplyMx(a, b, c);
    }
    double tEnd = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    tEnd *= 1000.0f; //ms
    printf("Exe time = %0.3fms\n", tEnd);

    //lag::prontMatrix(a);
    //lag::prontMatrix(b);
    //lag::prontMatrix(c);
    return 0;
}
