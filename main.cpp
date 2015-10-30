/* 
 * File:   main.cpp
 * Author: jnortiz
 *
 * Created on March 10, 2015, 2:12 PM
 */

#include <cstdlib>
#include <math.h>
#include <NTL/RR.h>
#include "HIBE.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

int main(void) {
        
    RR sigmaRR = to_RR(3.195); // Standard deviation
    RR c = to_RR(0); // Center of the distribution            
    int security_level;
    float tailcut;
        
    int security_level_id = 3;
    
    // Parameters from (Saarinen, 2015)
    switch(security_level_id) {
        case 1: {
            security_level = 100;
            tailcut = 8.1;
            break;
        }
        case 2: {
            security_level = 128;
            tailcut = 9.2;
            break;
        }
        case 3: {
            security_level = 192;
            tailcut = 11.4;
            break;
        }
        case 4: {
            security_level = 256;
            tailcut = 13.2;
            break;
        }
        default:
            break;
            
    }//end-switch
    
    
    long precision = security_level/2;
    int nRectangles = 63; // Parameter of Ziggurat algorithm
    int omega = precision; // Parameter of Ziggurat algorithm
    
    Samplers *samplers = new Samplers();
    
    timestamp_t averageZiggurat, averageKnuthYao;
    timestamp_t ts_start, ts_end;

    averageZiggurat = 0.0;
    averageKnuthYao = 0.0;

    Vec<int> ZigguratPoly, KnuthPoly;
    int nSamples = 1024; // #coefficients in the polynomial

    int nIterations = 1000;
    for(int i = 0; i < nIterations; i++) {

        cout << endl;

        ts_start = get_timestamp();            
        ZigguratPoly = samplers->PolyGeneratorZiggurat(nSamples, nRectangles, sigmaRR, omega, precision, tailcut); // Coefficients, rectangles, sigma, omega and precision
        ts_end = get_timestamp();            

        averageZiggurat += (ts_end - ts_start);

//        cout << ZigguratPoly << endl;

        ts_start = get_timestamp();                        
        KnuthPoly = samplers->PolyGeneratorKnuthYao(nSamples, precision, tailcut, sigmaRR, c); // Coefficients, precision, tailcut, and sigma
        ts_end = get_timestamp();            

        averageKnuthYao += (ts_end - ts_start);

//        cout << KnuthPoly << endl;

    }//end-for

    cout << "\nZiggurat average running time for " << nIterations << " iterations: " << (float)(averageZiggurat/((float)(nIterations)*1000000000.0)) << endl;
    cout << "Knuth-Yao average running time for " << nIterations << " iterations: " << (float)(averageKnuthYao/((float)(nIterations)*1000000000.0)) << endl;

    delete(samplers);
    
    return 0;
    
}