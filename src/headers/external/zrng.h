/*
 * 2017.06.06   - Created
 */

/* 
 * File:   zrng.h
 * Author: Ilamah, Osho
 * 
 * Header file for Ziggurat pseudo random number stream generator. The current implementation supports double precision floating point numbers.
 * 
 *  References: George Marsaglia, Wai Wan Tsang,
 *  The Ziggurat Method for Generating Random Variables, Journal of Statistical Software, Volume 5, Number 8, October 2000, seven pages.
 *
 * Created on 06 June 2017, 09:53
 */

#ifndef ZRNG_H
#define ZRNG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
    /**
     * Returns the current reading on the CPU clock.
     * @return  the current reading of the CPU clock, in seconds.
     */
    double cpu_time();

    //int rngSGaussian(int n, float* r, const float a, const float sigma);

    /**
     * Generates a set of normally distributed double precision floating point values using the Ziggurat method.
     * @param seed input seed
     * @param n number of unique elements in the stream
     * @param buffer an array of double of sufficient capacity to hold the contents of the stream
     * @param mu the desired mean
     * @param sigma the desired standard deviation
     * @return 
     */
    int rngDGaussian(uint32_t seed, int n, double* buffer, const double mu, const double sigma);

    /**
     * Generates a set of uniformly distributed double precision floating point valuesin the range [0,1].
     * @param seed
     * @param n
     * @param buffer
     * @return 
     */
    int rngDUniform(uint32_t seed, int n, double* buffer);

    //
    //    float r4_exp(uint32_t *jsr, uint32_t ke[256], float fe[256], float we[256]);
    //
    //    void r4_exp_setup(uint32_t ke[256], float fe[256], float we[256]);
    //
    //    float r4_uni(uint32_t *jsr);



#ifdef __cplusplus
}
#endif

#endif /* ZRNG_H */

