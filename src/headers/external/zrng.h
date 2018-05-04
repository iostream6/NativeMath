/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   zrng.h
 * Author: Ilamah, Osho
 *
 * Created on 06 June 2017, 09:53
 */

#ifndef ZRNG_H
#define ZRNG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

    double cpu_time();

    //int rngSGaussian(int n, float* r, const float a, const float sigma);

    int rngDGaussian(uint32_t seed, int n, double* buffer, const double mu, const double sigma);

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

