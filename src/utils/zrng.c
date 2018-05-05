/*
 * 2017.06.04  - Created, based on ideas from:: http://people.sc.fsu.edu/~jburkardt/c_src/ziggurat/ziggurat.html
 * 
 *  References: George Marsaglia, Wai Wan Tsang,
 *  The Ziggurat Method for Generating Random Variables, Journal of Statistical Software, Volume 5, Number 8, October 2000, seven pages.
 */

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../headers/external/zrng.h"

/**
 * Provides a congruential pseudo random number generator. (see https://en.wikipedia.org/wiki/Linear_congruential_generator)
 * and http://www.javamex.com/tutorials/random_numbers/java_util_random_algorithm.shtml#.WTQ2fuvyu5M
 * 
 * @param jcong the seed, which is updated on each call.
 * @return  the seed new value.
 */
uint32_t CONG_PRNG(uint32_t* jcong) {
    *jcong = 69069 * (*jcong) + 1234567;
    return *jcong;
}

/**
 * Returns the current reading on the CPU clock.
 * @return  the current reading of the CPU clock, in seconds.
 */
double cpu_time() {
    return (double) clock() / (double) CLOCKS_PER_SEC;
}

/**
 * Provides a "two multiply-with-carry" pseudo random number generator. (see https://en.wikipedia.org/wiki/Multiply-with-carry)
 * @param w an input seed which is updated on each call
 * @param z an input seed which is updated on each call
 * @return  the new MWC random value
 */
uint32_t MWC_PRNG(uint32_t* w, uint32_t* z) {
    *z = 36969 * (*z & 65535) + (*z >> 16);
    *w = 18000 * (*w & 65535) + (*w >> 16);
    return (*z << 16) + *w;
}

/**
 * Provides the 3-shift register (SHR3, XORSHIFT) pseudo random generator for integers. See https://en.wikipedia.org/wiki/Xorshift and
 * http://www.javamex.com/tutorials/random_numbers/xorshift.shtml#.WTQ2YOvyu5M
 * @param jsr the seed, which is updated on each call
 * @return the new SHR3 random integer value
 */
uint32_t SHR3_PRNG(uint32_t* jsr) {
    uint32_t value = *jsr;
    *jsr ^= *jsr << 13;
    *jsr ^= *jsr >> 17;
    *jsr ^= *jsr << 5;
    return value + *jsr;
}

/**
 * Provides the KISS99 pseudo random number generator (see https://eprint.iacr.org/2011/007.pdf). 
 * It essensentially combines three other prngs and updates their states concurrently
 * 
 * @param jcong an input seed which is updated on each call
 * @param jsr an input seed which is updated on each call
 * @param w an input seed which is updated on each call
 * @param z an input seed which is updated on each call
 * @return the new KISS random value
 */
uint32_t KISS99_PRNG(uint32_t* jcong, uint32_t* jsr, uint32_t* w, uint32_t* z) {
    return (MWC_PRNG(w, z) ^ CONG_PRNG(jcong)) +SHR3_PRNG(jsr);
}

/**
 * Provides a uniformly distributed single precision floating point value in the range [0,1]. This function is used by other more complex (e.g. gaussian distribution) generators
 * 
 * @param jsr the seed
 * @return 
 */
float r4_uni(uint32_t* jsr) {
    float value;
    uint32_t jsr_input = *jsr;
    *jsr ^= *jsr << 13;
    *jsr ^= *jsr >> 17;
    *jsr ^= *jsr << 5;
    value = fmod(0.5 + (float) (jsr_input + *jsr) / 65536.0 / 65536.0, 1.0);
    return value;
}

/**
 * Sets up parameters required by rngDGaussian (and rngSGaussian not yet implemeneted)
 * @param kn buffer for output
 * @param fn buffer for output
 * @param wn buffer for output
 */
void r4_nor_setup(uint32_t kn[128], float fn[128], float wn[128]) {
    double dn = 3.442619855899;
    int i;
    const double m1 = 2147483648.0;
    double q;
    double tn = 3.442619855899;
    const double vn = 9.91256303526217E-03;
    q = vn / exp(-0.5 * dn * dn);
    kn[0] = (uint32_t) ((dn / q) * m1);
    kn[1] = 0;
    wn[0] = (float) (q / m1);
    wn[127] = (float) (dn / m1);
    fn[0] = 1.0;
    fn[127] = (float) (exp(-0.5 * dn * dn));
    for (i = 126; 1 <= i; i--) {
        dn = sqrt(-2.0 * log(vn / dn + exp(-0.5 * dn * dn)));
        kn[i + 1] = (uint32_t) ((dn / tn) * m1);
        tn = dn;
        fn[i] = (float) (exp(-0.5 * dn * dn));
        wn[i] = (float) (dn / m1);
    }
    return;
}

/**
 * Generates a set of normally distributed double precision floating point values using the Ziggurat method.
 * @param seed input seed
 * @param n number of unique elements in the stream
 * @param buffer an array of double of sufficient capacity to hold the contents of the stream
 * @param mu the desired mean
 * @param sigma the desired standard deviation
 * @return 
 */
int rngDGaussian(uint32_t seed, int n, double* buffer, const double mu, const double sigma) {
    float fn[128];
    uint32_t kn[128];
    float wn[128];
    r4_nor_setup(kn, fn, wn);
    uint32_t state = seed;
    uint32_t *jsr = &state;
    for (int element = 0; element < n; element++) {
        int hz;
        uint32_t iz;
        const float r = 3.442620;
        float value;
        float x;
        float y;
        hz = (int) SHR3_PRNG(jsr);
        iz = (hz & 127);
        if (fabs(hz) < kn[iz]) {
            value = (float) (hz) * wn[iz];
        } else {
            for (;;) {
                if (iz == 0) {
                    for (;;) {
                        x = -0.2904764 * log(r4_uni(jsr));
                        y = -log(r4_uni(jsr));
                        if (x * x <= y + y) {
                            break;
                        }
                    }
                    if (hz <= 0) {
                        value = -r - x;
                    } else {
                        value = +r + x;
                    }
                    break;
                }
                x = (float) (hz) * wn[iz];
                if (fn[iz] + r4_uni(jsr) * (fn[iz - 1] - fn[iz])
                        < exp(-0.5 * x * x)) {
                    value = x;
                    break;
                }
                hz = (int) SHR3_PRNG(jsr);
                iz = (hz & 127);
                if (fabs(hz) < kn[iz]) {
                    value = (float) (hz) * wn[iz];
                    break;
                }
            }
        }
        buffer[element] = value;
    }
    if (mu != 0 || sigma != 1.0) {
        for (int element = 0; element < n; element++) {
            buffer[element] = (buffer[element] * sigma) + mu;
        }
    }
    return 0;
}

/**
 * Generates a set of uniformly distributed double precision floating point valuesin the range [0,1].
 * @param seed
 * @param n
 * @param buffer
 * @return 
 */
int rngDUniform(uint32_t seed, int n, double* buffer) {
    uint32_t state = seed;
    for (int element = 0; element < n; element++) {
        uint32_t intial_state = state;
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        buffer[element] = fmod(0.5 + (intial_state + state) / 65536.0 / 65536.0, 1.0);
    }
    return 0;
}

//
//
//  The below are currently considered to be unneeded
//  The below are currently considered to be unneeded
//  The below are currently considered to be unneeded
//
//

///**
// * Prints the current YMDHMS date as a time stamp.
// */
//void timestamp() {
//#define TIME_SIZE 40
//    static char time_buffer[TIME_SIZE];
//    const struct tm *tm;
//    size_t len;
//    time_t now;
//    now = time(NULL);
//    tm = localtime(&now);
//    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
//    printf("%s\n", time_buffer);
//    return;
//#undef TIME_SIZE
//}

///**
// * Sets up parameters required by r4_exp
// * @param ke  buffer for output
// * @param fe  buffer for output
// * @param we  buffer for output
// */
//void r4_exp_setup(uint32_t ke[256], float fe[256], float we[256]) {
//    double de = 7.697117470131487;
//    int i;
//    const double m2 = 2147483648.0;
//    double q;
//    double te = 7.697117470131487;
//    const double ve = 3.949659822581572E-03;
//    q = ve / exp(-de);
//    ke[0] = (uint32_t) ((de / q) * m2);
//    ke[1] = 0;
//    we[0] = (float) (q / m2);
//    we[255] = (float) (de / m2);
//    fe[0] = 1.0;
//    fe[255] = (float) (exp(-de));
//    for (i = 254; 1 <= i; i--) {
//        de = -log(ve / de + exp(-de));
//        ke[i + 1] = (uint32_t) ((de / te) * m2);
//        te = de;
//        fe[i] = (float) (exp(-de));
//        we[i] = (float) (de / m2);
//    }
//    return;
//}
//
///**
// * Provides exponentially distributed single precision floating point value, using the Ziggurat method. Before the first call to this function, 
// * the user must call r4_exp_setup to determine the values of KE, FE and WE
// * @param jsr the seed
// * @param ke input data computed by r4_exp_setup
// * @param fe input data computed by r4_exp_setup
// * @param we input data computed by r4_exp_setup
// * @return an exponentially distributed random value
// */
//float r4_exp(uint32_t *jsr, uint32_t ke[256], float fe[256], float we[256]) {
//    uint32_t iz;
//    uint32_t jz;
//    float value;
//    float x;
//    jz = SHR3_PRNG(jsr);
//    iz = (jz & 255);
//    if (jz < ke[iz]) {
//        value = (float) (jz) * we[iz];
//    } else {
//        for (;;) {
//            if (iz == 0) {
//                value = 7.69711 - log(r4_uni(jsr));
//                break;
//            }
//            x = (float) (jz) * we[iz];
//            if (fe[iz] + r4_uni(jsr) * (fe[iz - 1] - fe[iz]) < exp(-x)) {
//                value = x;
//                break;
//            }
//            jz = SHR3_PRNG(jsr);
//            iz = (jz & 255);
//            if (jz < ke[iz]) {
//                value = (float) (jz) * we[iz];
//                break;
//            }
//        }
//    }
//    return value;
//}
