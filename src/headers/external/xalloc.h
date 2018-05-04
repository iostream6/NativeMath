/*
 * 2017.04.23  - Created initial version as part of matvec.h
 * 2018.04.15  - Moved memory allocation function to separate file
 */

/* 
 * File:   xalloc.h
 * Author: Ilamah, Osho
 *
 * Created on 15 April 2018, 17:45
 */

#ifndef XALLOC_H
#define XALLOC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

    /**
     * Allocates a bit aligned memory buffer for data elements of the same type
     * @param num the number of elements expected to be stored in the memory buffer
     * @param size the number of bytes of each type of the stored element as given by sizeof
     * @param alignment the number of bits to align to 
     * @return a pointer to the allocated bit aligned memory buffer that has capacity of (num*size) bytes
     */
    void* aligned_calloc(size_t num, size_t size, size_t alignment);


    /**
     * Frees a bit aligned memory buffer that was previously allocated by aligned_calloc
     * @param p pointer to the previously allocated bit aligned memory buffer that is to be freed
     */
    void aligned_free(void* p);

#ifdef __cplusplus
}
#endif

#endif /* XALLOC_H */

