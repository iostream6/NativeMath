/*
 * 2017.04.23  - Created initial version as part of matvec.h
 * 2018.04.15  - Moved memory allocation function to separate file
 */

#include <stdlib.h>

#include "../headers/external/xalloc.h"

//https://sites.google.com/site/ruslancray/lab/bookshelf/interview/ci/low-level/write-an-aligned-malloc-free-function
//https://flyeater.wordpress.com/2010/11/29/memory-allocation-and-data-alignment-custom-mallocfree/
//http://cottonvibes.blogspot.co.uk/2011/01/dynamically-allocate-aligned-memory.html

/**
 * Allocates a bit aligned memory buffer for data elements of the same type
 * @param num the number of elements expected to be stored in the memory buffer
 * @param size the number of bytes of each type of the stored element as given by sizeof
 * @param alignment the number of bits to align to 
 * @return a pointer to the allocated bit aligned memory buffer that has capacity of (num*size) bytes
 */
void* aligned_calloc(size_t num, size_t size, size_t alignment) {
    if (alignment <= 0 || alignment % 2 != 0) {
        //force to be non zero positive multiple of 2; we just use 64 as is done in MKL
        alignment = 64;
    }
    void* block_address; // original oversized block
    void** aligned_address; // aligned data block
    int offset = alignment - 1 + sizeof (void*);
    //rather than malloc and memset (which allegedly is slow see https://developer.apple.com/library/content/documentation/Performance/Conceptual/ManagingMemory/Articles/MemoryAlloc.html),
    //we want to use calloc - which ensures initialization and is fast(er), but calloc expects num,size so we calculate the extra nums in offset and add that on - may be slightly oversized, but wouldnot matter much
    int offset_num = (offset / size) + 1;
    if ((block_address = (void*) calloc((num + offset_num), size)) == NULL) {//use calloc to ensure it is all clear!! And Also calloc is lazy initialized 
        return NULL;
    }
    aligned_address = (void**) (((size_t) (block_address) + offset) & ~(alignment - 1)); //rounding down
    aligned_address[-1] = block_address;
    return aligned_address;
}

/**
 * Free a bit aligned memory buffer that was previously allocated by aligned_calloc
 * @param p pointer to the previously allocated bit aligned memory buffer that is to be freed
 */
void aligned_free(void* p) {
    free(((void**) p)[-1]);
}