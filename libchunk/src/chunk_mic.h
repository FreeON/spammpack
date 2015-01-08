#ifndef __CHUNK_MIC_H
#define __CHUNK_MIC_H

void *
malloc_huge_pages (size_t size);

void
free_huge_pages (void *ptr);

#endif
