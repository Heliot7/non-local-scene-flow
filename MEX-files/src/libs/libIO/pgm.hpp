#ifndef _PGM_H_
#define _PGM_H_

#include <iostream>

#define HI(num)	(((num) & 0x0000FF00) >> 8)
#define LO(num)	((num) & 0x000000FF)

typedef struct _PGMData {
    int row;
    int col;
    int max_gray;
    int **matrix;
} PGMData;

PGMData* readPGM(const char *file_name, PGMData *data);
void writePGM(const char *filename, const PGMData *data);

/* Private functions */
void SkipComments(FILE *fp);
int **allocate_dynamic_matrix(int row, int col);
void deallocate_dynamic_matrix(int **matrix, int row);

#endif
