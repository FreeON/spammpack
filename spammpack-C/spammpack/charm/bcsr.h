/** @file
 *
 * The header file for the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __BCSR_H
#define __BCSR_H

class BCSR
{
  private:

    int NSMat;
    int NAtoms;
    int numberNonZero;
    int numberBlocks;

    int *rowPointer;
    int *colPointer;
    int *blockPointer;
    double *matrix;

  public:

    BCSR (char *filename);
    void toDense (int *N, double **PDense);
    void toStr (void);
};

#endif
