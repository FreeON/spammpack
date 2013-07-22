/** @file
 *
 * The header file for the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLY_H
#define __MULTIPLY_H

#include "multiply.decl.h"

class MultiplyElement : public CBase_MultiplyElement
{
  public:

    MultiplyElement ();
    MultiplyElement (CkMigrateMessage *msg);
};

class Multiply : public CBase_Multiply
{
  private:

    CProxy_MultiplyElement convolution;

  public:

    Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C);
    void multiply (CkCallback &cb);
};

#endif
