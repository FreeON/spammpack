/** @file
 *
 * The header file for the Add class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __ADD_H
#define __ADD_H

#include "add.decl.h"

/** An addition. */
class Add : public CBase_Add
{
  private:

    /** Matrix A. */
    CProxy_Matrix A;

    /** Matrix B. */
    CProxy_Matrix B;

  public:

    Add (CProxy_Matrix A, CProxy_Matrix B);
    void add (double alpha);
};

#endif
