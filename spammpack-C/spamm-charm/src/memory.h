/** @file
 *
 * The header file for the Memory class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MEMORY_H
#define __MEMORY_H

/** A Memory object, containing information on memory useage.
 */
class Memory
{
  public:
    static int get_virtual (void);
};

#endif
