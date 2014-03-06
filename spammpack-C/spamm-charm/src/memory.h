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
  private:
    static void parse_proc (int *VmSize, int *VmPeak);

  public:
    static int get_virtual (void);
    static int get_peak_virtual (void);
};

#endif
