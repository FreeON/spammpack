#include <stdio.h>
#include <string.h>

/** Convert a pointer address to a string.
 *
 * @param ptr The pointer.
 * @param length The maximum size of string.
 * @param string The string (must have enough space to store the string).
 */
void pointer_to_string(void *ptr, int *length, char *string)
{
    int i;

    snprintf(string, *length, "%p", ptr);
    for(i = strlen(string); i < *length; i++) {
        string[i] = ' ';
    }
}
