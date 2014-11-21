/** @file */

#ifndef __TYPED_FUNC_H
#define __TYPED_FUNC_H

#define CONCAT_2(a, b) a ## _ ## b

#define TYPED_FUNCTION(base, type) CONCAT_2(base, type)

#endif
