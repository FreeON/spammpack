/* vim: syntax=c
 */

#define COMPUTE_LOOP_BODY(c_row) temp_c##c_row += temp_a*b[comp+c_row*16]
#define PREFIX_LOOP_BODY(c_row) float temp_c##c_row = 0.0f
#define POSTFIX_LOOP_BODY(c_row) foreach (c_col = 0 ... 16) { c[c_row*16+c_col] += temp_c##c_row; }

export void sgemm_ispc_16x16_ (uniform float* uniform c, uniform float* uniform a, uniform float* uniform b)
{
  float temp_c[16];
  PREFIX_LOOP_BODY(0);
  PREFIX_LOOP_BODY(1);
  PREFIX_LOOP_BODY(2);
  PREFIX_LOOP_BODY(3);
  PREFIX_LOOP_BODY(4);
  PREFIX_LOOP_BODY(5);
  PREFIX_LOOP_BODY(6);
  PREFIX_LOOP_BODY(7);
  PREFIX_LOOP_BODY(8);
  PREFIX_LOOP_BODY(9);
  PREFIX_LOOP_BODY(10);
  PREFIX_LOOP_BODY(11);
  PREFIX_LOOP_BODY(12);
  PREFIX_LOOP_BODY(13);
  PREFIX_LOOP_BODY(14);
  PREFIX_LOOP_BODY(15);

  for (uniform int comp = 0; comp < 16; comp++) {
    foreach (c_col = 0 ... 16) {

      float temp_a = a[c_col+comp*16];

      COMPUTE_LOOP_BODY(0);
      COMPUTE_LOOP_BODY(1);
      COMPUTE_LOOP_BODY(2);
      COMPUTE_LOOP_BODY(3);
      COMPUTE_LOOP_BODY(4);
      COMPUTE_LOOP_BODY(5);
      COMPUTE_LOOP_BODY(6);
      COMPUTE_LOOP_BODY(7);
      COMPUTE_LOOP_BODY(8);
      COMPUTE_LOOP_BODY(9);
      COMPUTE_LOOP_BODY(10);
      COMPUTE_LOOP_BODY(11);
      COMPUTE_LOOP_BODY(12);
      COMPUTE_LOOP_BODY(13);
      COMPUTE_LOOP_BODY(14);
      COMPUTE_LOOP_BODY(15);
    }
  }

  POSTFIX_LOOP_BODY(0);
  POSTFIX_LOOP_BODY(1);
  POSTFIX_LOOP_BODY(2);
  POSTFIX_LOOP_BODY(3);
  POSTFIX_LOOP_BODY(4);
  POSTFIX_LOOP_BODY(5);
  POSTFIX_LOOP_BODY(6);
  POSTFIX_LOOP_BODY(7);
  POSTFIX_LOOP_BODY(8);
  POSTFIX_LOOP_BODY(9);
  POSTFIX_LOOP_BODY(10);
  POSTFIX_LOOP_BODY(11);
  POSTFIX_LOOP_BODY(12);
  POSTFIX_LOOP_BODY(13);
  POSTFIX_LOOP_BODY(14);
  POSTFIX_LOOP_BODY(15);
}
