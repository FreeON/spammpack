  .text
  .global jump_table_asm
  .type jump_table_asm, @function

jump_table_asm:
  # Push stack pointer so we can make room for local storage.
  push %rax

  mov $0x02, %rax # Move index into rax; 2 is supposed to end up at label_02.

#define SHARED
#ifdef SHARED
  lea 0(,%rax, 4), %rdx
  lea table(%rip), %rax
  mov (%rdx, %rax), %edx
  movslq %edx, %rdx
  lea table(%rip), %rax
  lea (%rdx, %rax), %rax
  jmp *%rax

  .section .rodata
  .balign 4
table:
  .long label_00-table
  .long label_01-table
  .long label_02-table
  .long label_03-table
#else
  jmp *table(, %rax, 4)

  .section .rodata
  .balign 4
table:
  .long label_00
  .long label_01
  .long label_02
  .long label_03
#endif

  .text
label_00:
  jmp done

label_01:
  jmp done

label_02:
  jmp done

label_03:
  jmp done

done:
  pop %rax
  ret

  .size jump_table_asm, .-jump_table_asm
