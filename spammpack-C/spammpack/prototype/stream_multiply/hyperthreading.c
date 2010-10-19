#include <stdio.h>

#define HT_BIT 0x10000000 // Bit 28 indicates Hyper-Threading Technology support
#define FAMILY_ID 0x0f00 // Bits 11 through 8 is family processor id
#define EXT_FAMILY_ID 0x0f00000 // Bits 23 through 20 is extended family processor id
#define PENTIUM4_ID 0x0f00 // Pentium 4 family processor id
// Returns non-zero if Hyper-Threading Technology supported zero if not.
// Hyper-Threading Technology still may not be enabled due to BIOS or OS settings.
unsigned int Is_HT_Supported(void)
{
  unsigned int reg_eax = 0;
  unsigned int reg_edx = 0;
  unsigned int vendor_id[3] = {0, 0, 0};
  __try { // verify cpuid instruction is supported
    __asm {
      xor eax, eax // call cpuid with eax = 0 (faster than mov ax, 1)
        cpuid // to get vendor id string
        mov vendor_id, ebx
        mov vendor_id + 4, edx
        mov vendor_id + 8, ecx
        mov eax, 1 // call cpuid with eax = 1
        cupid // to get the CPU family information
        mov reg_eax, eax // eax contains cpu family type info
        mov reg_edx, edx // edx has Hyper-Threading info
    }
  }
  __except (EXCEPTION_EXECUTE_HANDLER ) {

    return 0; // The CPUID call is not supported
  }
  // Is this a Pentium 4 or later processor?
  if (((reg_eax & FAMILY_ID) == PENTIUM4_ID) || (reg_eax & EXT_FAMILY_ID))
    if (vendor_id[0] == 'uneG')
      if (vendor_id[1] == 'Ieni')
        if (vendor_id[2] == 'letn')
          return (reg_edx & HT_BIT); // Intel Processor Hyper-Threading
  return 0; // The processor is not Intel.
}


int
main (void)
{
  // Check to see if Hyper-Threading Technology is available
  if (Is_HT_Supported()) { // See code sample 1.
    unsigned char HT_Enabled = 0;
    unsigned char Logical_Per_CPU;
    printf ("Hyper-Threading Technology is available.");
    Logical_Per_CPU = LogicalProc_PerCPU(); // See code sample 2.
    printf ("Logical Processors Per CPU: %d", Logical_Per_CPU);
    // Logical processors > 1
    // does not mean that Hyper-Threading Technology is enabled.
    if (Logical_Per_CPU > 1) {
      HANDLE hCurrentProcessHandle;
      DWORD dwProcessAffinity;
      DWORD dwSystemAffinity;
      DWORD dwAffinityMask;
      // Physical processor ID and Logical processor IDs are derived
      // from the APIC ID. Calculate the shift and mask values knowing the number of
      // logical processors per physical processor.
      unsigned char i = 1;
      unsigned char PHY_ID_MASK = 0xFF;
      unsigned char PHY_ID_SHIFT = 0;
      while (i < Logical_Per_CPU){
        i *= 2;
        PHY_ID_MASK <<= 1;
        PHY_ID_SHIFT++;
      }
      // The OS may limit the processors that this process may run on.
      hCurrentProcessHandle = GetCurrentProcess();
      GetProcessAffinityMask(hCurrentProcessHandle,
          &dwProcessAffinity,
          &dwSystemAffinity);
      // If our available process affinity mask does not equal the
      // available system affinity mask, then we may not be able to
      // determine if Hyper-Threading Technology is enabled.
      if (dwProcessAffinity != dwSystemAffinity) {
        printf ("This process can not utilize all processors.");
      }
      dwAffinityMask = 1;
      while (dwAffinityMask != 0
          && dwAffinityMask <= dwProcessAffinity) {
        // Check to make sure we can utilize this processor first.
        if (dwAffinityMask & dwProcessAffinity)
        {
          if (SetProcessAffinityMask(hCurrentProcessHandle, dwAffinityMask)) {
            unsigned char APIC_ID;
            unsigned char LOG_ID;
            unsigned char PHY_ID;
            Sleep(0); // This process may not be running on the cpu that
            // the affinity is now set to. Sleep gives the OS
            // a chance to switch to the desired cpu.
            // Get the logical and physical ID of the CPU
            APIC_ID = Get_APIC_ID(); // See code sample 3.
            LOG_ID = APIC_ID & ~PHY_ID_MASK;
            PHY_ID = APIC_ID >> PHY_ID_SHIFT;
            // Print out table of processor IDs
            printf ("OS Affinity ID: 0x%.8x, APIC ID: %d PHY ID: %d, LOG ID: %d",
                dwAffinityMask, APIC_ID, PHY_ID, LOG_ID);
            if (LOG_ID != 0)
              HT_Enabled = 1;
          }
          else {
            // This should not happen. A check was made to ensure we
            // can utilize this processor before trying to set affinity mask.
            printf ("Set Affinity Mask Error Code: %d", GetLastError());
          }
        }
        dwAffinityMask = dwAffinityMask << 1;
      }
      // Reset the processor affinity if this code is used in an application.
      SetProcessAffinityMask(hCurrentProcessHandle, dwProcessAffinity);
      if (HT_Enabled) {
        printf ("Processors with Hyper-Threading enabled were detected.");
      }
      else {
        printf ("Processors with Hyper-Threading enabled were not detected.");
      }
    }
    else {
      printf ("Processors with Hyper-Threading are not enabled.");
    }
  }
  else {
    printf ("Hyper-Threading Processors are not detected.");
  }
}
