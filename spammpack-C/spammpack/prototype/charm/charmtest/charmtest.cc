#include "testmodule.decl.h"

Main::Main (CkArgMsg *msg)
{
  CkPrintf("Hello World!\n");
  CkExit();
}

Main::Main (CkMigrateMessage *msg) {}

#include "testmodule.def.h"
