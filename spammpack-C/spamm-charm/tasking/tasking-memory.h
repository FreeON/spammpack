#ifndef __TASKIN_MEMORY_H
#define __TASKIN_MEMORY_H

#include <stdarg.h>

#ifdef INTEL_TASK_API
#include <ittnotify.h>
static __itt_domain *domain;
#endif

void intel_task_domain_create (void)
{
#ifdef INTEL_TASK_API
  domain = __itt_domain_create("tasking domain");
#endif
}

void intel_task_begin (const char *format, ...)
{
#ifdef INTEL_TASK_API
  va_list ap;
  char task_name[100];
  va_start(ap, format);
  vsnprintf(task_name, 100, format, ap);
  __itt_string_handle *leaf_task_handle = __itt_string_handle_create(task_name);
  __itt_task_begin(domain, __itt_null, __itt_null, leaf_task_handle);
  va_end(ap);
#endif
}

void intel_task_end (void)
{
#ifdef INTEL_TASK_API
  __itt_task_end(domain);
#endif
}

#endif
