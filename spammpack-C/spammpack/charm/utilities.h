#define LOG_INFO(format, ...)  logging::log(logging::INFO, __FILE__, __LINE__, format, ##__VA_ARGS__)
#define LOG_ERROR(format, ...) logging::log(logging::ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__)

#include <string>

class logging
{
  public:

    static const int INFO  = 0;
    static const int ERROR = 1;

    static void log (const int level, const char *const filename, const int line, const char *const format, ...);
};
