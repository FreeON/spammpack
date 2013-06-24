#define LOG_DEBUG(format, ...) logging::log(logging::DEBUG, __FILE__, __LINE__, format, ##__VA_ARGS__)
#define LOG_INFO(format, ...)  logging::log(logging::INFO,  __FILE__, __LINE__, format, ##__VA_ARGS__)
#define LOG_ERROR(format, ...) logging::log(logging::ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__)

#include <string>

class logging
{
  private:

    int loglevel;

  public:

    static const int DEBUG = 0;
    static const int INFO  = 1;
    static const int ERROR = 2;

    logging (const int loglevel);

    static void log (const int level, const char *const filename, const int line, const char *const format, ...);
    static void printBlock (const int level, const int chunksize, const std::string name, const float *const A);
};
