#define LOG_INFO(format, ...) logging::log(logging::INFO, __FILE__, __LINE__, format, ##__VA_ARGS__)

class logging
{
  public:

    static const int INFO = 0;
    static void log (const int level, const char *const filename, const int line, const char *const format, ...);
};
