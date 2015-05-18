#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

typedef long long int64; 
typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

int64 GetTimeMs64();