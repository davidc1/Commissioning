#include <sys/time.h>
#include <time.h>
class Watch {
public:
  Watch(){}
  ~Watch(){}

  void Start() {
    // Get current wall time
    struct timeval current_time;
    gettimeofday(&current_time,NULL);
    _wall_time_start = (double)current_time.tv_sec + (double)current_time.tv_usec * 1.e-6;
    // Get current cpu time
    _cpu_time_start = (double)(clock());
  }
  double WallTime() {
    // Get current wall time
    struct timeval current_time;
    gettimeofday(&current_time,NULL);
    double now = (double)current_time.tv_sec + (double)current_time.tv_usec * 1.e-6;
    // Return diff
    return (now - _wall_time_start);
  }
  double CPUTime() {
    // Get cpu time
    double now = (double)(clock());
    // Return diff
    return (now - _cpu_time_start)/CLOCKS_PER_SEC;
  }
private:
  double _cpu_time_start;
  double _wall_time_start;
};

