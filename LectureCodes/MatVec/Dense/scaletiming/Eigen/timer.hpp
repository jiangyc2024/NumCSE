# ifndef TIMER_HPP
# define TIMER_HPP

# include <chrono>

/* USAGE: Timer t;
 *        t.start(); ... thing to measure ...; t.end();
 *        double duration = t.duration();                */

class Timer{
  typedef std::chrono::high_resolution_clock clock;
  typedef std::chrono::nanoseconds prec;
  const unsigned int divisor = 1e9;
  
  private:
    clock::time_point t_start;
    clock::time_point t_end;

  public:
   void start(){ t_start = clock::now(); }
   void stop(){ t_end = clock::now(); }

   double duration()
   {  
      auto dur = std::chrono::duration_cast<prec>(t_end - t_start);
      return double(dur.count())/divisor;
   }
};

# endif
