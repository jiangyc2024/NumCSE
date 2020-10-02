#include <iostream>
#include <time.h> //header for clock()

/*
 * simple Timer Class, using clock()-command from the time.h (should work on all platforms)
 * this class will only report the cputime (not walltime)
 */
class simpleTimer {
public:
	simpleTimer():time(0),bStarted(false)
	{}
	void start()
	{
		time=clock();
		bStarted=true;
	}
	double getTime()
	{
		assert(bStarted);
		return(clock()-time)/(double)CLOCKS_PER_SEC;;
	}
	void reset()
	{
		time=0;
		bStarted=false;
	}
private:
	double time;
	bool bStarted;
};
