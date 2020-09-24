/*
 * unixTimer Class using times()-command from the unixbased times.h
 * this class will report the user, system and real time.
 */
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
class unixtimer {
public:
	unixtimer():utime(0),stime(0),rtime(0),bStarted(false)
	{}
	
  void start() {rt0=times(&t0);	bStarted=true;	}
	
	double stop() {		
		tms t1;
		long rt1;
		assert(bStarted);
		rt1=times(&t1);
		utime=((double)(t1.tms_utime-t0.tms_utime))/ CLOCKS_PER_SEC*10000;
		stime=((double)(t1.tms_stime-t0.tms_stime))/ CLOCKS_PER_SEC*10000;
		rtime=((double)(rt1-rt0))/ CLOCKS_PER_SEC*10000;
		bStarted=false;
		return rtime;
	}
	
	double user() { assert(!bStarted); return utime;}
	double system(){assert(!bStarted); return stime;}
	double real(){assert(!bStarted); return rtime;}
	
private:
	double utime,stime,rtime;
	tms t0;
	long rt0;
	bool bStarted;
};

