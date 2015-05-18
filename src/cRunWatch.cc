/*

Copyright (c) 2010 James Bremner, Raven's Point Consulting, ravenspoint@yahoo.com

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "stdafx.h"
#ifndef raven_set_crunwatch
#include "ravenset/cRunWatch.h"
#endif
using namespace std;
using namespace boost::accumulators;

namespace raven {
	namespace set {

		int cRunWatch::myFlagProfiling = 0;
		map < string, accumulator_set<double, stats<tag::variance(lazy)> > > cRunWatch::myMap;

		boost::mutex cRavenProfileMutex;

		cRunWatch::cRunWatch( const char* name )
		{
			if( myFlagProfiling ) {
				myFlagTimer = 1;
				myName = std::string( name );
				QueryPerformanceCounter( (LARGE_INTEGER *)&myTimeStart );
			} else {
				myFlagTimer = 0;
			}
		}


		/**

		Destructor

		If the timer is running, stop timer and store in accumulator

		*/
		cRunWatch::~cRunWatch()
		{
			// check timer is running
			if( ! myFlagTimer )
				return;

			// get duration of this scope
			__int64 t;
			QueryPerformanceCounter( (LARGE_INTEGER *)&t );
			t -= myTimeStart;

			// ensure that no other thread is updating accumulors
			boost::mutex::scoped_lock
				lock(cRavenProfileMutex);

			// find accumulator for this scope
			map_iter_t p = myMap.find( myName );
			if( p == myMap.end() ) {
				// this is the first time this scope has run
				p = myMap.insert( pair<string,acc_t > (myName, acc_t()) ).first;
			}
			// add the time of running to the accumulator for this scope
			(p->second)( (double) t );

		}
		/**

		Extract statistics from accumulator

		@param[in] f clock frequency
		@param[in] p iterator into map of accumumators

		*/
		void cRunWatch::Tally(__int64& f, map_iter_t& p)
		{
			aver	= ( (float) mean(p->second) / f );
			//stdev	= ( (float) sqrt( ((double) variance(p->second))  ) / f );
			count	= ( boost::accumulators::count(p->second) );
			myTotalTime = count * aver;

		}

		FILE * cRunWatch::fp;

		void cRunWatch::Print() const
		{ 
			if( ! fp )
				printf("%25s %8d\t%f\t%f\n",
				myName.c_str(), count, aver, count*aver);
			else
				fprintf(fp,"%25s %8d\t%f\t%f\n",
				myName.c_str(), count, aver, count*aver);
		}

		/**

		Generate profile report

		@param[in] filename Name of result log file to create, default standard out

		This is a static function, so it can be invoked without an instance of the profiler.


		*/
		void cRunWatch::Report( const wchar_t * filename )
		{
			if( ! myFlagProfiling ) {
				return;
			}

			// number of ticks per second
			__int64 f;
			QueryPerformanceFrequency( (LARGE_INTEGER *)&f );


			// extract final statistivs for all scopes
			multiset < cRunWatch,less<cRunWatch> > sdata;
			for( map_iter_t p = myMap.begin();
				p != myMap.end(); p++ )
			{
				sdata.insert( cRunWatch( p, f ) );
			}

			// print
			fp = 0;
			if( filename )
				if( _wfopen_s(&fp,filename,L"w") )
					return;
			if( ! fp ) {
				printf("raven::set::cRunWatch code timing profile\n");
				printf("%25s   Calls\tMean (secs)\tTotal\n","Scope");
			} else {
				fprintf(fp,"raven::set::cRunWatch code timing profile\n");
				fprintf(fp,"%25s   Calls\tMean (secs)\tTotal\n","Scope");
			}
			for_each( sdata.begin(), sdata.end(),
				boost::mem_fn( &cRunWatch::Print ) );
			if( fp )
				fclose(fp);
		}

	}
}
