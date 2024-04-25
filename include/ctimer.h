//--------------------------------------------------------------------
//
//  FILE:  ctimer.h
//
//  Class for time testing
//
//  created: Sat Feb  4 12:37:18 MST 2012
//  updated:
//
//  Copyrigth 2002-2016 by Edmanuel Torres
//                      eeetorres@gmail.com
//--------------------------------------------------------------------
/*
 * MIT License
 *
 * Copyright (c) 2018 Edmanuel Torres
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Insert a declaration for system_time here.

#ifndef _CTIME_H_
#define _CTIME_H_

#include <ctime>
#include <iostream>

typedef clock_t system_time;
//typedef time_t system_time;

class CTimer
{
public:
    // Start and stop the timer
    void start(void)
    {
        start_time=clock();
    };
    //void start(void){ time(&start_time);};
    void stop(void)
    {
        stop_time=clock();
    };
    //void stop(void){  time(&stop_time);};

    // Compute the elapsed time (in seconds)
    double get_elapsed_time () const
    {
        return  difftime(stop_time,start_time )/CLOCKS_PER_SEC;
        //return  (double)difftime(stop_time,start_time );
    };
    // Display the elapsed time (in seconds)
    void show() const
    {
        double t =  difftime(stop_time,start_time )/CLOCKS_PER_SEC;
        std::cout<<" elapsed time = "<<t<<" s"<<std::endl;
        //return  (double)difftime(stop_time,start_time );
    };
    std::string GetDate(void)
    {
        time_t rawtime;
        char * _st, tmp_buffer[256];
        std::string stime;
        time (&rawtime);
        sprintf(tmp_buffer,"%s",ctime(&rawtime));
        _st = strtok(tmp_buffer,"\n");
        stime = _st;
        return stime;
    };

private:
    system_time start_time,   // Time that the timer was started
                stop_time;    // Time that the timer was stopped
};

#endif

//END
