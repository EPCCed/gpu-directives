#ifndef TIMER_H
#define TIMER_H


#include <iostream>
#include <omp.h>

struct timer
{
    timer(std::string name) : _name(name), _time(0),_start(0),_niterations(0) {}
    void start() { _start=omp_get_wtime(); }
    void stop() {_time += omp_get_wtime() - _start;_start=0;_niterations+=1;}

    auto get_time() const {return _time;}
    auto get_name() const {return _name;}
    auto get_niterations() const {return _niterations;}


    private:
    std::string _name;
    double _start;
    double _time;
    size_t _niterations;

};

#endif

std::ostream& operator<<(std::ostream & os, const timer & t);

