#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>

class timer {
public:
    using clock = std::chrono::high_resolution_clock;
    using time_point = std::chrono::time_point<clock>;
    using duration = std::chrono::duration<double>;
    timer() { reset(); }
    timer & reset() {
        _time_span = duration(0.0);
        _time_start = clock::now();
        return *this;
    }
    double operator()() const { return time(); }
    double time() const { 
        return _time_span.count(); 
    }
    timer & start() {
        _time_start = clock::now();
        return *this;
    }
    double stop() {
        _time_end = clock::now();
        _time_span += _time_end - _time_start;
        return time();
    }
    friend std::ostream & operator<<(std::ostream & o, timer & t) {
        return o << "Time: " << t.stop();
    }
private:
    time_point _time_start, _time_end;
    duration _time_span;
};

#endif
