#pragma once
#ifndef FEM_TIMER_H
#define FEM_TIMER_H
#include <windows.h>

class HighResolutionTimer
{
public:
    virtual void set_start() = 0;
    virtual void set_end() = 0;
    virtual float get_millisecond() = 0;
};

class HighResolutionTimerForWin : public HighResolutionTimer
{
public:

    HighResolutionTimerForWin() {
        QueryPerformanceFrequency(&freq_);
        start_.QuadPart = 0;
        end_.QuadPart = 0;
    }

    void set_start() {
        QueryPerformanceCounter(&start_);
    }

    void set_end() {
        QueryPerformanceCounter(&end_);
    }

    float get_millisecond() {
        return static_cast<float>((end_.QuadPart - start_.QuadPart) * 1000 / (float)freq_.QuadPart);
    }

private:
    LARGE_INTEGER freq_;
    LARGE_INTEGER start_, end_;
};


#endif
