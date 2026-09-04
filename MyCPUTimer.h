#pragma once

#include <windows.h>
#include <iostream>
#include <chrono>

class MyCpuTimer {
public:
    enum Mode { Process, Thread };

    MyCpuTimer(Mode mode = Process) : mode_(mode) {
        start_ = current();

        start__ =
            std::chrono::high_resolution_clock::now();
    }

    void reset() {
        start_ = current();
    }

    double elapsed() const {
        //return current() - start_;

        auto end =  std::chrono::high_resolution_clock::now();

          double milliseconds =
            std::chrono::duration<double, std::milli>(end - start__).count();

          return milliseconds;
    }

private:
    Mode mode_;
    double start_;

    std::chrono::steady_clock::time_point start__;

    double current() const {
        FILETIME createTime, exitTime, kernelTime, userTime;
        BOOL ok = FALSE;

        if (mode_ == Process)
            ok = GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime);
        else
            ok = GetThreadTimes(GetCurrentThread(), &createTime, &exitTime, &kernelTime, &userTime);

        if (ok) {
            ULARGE_INTEGER k, u;
            k.LowPart = kernelTime.dwLowDateTime;
            k.HighPart = kernelTime.dwHighDateTime;

            u.LowPart = userTime.dwLowDateTime;
            u.HighPart = userTime.dwHighDateTime;

            return (k.QuadPart + u.QuadPart) * 1e-7; // seconds
        }

        return -1.0; // error
    }
};

 