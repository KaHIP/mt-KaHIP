#pragma once

#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>

static std::mutex time_mutex;
#define SYNC std::lock_guard<std::mutex> guard(time_mutex)

#define TIME

#if defined(TIME)

#define GET_MACRO12(_1, _2, NAME, ...) NAME
#define CLOCK_END(...) GET_MACRO12(__VA_ARGS__, CLOCK_END2, CLOCK_END1)(__VA_ARGS__)

#define CLOCK \
std::chrono::high_resolution_clock::now()

#define CLOCK_START \
std::chrono::time_point<std::chrono::high_resolution_clock> __begin; \
do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define CLOCK_END_THREAD_SAFE(message) \
{SYNC; std::cerr << "thread " << std::this_thread::get_id() << ": " << message << '\t' << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - __begin).count() \
<< std::endl;} while (false)

#define CLOCK_END1(message) \
{std::cout << message << '\t' << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - __begin).count() \
<< std::endl;} while (false)

#define CLOCK_END2(message, ratio) \
{std::cout << message << '\t' << std::chrono::duration<double, ratio>(std::chrono::high_resolution_clock::now() - __begin).count() \
<< std::endl;} while (false)

#define CLOCK_END_TIME \
std::chrono::duration<double>(CLOCK - __begin).count()

#define CLOCK_START_N do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define PRINT_TIME(message, time) \
{SYNC std::cout << message << '\t' << std::chrono::duration<double>(time).count() \
<< std::endl;} while (false)

#else
#define CLOCK
#define CLOCK_START
#define CLOCK_END(mes)
#define CLOCK_END_TIME
#define CLOCK_START_N
#define CLOCK_END_N(mes)
#endif