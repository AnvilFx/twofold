#pragma once

#include <cstdio>
#include <vector>

struct TestCase {
    const char* name;
    void (*fn)();
};

std::vector<TestCase>& testRegistry();
void testFailed(const char* file, int line, const char* expr);

#define TEST(name)                                                                     \
    static void name();                                                                \
    static const bool name##_registered = (testRegistry().push_back({ #name, name }), true); \
    static void name()

#define CHECK(cond)                                    \
    do {                                               \
        if (!(cond)) {                                 \
            testFailed(__FILE__, __LINE__, #cond);     \
            return;                                    \
        }                                              \
    } while (0)

#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))
