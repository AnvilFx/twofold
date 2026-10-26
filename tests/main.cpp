#include <cstring>

#include "test.h"

namespace {
bool g_failed = false;
}

std::vector<TestCase>& testRegistry() {
    static std::vector<TestCase> tests;
    return tests;
}

void testFailed(const char* file, int line, const char* expr) {
    std::printf("  %s:%d: CHECK(%s) failed\n", file, line, expr);
    g_failed = true;
}

int main(int argc, char** argv) {
    const char* filter = argc > 1 ? argv[1] : nullptr;
    int failures = 0, run = 0;
    for (const TestCase& t : testRegistry()) {
        if (filter && !std::strstr(t.name, filter)) { continue; }
        g_failed = false;
        t.fn();
        ++run;
        std::printf("%s %s\n", g_failed ? "FAIL" : "ok  ", t.name);
        if (g_failed) { ++failures; }
    }
    std::printf("%d of %d passed\n", run - failures, run);
    return failures == 0 ? 0 : 1;
}
