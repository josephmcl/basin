/* profile.h */
#pragma once

#define DECLARE_PROFILE(x) void x(void);
#define DEFINE_PROFILE(x) void x(void) {
#define ENDDEF_PROFILE }

namespace test {

    template<typename TestName>
    void test(TestName f) {
        // TODO: Add TSC timer. 
        f(); }
}