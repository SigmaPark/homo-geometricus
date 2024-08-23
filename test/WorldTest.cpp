/*  SPDX-FileCopyrightText: (c) 2019 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#include "WorldTest.hpp"
#include <iostream>


static void test() noexcept(false);


template<class A, class...ARGS>
static void log_message(A&& a, ARGS&&...args) noexcept
{
    std::wcout << static_cast<A&&>(a);
    
    if constexpr( sizeof...(ARGS) != 0 )
        log_message( static_cast<ARGS&&>(args)... );
}


static void system_pause() noexcept
{
    std::wcout << L"Press Enter key to continue. . .";
    std::cin.get();   
}


bool wt::Tests(wchar_t const* const module_title) noexcept
{
    ::log_message
    (   L"//========//========//========//========//=======#\n"
    ,   module_title, L" test Start\n"
    );

    try
    {
        ::test();

        ::log_message
        (   module_title, L" test Complete\n"
        ,   L"//========//========//========//========//=======#\n"
        );

        ::system_pause();
        
        return true;
    }
    catch(...)
    {
        std::wcout << L"Unexpected Error!\n";
    }

    ::log_message
    (   module_title, L" test Failed"
    ,   L"//========//========//========//========//=======#\n"
    );

    ::system_pause();

    return false;
}
//========//========//========//========//=======#//========//========//========//========//=======#


#include "Hamilton/Test_Hamilton.hpp"
#include "Decomposition/Test_Decomposition.hpp"
#include "Euclid/Test_Euclid.hpp"
#include "Quaternion/Test_Quaternion.hpp"
#include "Affine/Test_Affine.hpp"


void test() noexcept(false)
{
    s3d::spec::Test_Hamilton::test();
    s3d::spec::Test_Decomposition::test();
    s3d::spec::Test_Euclid::test();
    s3d::spec::Test_Quaternion::test();
    s3d::spec::Test_Affine::test();
}


int main(int const, char const* const[])
{
    wchar_t const* const os   
#if defined(_WINDOWS_SOLUTION_)
    =   L"Windows";
#elif defined(_POSIX_SOLUTION_)
    =   L"POSIX";
#else
    =   L"Others";
#endif

    return wt::Tests(os),  (int)0;
}