/*  SPDX-FileCopyrightText: (c) 2020 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#include <cfloat>
#include <cmath>
#include <limits>
#include "Type_Analysis/Type_Analysis.hpp"


//	C++17 or higher version of language support is required.


namespace s3d::spec
{
	
	enum class _Equiv_Number_Tag;

	template<class TAG>  
	struct _Equivalent;

}


template<>
struct s3d::spec::_Equivalent<s3d::spec::_Equiv_Number_Tag>
{
private:
	template<class T>
	static auto constexpr _epsilon() noexcept
	{
		static_assert(std::numeric_limits<T>::is_iec559, "NOT a Real Number.");

		auto constexpr res 
		=	sgm::is_Same<T, double>::value ? T(DBL_EPSILON)
		:	sgm::is_Same<T, long double>::value ? T(LDBL_EPSILON)
		:	/* otherwise */ T(FLT_EPSILON);

		return res;
	}

public:
	template<class L, class R, class...TYPES>
	static bool calc(L Lhs, R rhs, TYPES...args)
	{
		if constexpr( sizeof...(TYPES) > 0 )
			return calc(Lhs, rhs) && calc(Lhs, args...);
		else if constexpr(std::numeric_limits<L>::is_integer || sgm::is_Pointer<L>::value)
			return Lhs == rhs;
		else if constexpr (std::numeric_limits<L>::is_iec559)
			return std::abs(Lhs - rhs) < static_cast<L>(1e3) * _epsilon<L>();
		else
			return sgm::Compile_Fails(); // "no method to compare them."
	}
};
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#