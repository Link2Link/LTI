#include "doctest/doctest.h"
#include "LTI/TrackingDifferentiator.hpp"
#include "iostream"


TEST_CASE("State Space Test : TD Order 5"){
	using namespace LTI;
	using namespace std;
	LTD<5> TD(1);
	double input = 0;
	for (int k = 0 ; k < 2000 ; ++k)
	{
		if (TD.ss.t >= 1)
			input = 10;
		TD.setU(input);
		TD.step(1E-3);
	}

	double matlab_result = 9.9967725340835;
	CHECK(TD.y(0) == doctest::Approx(matlab_result).epsilon(1E-3));
}


TEST_CASE("State Space Test : TD Order 3"){
	using namespace LTI;
	using namespace std;
	LTD<3> TD(1);
	double input = 0;
	for (int k = 0 ; k < 2000 ; ++k)
	{
		if (TD.ss.t >= 1)
			input = 10;
		TD.setU(input);
		TD.step(1E-3);

//		cout << "t : " << TD.ss.t << " y(0) " << TD.y(0) << endl;
	}

	double matlab_result = 1.2739121;
	CHECK(TD.y(0) == doctest::Approx(matlab_result).epsilon(1E-3));
}
