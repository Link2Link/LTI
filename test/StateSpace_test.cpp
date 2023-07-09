#include "doctest/doctest.h"
#include "LTI/StateSpace.hpp"
#include "iostream"


TEST_CASE("State Space Test : Second Order Integrator"){
	using namespace LTI;
	StateSpace<2, 1, 1> ss;
	ss.A << 0, 1,
		0, 0;
	ss.B << 0, 1;
	ss.C << 1, 0;

	ss.x(0) = 0;
	ss.u(0) = 1;

	State<1> y;

	for (; ss.t < 5; )
	{
		ss.step(1E-3);
//        std::cout << " t: " << ss.t << " y: "  << ss.y << std::endl;
	}
	y(0) = 0.5* pow(ss.t, 2);   // 解析表达式 y = 1/2*t^2
	CHECK(ss.y(0) == doctest::Approx(y(0)).epsilon(1E-3));
}

TEST_CASE("State Space Test : Linear Check"){
	using namespace LTI;
	StateSpace<2, 1, 1> ss;
	ss.A << -1, -2,
		-3, -4;
	ss.B << -1, 1;
	ss.C << 1, 0;

	ss.x(0) = 0;

	auto ss2 = ss;
	auto ss3 = ss;

	for (; ss.t < 5; )
	{
		ss.u(0) = sin(ss.t);
		ss2.u(0) = cos(ss2.t);
		ss3.u(0) = sin(ss.t) + 2 * cos(ss2.t);

		ss.step(1E-3);
		ss2.step(1E-3);
		ss3.step(1E-3);
	}
	CHECK(ss.y(0)+ss2.y(0)*2 == doctest::Approx(ss3.y(0)).epsilon(1E-3));    // Check linearity

//    std::cout << ss.y(0)+ss2.y(0)*2 << " |||| "  << ss3.y(0) << std::endl;

}

TEST_CASE("State Space Test : Catch Test"){
	using namespace LTI;
	using namespace std;
	StateSpace<2, 1, 1, 10> ss;
	ss.A << 0, 1,
		0, 0;
	ss.B << 0, 1;
	ss.C << 1, 0;

	ss.x(0) = 0;
	ss.u(0) = 1;

	CHECK((ss.x_catch[0] == Eigen::Vector2d::Zero()));
	CHECK((ss.x_catch[9] == Eigen::Vector2d::Zero()));

	CHECK((ss.u_catch[9] == Eigen::Vector<double, 1>::Zero()));
	CHECK((ss.y_catch[9] == Eigen::Vector<double, 1>::Zero()));
	CHECK(ss.t_catch[9] == 0);

	State<1> y_last;
	State<2> x_last;
	State<1> u_last;
	double dt = 1E-3;

	for (; ss.t < 5; )
	{
		y_last = ss.y;
		x_last = ss.x;
		u_last = ss.u;

		ss.step(dt);
	}

	CHECK(ss.t_catch[0] == doctest::Approx(ss.t).epsilon(1E-6));    // Check linearity
	CHECK(ss.t_catch[9] == doctest::Approx(ss.t-9*dt).epsilon(1E-6));    // Check linearity

	CHECK(ss.x_catch[1] == x_last);
	CHECK(ss.y_catch[1] == y_last);
	CHECK(ss.u_catch[1] == u_last);

}

TEST_CASE("State Space Test : Catch Default Test"){
	using namespace LTI;
	using namespace std;
	StateSpace<2, 1, 1> ss;
	ss.A << 0, 1,
		0, 0;
	ss.B << 0, 1;
	ss.C << 1, 0;

	ss.x(0) = 0;
	ss.u(0) = 1;

	CHECK((ss.x_catch[0] == Eigen::Vector2d::Zero()));
	CHECK((ss.x_catch[1] == Eigen::Vector2d::Zero()));

	CHECK((ss.u_catch[1] == Eigen::Vector<double, 1>::Zero()));
	CHECK((ss.y_catch[1] == Eigen::Vector<double, 1>::Zero()));
	CHECK(ss.t_catch[1] == 0);

	State<1> y_last;
	State<2> x_last;
	State<1> u_last;
	double dt = 1E-3;

	for (; ss.t < 5; )
	{
		y_last = ss.y;
		x_last = ss.x;
		u_last = ss.u;

		ss.step(dt);
	}

	CHECK(ss.t_catch[0] == doctest::Approx(ss.t).epsilon(1E-6));    // Check linearity
	CHECK(ss.t_catch[1] == doctest::Approx(ss.t-dt).epsilon(1E-6));    // Check linearity

	CHECK(ss.x_catch[1] == x_last);
	CHECK(ss.y_catch[1] == y_last);
	CHECK(ss.u_catch[1] == u_last);
}


TEST_CASE("State Space Test : Catch Index Safety"){
	using namespace LTI;
	using namespace std;
	StateSpace<2, 1, 1> ss;
	ss.A << 0, 1,
		0, 0;
	ss.B << 0, 1;
	ss.C << 1, 0;

	ss.x(0) = 0;
	ss.u(0) = 1;

	double dt = 1E-3;

	State<2> x_last;
	for (; ss.t < 5; )
	{
		x_last = ss.x;
		ss.step(dt);
	}

	CHECK(ss.getState(1) == ss.getState(-1));
	CHECK(ss.getState(1) == ss.getState(2));
	CHECK(ss.getState(10) == ss.getState(20));
	CHECK(x_last == ss.getState(-1));
	CHECK(ss.getState(-10) == ss.getState(10));

}