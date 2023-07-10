#ifndef LTI_PROJECT_INCLUDE_LTI_TRACKINGDIFFERENTIATOR_HPP_
#define LTI_PROJECT_INCLUDE_LTI_TRACKINGDIFFERENTIATOR_HPP_

#include "StateSpace.hpp"
namespace LTI
{

	namespace inner
	{
		static inline int fact(int N)
		{
			int result = 1;
			for (int k = 1; k <= N; ++k)
			{
				result *= k;
			}
			return result;
		}

		static inline int nchoosek(int n, int k)
		{
			return fact(n) / (fact(k) * fact(n-k));
		}

	}


	//LinearTrackingdifferentiator
	template <size_t Order=3>
	class LTD
	{
	 public:
		StateSpace<Order, 1, Order> ss;
		double freqHz;
		State<Order> y;

		void setFreqHz(double freqHz)
		{
			this->freqHz = freqHz;
			updateABC();
		}

		LTD()
		{
			setFreqHz(10);
		}
		LTD(double freqHz)
		{
			setFreqHz(freqHz);
		}

		void updateABC()
		{
			int N = Order;
			double wc = 2*EIGEN_PI*freqHz;
			double w = wc / sqrt(pow(2.0, 1.0/Order) - 1);
			Eigen::Matrix<double, Order, Order> A = Eigen::Matrix<double, Order, Order>::Zero();
			for (int k = 0; k < Order-1; ++k)
			{
				A(k, k+1) = 1;
			}

			for (int k = 0; k < Order; ++k)
			{

				A(Order-1, k) = -inner::nchoosek(N, k) * pow(w, N-k);
			}

			Eigen::Matrix<double, Order, 1> B = Eigen::Matrix<double, Order, 1>::Zero();
			B(Order-1) = pow(w, Order);


			Eigen::Matrix<double, Order, Order> C = Eigen::Matrix<double, Order, Order>::Identity();

			ss.A = A;
			ss.B = B;
			ss.C = C;
		}

		void setx0(double x0)
		{
			ss.x.Zero();
			ss.x(0) = x0;
		}

		void setU(double u)
		{
			ss.u(0) = u;
		}

		void step(double dt)
		{
			ss.step(dt);
			y = ss.y;
		}

		void reset(double x0 = 0)
		{
			ss.cleanCatch();
			setx0(x0);
		}


	};
}


#endif //LTI_PROJECT_INCLUDE_LTI_TRACKINGDIFFERENTIATOR_HPP_
