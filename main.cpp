#include <iostream>
#include <array>
#include "LTI/TrackingDifferentiator.hpp"
using namespace std;
using namespace LTI;


int main()
{
	Eigen::Matrix3d C;
	cout << C << endl;
	C.setIdentity();
	cout << C << endl;


	cout << inner::nchoosek(5,3) << endl;

	return 0;

}