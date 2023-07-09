#include <iostream>
#include <array>
#include "LTI/TrackingDifferentiator.hpp"
using namespace std;
using namespace LTI;


int main()
{
	LTD<5> TD(2);

	cout << inner::nchoosek(5,3) << endl;

	return 0;

}