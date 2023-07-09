#include <iostream>
#include <array>
using namespace std;
int main()
{

	constexpr int CatchLength = 10;
	std::array<double, CatchLength> Q;
	for (int k=0; k<CatchLength; ++k)
	{
		Q[k] = k+1;
	}


	for (int k=0; k<CatchLength-1; ++k)
	{
		Q[k] = Q[k+1];
	}

	for (auto item : Q)
	{
		cout << item << endl;
	}


	return 0;

}