#include <iostream>
#include <array>
#include "LTI2/matrix/math.hpp"
using namespace std;
using namespace matrix;


int main()
{
    Matrix<double, 3, 4> A;
    Matrix<double, 4, 5> B;

    A(0,0) = 1;
    A(1,1) = 1;
    A(2,2) = 1;

    B(0,0) = 1;
    B(1,1) = 2;
    B(2,2) = 3;
    B(3,3) = 4;


    auto C = A*B;

    std::cout << C << std::endl;

    return 0;

}