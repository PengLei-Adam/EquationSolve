#include "EquationSolve.h"
#include <iostream>
using namespace std;
float myfunc( float x ) {
	return ( x*x - 13 );
}

int main() {
	EqSolve<float> solve1(myfunc, 0.4, 70.0, 1.0e-6 );
	cout << solve1.bisect() <<'\n'
		<< solve1.regfalse() <<'\n'
		<< solve1.ridder() <<'\n'
		<< solve1.false_quad() <<endl;
	return 0;
}