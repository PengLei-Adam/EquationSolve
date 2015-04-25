
/************************多种方法求解方程*************************/
/********************仅限区间内单调多项式方程*********************/
/************************二分法: bisect()*************************/
/************************试位法：regfalse()***********************/
/************************Ridders: ridder()************************/
/*****************反二次插值试位法：false_quad()******************/

#ifndef _EQUATIONSOLVE_H_
#define _EQUATIONSOLVE_H_

#include <math.h>

template <typename T>
class EqSolve {
	typedef T (*FSN ) ( T );		//函数指针类型模板
private:
	T aa;
	T bb;	//二分法解的区间[aa,bb]
	float eps;	//允许误差
	T _root;
	FSN func;	//解法实现

public:
	EqSolve(): aa(0), bb(10), eps(1), _root(1), func(NULL) {};
	EqSolve( FSN ,const T & ,const T & , const float & );
	~EqSolve(){};
	T root() const;		//返回根

	T bisect();			//二分法
	T regfalse();		//试位法
	T ridder();			//Ridders
	T false_quad();		//反二次插值试位法

};
template <typename T>
EqSolve<T>::EqSolve( FSN f,const T & lt,const T & rt, const float & e )  {
	aa = lt;
	bb = rt;
	eps = e;
	func = f;

}
template <typename T>		//二分法实现
T EqSolve<T>::bisect () {
	T a=aa, b=bb, c, fa, fb, fc;
	while ( fabs(a-b) > eps ) {
		c = (a+b) /2.0;
		fa = (*func) (a);
		fb = (*func) (b);
		fc = (*func) (c);
		if( fa*fc < 0 ) b = c;
		else  a = c;
	}
	return _root = c;
}

template <typename T>		
T EqSolve<T>::root() const {  return _root; }

template <typename T>		//试位法实现
T EqSolve<T>::regfalse() {
	T a=aa, b=bb, c, fa, fb, fc, fmin;
	do {
		fa = (*func) (a);
		fb = (*func) (b);
		c = b - (b-a)*fb/(fb-fa);
		fc = (*func) (c);
		if( fa*fc < 0 ) b = c;
		else  a = c;
		fmin =fabs( fabs(fa) < fabs(fb) ? fa : fb );		//fmin取fa，fb较小绝对值
	} while( fabs( fabs(fc) -fmin) > eps );					//插值函数值fc与fmin之差小于eps，则停止
	return _root = c;
}

template <typename T>		//Ridders实现
T EqSolve<T>::ridder() {
	T x0=aa, x2=bb, f0, f1, f2,
		x1, x3, d0,
		a, b, u, v;
	T delmin;			//新点与3个原始点的最小距离
	x1 = ( x0 + x2 )/2.0;	//确定一个中点
	while( fabs( x0 - x2 ) > eps ) {
		f0 = (*func)(x0);
		f1 = (*func)(x1);
		f2 = (*func)(x2);
		d0 = x1 - x0;
		a = (f0-f1)/(f1-f2);
		b = (f0-f1)/(f1-a*f2);
		u = (b-1)/(b+1);
		v = (a-1)/(a+1);
		x3 = x1 + d0*u*(3+u*u)/(v*(3+v*v));
		//判断最小距离
		delmin = fabs( x0-x3 );
		if( fabs( x1-x3 ) < delmin ) delmin = fabs( x1-x3 );
		if( fabs( x2-x3 ) < delmin ) delmin = fabs( x2-x3 );

		x1 = x3;
		x0 = x1 - delmin;
		x2 = x1 + delmin;
	}
	return _root = x3;
}
template <typename T>	//反二次插值试位法
T EqSolve<T>::false_quad()  {
	T a=aa, b=bb,c,d,
		fa, fb, fc, fd;
	while( fabs(a-b) > eps ) {
		c = (a+b)/2;
		fa = (*func)(a);
		fb = (*func)(b);
		fc = (*func)(c);
		if( fa != fc && fb != fc ){
			d = fa*fb*c/((fa-fc)*(fb-fc)) + fa*fc*b/((fa-fb)*(fc-fb))
				+ fb*fc*a/((fb-fa)*(fc-fa));
		}
		else
			d = fa*b/(fa-fb) + fb*a/(fb-fa);
		fd = (*func)(d);
		b = d;
		if( fd*fc < 0 ) a = c;
		else if( fd*fb < 0 ) a = b;
	}
	return _root = d;
}

#endif