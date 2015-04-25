
/************************���ַ�����ⷽ��*************************/
/********************���������ڵ�������ʽ����*********************/
/************************���ַ�: bisect()*************************/
/************************��λ����regfalse()***********************/
/************************Ridders: ridder()************************/
/*****************�����β�ֵ��λ����false_quad()******************/

#ifndef _EQUATIONSOLVE_H_
#define _EQUATIONSOLVE_H_

#include <math.h>

template <typename T>
class EqSolve {
	typedef T (*FSN ) ( T );		//����ָ������ģ��
private:
	T aa;
	T bb;	//���ַ��������[aa,bb]
	float eps;	//�������
	T _root;
	FSN func;	//�ⷨʵ��

public:
	EqSolve(): aa(0), bb(10), eps(1), _root(1), func(NULL) {};
	EqSolve( FSN ,const T & ,const T & , const float & );
	~EqSolve(){};
	T root() const;		//���ظ�

	T bisect();			//���ַ�
	T regfalse();		//��λ��
	T ridder();			//Ridders
	T false_quad();		//�����β�ֵ��λ��

};
template <typename T>
EqSolve<T>::EqSolve( FSN f,const T & lt,const T & rt, const float & e )  {
	aa = lt;
	bb = rt;
	eps = e;
	func = f;

}
template <typename T>		//���ַ�ʵ��
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

template <typename T>		//��λ��ʵ��
T EqSolve<T>::regfalse() {
	T a=aa, b=bb, c, fa, fb, fc, fmin;
	do {
		fa = (*func) (a);
		fb = (*func) (b);
		c = b - (b-a)*fb/(fb-fa);
		fc = (*func) (c);
		if( fa*fc < 0 ) b = c;
		else  a = c;
		fmin =fabs( fabs(fa) < fabs(fb) ? fa : fb );		//fminȡfa��fb��С����ֵ
	} while( fabs( fabs(fc) -fmin) > eps );					//��ֵ����ֵfc��fmin֮��С��eps����ֹͣ
	return _root = c;
}

template <typename T>		//Riddersʵ��
T EqSolve<T>::ridder() {
	T x0=aa, x2=bb, f0, f1, f2,
		x1, x3, d0,
		a, b, u, v;
	T delmin;			//�µ���3��ԭʼ�����С����
	x1 = ( x0 + x2 )/2.0;	//ȷ��һ���е�
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
		//�ж���С����
		delmin = fabs( x0-x3 );
		if( fabs( x1-x3 ) < delmin ) delmin = fabs( x1-x3 );
		if( fabs( x2-x3 ) < delmin ) delmin = fabs( x2-x3 );

		x1 = x3;
		x0 = x1 - delmin;
		x2 = x1 + delmin;
	}
	return _root = x3;
}
template <typename T>	//�����β�ֵ��λ��
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