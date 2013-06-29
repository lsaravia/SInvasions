#ifndef __SMATTPL_HPP
#define __SMATTPL_HPP

#include <iostream>
#include <process.h>

//#ifdef __GNUC__
//#include "fortify.h"
//#endif

#define INTEGER int


template <class Type> class simplmat
	{
	INTEGER rows,cols;
	Type *x;
	Type dummy;
	public:
	simplmat() {rows=0; cols=0; x=NULL; dummy=0; };
	~simplmat() { if( x!=NULL) delete []x; };

	simplmat(INTEGER r, INTEGER c=1);

	INTEGER getRows() {return rows;};
	INTEGER getCols() {return cols;};

	void resize(INTEGER r, INTEGER c=1);
	void resize(INTEGER r, INTEGER c, Type fval);
	Type & elem(INTEGER i, INTEGER j=0)
	{
	#ifdef NOWRAP
		if ( i < 0	||	rows <= i  ||  j < 0	||	cols <= j )
			return dummy;

//         {
//         cout << i << "\t" << j << "\t" << dummy << endl;
//			return dummy;
//         }
	#endif
	#ifdef _DEBUG
		if ( i < 0	||	rows <= i  ||  j < 0	||	cols <= j )
			error("Error, simplmat, check2, Index out of range");
	#endif

		return x[i + rows*j ];
	};

	Type & operator()(INTEGER i){ return elem(i,0); };

	Type & operator()(INTEGER i, INTEGER j){ return elem(i,j); };

	void fill(Type fval);
	void error(const char* msg);
	};

template <class Type> simplmat<Type>::simplmat(INTEGER r, INTEGER c)
	{
	if (r<=0) error("Error, simplmat, simplmat, Number of rows not positive");
	if (c<=0) error("Error, simplmat, simplmat, Number of columns not positive");
	dummy=0;
	rows=r;
	cols=c;
	x = new Type[r*c];
	if (x == 0) error("Error, simplmat, simplmat, Operator new failed");
	}

template <class Type> void simplmat<Type>::resize(INTEGER r, INTEGER c)
	{
	if( r==rows && c==cols )
		return;

	if (x != NULL)
		delete x;
	rows=r;
	cols=c;
	x = new Type[r*c];
	if (x == NULL)
		error ("Error, simplmat, resize, Operator new failed");
	}

template <class Type> void simplmat<Type>::resize(INTEGER r, INTEGER c, Type fval)
	{
	resize(r,c);
	fill(fval);
	}

template <class Type> void simplmat<Type>::error(const char* msg)
	{
	cerr << msg << "\n";
	exit(1);
	}

template <class Type> void simplmat<Type>::fill(Type fval)
	{
	INTEGER len = rows*cols;
	Type* top = &(x[len]);
	Type* t = x;
	while (t < top) *t++ = fval;
	}

#endif
