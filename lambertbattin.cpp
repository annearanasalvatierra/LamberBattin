/*
                          LAMBERBATTIN

   this subroutine solves Lambert's problem using Battins method. The method
   is developed in Battin (1987). It uses contiNued fractions to speed the
   solution and has several parameters that are defined differently than
   the traditional Gaussian technique.

 Inputs:         Description                    Range/Units
   ro          - IJK Position vector 1          m
   r           - IJK Position vector 2          m
   dm          - direction of motion            'pro','retro'
   Dtsec       - Time between ro and r          s

 OutPuts:
   vo          - IJK Velocity vector            m/s
   v           - IJK Velocity vector            m/s

 Reference:
 Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill, New York; 3rd edition(2007).
 
 Last modified:   2015/08/12   M. Mahooti  
*/

#include <iostream>
#include <iomanip>
#include <cstring>
#include <math.h>
#include "vector.h"
#include "seebatt.h"
#include "seebattk.h"
#include "Lambertbattin.h"

using namespace std;

void LAMBERTBATTIN(double ro[], double r[], char *dm, double Dtsec, double vo[], double v[])
{
	double small, mu, y1, magr, magro, CosDeltaNu, rcrossr[3], magrcrossr, SinDeltaNu,
		DNu, RoR, eps, tan2w, rp, L, m, x, xn, chord, s, lim1, tempx, Denom, h1, h2, b, u, k2, y,
		a, arg1, arg2, AlpH, BetH, DH, F, GDot, G, Sinv, Cosv, AlpE, BetE, am, ae, be, tm, DE;
	int n, Loops;
	small = 0.000001;
	mu=3.986004418e14; 
	y1=0;
	magr=norm(r);
	magro=norm(ro);
	CosDeltaNu= dot(ro,r)/(magro*magr);
	cross(rcrossr, n, ro, r, 3, 3); 
	magrcrossr=norm(rcrossr);
	if(strcmp(dm,"pro")==0)
	{
		SinDeltaNu=magrcrossr/(magro*magr);
	}else{
		SinDeltaNu=-magrcrossr/(magro*magr);
	}
	DNu = atan2(SinDeltaNu,CosDeltaNu);
	
	if(DNu < 0.0)
		DNu = 2.0*PI_+DNu;
	
	RoR = magr/magro;
	eps = RoR-1.0;
	tan2w=0.25*eps*eps/(sqrt(RoR)+RoR*(2.0 +sqrt(RoR)));
	rp=sqrt(magro*magr)*(pow(cos(DNu*0.25),2)+tan2w);
	
	if(DNu < PI_)
	{
		L=(pow(sin(DNu*0.25),2) +tan2w)/(pow(sin(DNu*0.25),2) +tan2w+cos(DNu *0.5));
	}else
	{
		L=(pow(cos(DNu*0.25),2) +tan2w - cos(DNu*0.5))/(pow(cos(DNu*0.25),2)+tan2w);
	}
	
	m=mu*Dtsec*Dtsec/(8.0*rp*rp*rp);
	x=10.0;
	xn=L;
	chord=sqrt(magro*magro + magr*magr -2.0*magro*magr*cos(DNu));
	s= (magro+magr+chord)*0.5;
	lim1=sqrt(m/L);
	Loops=1;
	while(1)
	{
		x=xn;
		tempx=seebatt(x);
		Denom=1.0/((1.0+2.0*x+L)*(4.0*x+tempx*(3.0+x)));
		h1 = pow(L+x,2)*(1.0+3.0*x+tempx)*Denom;
		h2 = m*(x-L+tempx)*Denom;
		 // ----------------------- EvalUAR: CÚBICA ------------------
		b = 0.25*27.0*h2 / pow((1.0+h1),3 );
		if (b < -1.0) 
		{
			xn = 1.0 - 2.0*L;
		}else
		{
			if (y1 > lim1)
			{
				xn = xn * (lim1/y1);
			}else
			{
				u = 0.5*b / ( 1.0 + sqrt( 1.0 + b ) );            
				k2 = seebattk(u);
				y = ( ( 1.0+h1 ) / 3.0 )*( 2.0 + sqrt( 1.0+b )/( 1.0+2.0*u*k2*k2 ) );
				xn= sqrt( pow((1.0-L)*0.5,2 ) + m/(y*y) ) - ( 1.0+L )*0.5;
			}
		}
		Loops = Loops + 1;
		y1=sqrt(m/((L+x)*(1.0+x)) );
		if ((fabs(xn-x) < small) && (Loops > 30))
		{
			break;
		}
	}

	a=  mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );
	
	// ------------------ Find Eccentric anomalies -----------------
	// ------------------------ HIPERBÓLICA -------------------------
	if ( a < -small )
	{
		arg1 = sqrt(s/(-2.0*a));
		arg2 = sqrt((s-chord)/(-2.0*a));
		// ------- Evaluate f and g functions --------
		AlpH = 2.0*asinh(arg1);
		BetH = 2.0*asinh(arg2);
		DH   = AlpH-BetH;
		F    = 1.0-(a/magro)*(1.0-cosh(DH));
		GDot = 1.0-(a/magr)*(1.0-cosh(DH));
		G    = Dtsec-sqrt(-a*a*a/mu)*(sinh(DH)-DH);
	}else{
		// ------------------------ ELÍPTICA ---------------------
		if ( a > small ){
			arg1 = sqrt( s / ( 2.0*a ) );
			arg2 = sqrt( ( s-chord ) / ( 2.0*a ) );
			Sinv = arg2;
			Cosv = sqrt( 1.0 - (magro+magr-chord)/(4.0*a) );
			BetE = 2.0*asin(Sinv);
			if ( DNu > PI_ ){
				BetE= -BetE;
			}
			Cosv= sqrt( 1.0 - s/(2.0*a) );
			Sinv= arg1;
			am  = s*0.5;
			ae  = PI_;
			be  = 2.0*asin( sqrt( (s-chord)/s ) );
			tm  = sqrt(am*am*am/mu)*(ae - (be-sin(be)));
			if ( Dtsec > tm )
			{
				AlpE= 2.0*PI_-2.0*asin( Sinv );
			}else
			{
				AlpE= 2.0*asin( Sinv );
			}
			DE  = AlpE - BetE;
			F   = 1.0 - (a/magro)*(1.0 - cos(DE) );
			GDot= 1.0 - (a/magr)* (1.0 - cos(DE) );
			G   = Dtsec - sqrt(a*a*a/mu)*(DE - sin(DE));
		}else{
			// --------------------- PARABÓLICA ---------------------
			arg1 = 0.0;
			arg2 = 0.0;
			cout<<" a parabolic orbit "<<endl;
		}
	}

	for (int i=0; i<3;i++)
	{
		vo[i]= ( r[i] - F*ro[i] )/G;
		v[i] = ( GDot*r[i] - ro[i] )/G;
	}
}
