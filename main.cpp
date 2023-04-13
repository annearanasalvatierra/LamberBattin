#include <iostream>
#include <iomanip>
#include <cstring>
#include <math.h>
#include <cmath>
#include "vector.h"
#include "seebatt.h"
#include "seebattk.h"
#include "Lambertbattin.h"

using namespace std;

int main(){

	//Vector:
	double v[] = {1, 1, 1}, v1[] = {1, 2, 3}, v2[] = {1, 2, 4}, v3[] = {2, -1, 0};

	if((norm(v, 3) - sqrt(3.0)) < pow(10,-12))
		cout << "norm(v, 3) es correcto" << endl;
	else 
		cout << "Falla: norm(v, 3)" << endl;

	if((dot(v, v1, 3, 3) - 6) < pow(10,-12))
		cout << "dot(v, v1, 3, 3) es correcto" << endl;
	else 
		cout << "Falla: dot(v, v1, 3, 3)" << endl;

	int n;
	cross(v, n, v1, v2, 3, 3);
	if(v[0] == 2 and v[1] == -1 and v[2] == 0)
		cout << "cross(v, n, v1, v2, 3, 3) es correcto" << endl;
	else 
		cout << "Falla: cross(v, n, v1, v2, 3, 3)" << endl;

	//Seebatt:
	if(seebatt(0.0)-5.0<=pow(10,-14))
		cout<<"seebatt(0.0) es correcto" <<endl;
	else
		cout<<"Falla: seebatt(0.0)" <<endl;

	if(seebatt(1.0)-6.06251330587321<=pow(10,-14))
		cout<<"seebatt(1.0) es correcto" <<endl;
	else
		cout<<"Falla: seebatt(1.0)" <<endl;

	if(seebatt(10000.0)-257.144804617101<=pow(10,-12))
		cout<<"seebatt(10000.0) es correcto" <<endl;
	else
		cout<<"Falla: seebatt(10000.0)" <<endl;

	//Seebattk:
	if( seebattk(0.0)-0.333333333333333<= pow(10,-14))
		cout<<"seebattk(0.0) es correcto"<<endl;
	else
		cout<<"Falla: seebattk(0.0)"<<endl;
	
	if( seebattk(1.0)-0.290322580645161<= pow(10,-14))
		cout<<"seebattk(1.0) es correcto"<<endl;
	else
		cout<<"Falla: seebattk(1.0)"<<endl;
	
	if( seebattk(10000.0)-0.00022484822744645<= pow(10,-12))
		cout<<"seebattk(10000.0) es correcto"<<endl;
	else
		cout<<"Falla: seebattk(10000.0)"<<endl;

	//Main:
	double r1[3] = {20.0e6,20.0e6,0}, r2[3] = {-20.0e6,10.0e6,0}, tof;
	char dm[] = "retro";
	
	tof=1.0*86400;
	
	LAMBERTBATTIN(r1, r2, dm, tof, v1, v2);
	
	cout<<"("<<v1[0]<<","<<v1[1]<<","<<v1[2]<<")"<<endl;  
	cout<<"("<<v2[0]<<","<<v2[1]<<","<<v2[2]<<")"<<endl;  
    return 0;	
}