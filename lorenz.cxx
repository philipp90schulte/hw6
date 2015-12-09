/**
 * Homework 6
 * Author: Philipp Schulte
 * Date: 08.12.2015
 * 
 * Description:
 * Implement the classic fourth-order Runge-Kutta scheme to integrate the Lorenz model from t=0 to t=100 with the initial 
 * conditions x(0) = 1, y(0) = 1, z(0)=1 and a=10, b= 28, c=8/3. Try different step sizes dt and compare the results, 
 * e.g. plot x(t) for different dt.
 * 
 * You can plot the three-dimensional trajectory in gnuplot via gnuplot> splot "data" u 2:3:4
 * if you have x,y,z stored in the second, third and fourth column, respectively, of your file. This should give you 
 * nice looking graphs, however they are not really suitable to compare results.
 * 
 * When you submit your code, include some nice plots of your results.
 * 
 * 
 */

// Imports and name space description
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


// Section to define the needed functions
void calcKvalues(double* k, const int a, const int b, const double c, const double x, const double y, const double z);



// main function
int main() {
	// Define variables for calculation
	const int a = 10, b = 28;
	const double c = 8.0/3.0;
	
	// Define the time variables
	const double tend = 100.0, N = 10000.0;
	double dt = tend/N;
	
	// Define the arrays which save the K values and the last calculated x, y, and z values
	double x, y, z;
	double k1[3], k2[3], k3[3], k4[3];
	double xtemp, ytemp, ztemp;
	
	// Open the stream for saving data
	string outname = "Lorenzresults.txt";
	ofstream out(outname.c_str());
	
	//Set the initial values for t = 0
	x = 1;
	y = 1;
	z = 1;
	double t = 0;
	
	out << t << "\t" << x << "\t" << y << "\t" << z << endl;
	
	for (int i = 1; i <= N; i++) {
		
		t = i * dt;
		
		// Calculate K1
		calcKvalues(k1, a, b, c, x, y, z);
				
		
		// Calculate K2
		xtemp = x + dt * 0.5 * k1[0];
		ytemp = y + dt * 0.5 * k1[1];
		ztemp = z + dt * 0.5 * k1[2];
		calcKvalues(k2, a, b, c, xtemp, ytemp, ztemp);
		
				
		// Calculate K3
		xtemp = x + dt * 0.5 * k2[0];
		ytemp = y + dt * 0.5 * k2[1];
		ztemp = z + dt * 0.5 * k2[2];
		calcKvalues(k3, a, b, c, xtemp, ytemp, ztemp);
						
						
		// Calculate K4
		xtemp = x + dt * k3[0];
		ytemp = y + dt * k3[1];
		ztemp = z + dt * k3[2];		
		calcKvalues(k4, a, b, c, xtemp, ytemp, ztemp); 
		
		
		// Calculate the new x, y, and z values
		// Example: x_i+1 = x_i + dt/6 * (K1 + 2*K2 + 2*K3 + K4)
		x = x + (dt/6) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
		y += (dt/6) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
		z += (dt/6) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
		
		out << t << "\t" << x << "\t" << y << "\t" << z << endl;
		
	}
	// Close the stream to the file
	out.close();
	
	return 0;
}

// Section to implement the functions
void calcKvalues(double* k, const int a, const int b, const double c, const double x, const double y, const double z) {
	k[0] = a*(y-x);	
	k[1] = x*(b-z)-y;
	k[2] = x*y-c*z;
}