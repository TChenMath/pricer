#ifndef DIVPRICER_HPP
#define DIVPRICER_HPP

#include <cmath>
#include <vector>
#include <functional>

struct Math {

	static int newtonSolve(std::function< double(double) > f,
		std::function< double(double) > fp, const int maxIter,
		const double tol, double x0, double& xstar);

	static double normalCDF(double x) {
		return std::erfc(-x / 1.41421356237) / 2;
	}
};

struct BlackPricer {
	BlackPricer(double K, double sigma, double fwd, double df) : 
		K_(K), sigma_(sigma), fwd_(fwd), df_(df) {}

	double tv();

	double K_, sigma_, fwd_, df_;
};

struct MonteCarloPricer {
	MonteCarloPricer(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds, size_t npath) : 
	K_(K), vol_(vol), S0_(S0), T_(T), r_(r), netrate_(netrate), 
	ts_(ts), ds_(ds), npath_(npath) {}

	double tv();

	double K_, vol_, S0_, T_, r_, netrate_;
	std::vector<double> ts_, ds_;	// TBD: can be const reference	
	size_t npath_;
};

struct EuropeanDivPricer1 {

	EuropeanDivPricer1(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds);

	double tv();

	double K_, vol_, S0_, T_, r_, netrate_;
	std::vector<double> ts_, ds_;	// TBD: can be const reference
	double d1_, fwd_;
	double tv0_;
};

#endif