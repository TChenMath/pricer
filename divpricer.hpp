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
	BlackPricer(int cp, double K, double sigma, double fwd, double df) : 
		cp_(cp), K_(K), sigma_(sigma), fwd_(fwd), df_(df) {}

	double tv();

	double K_, sigma_, fwd_, df_;
	int cp_;
};

struct MonteCarloPricer {
	MonteCarloPricer(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds, size_t npath) : 
	cp_(cp), K_(K), vol_(vol), S0_(S0), T_(T), r_(r), netrate_(netrate), 
	ts_(ts), ds_(ds), npath_(npath) {}

	double tv();

	double K_, vol_, S0_, T_, r_, netrate_;
	int cp_;
	std::vector<double> ts_, ds_;	// TBD: can be const reference	
	size_t npath_;
};

// base class
struct EuropeanDivPricer {

	struct Result {
		double tv_, delta_, gamma_, vega_, theta_;
	};

	EuropeanDivPricer(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds);

	virtual double tv() = 0;

	virtual Result results() = 0;

	virtual void initialize();

	void capFloorTV();
	
	void applyParityIfNeeded();

	double K_, vol_, S0_, T_, r_, netrate_;
	std::vector<double> ts_, ds_;	// TBD: can be const reference
	// d1 = vol * sqrtt - zstar
	double zstar_, fwd_, df_;
	double tv0_;
	int cp_;	// call +1, put -1
};

struct EuropeanDivPricer1 : public EuropeanDivPricer {

	EuropeanDivPricer1(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds);

	double tv() override;

	Result results() override;

	void initialize() override;

	// d1 = vol * sqrtt - zstar
	double d1_;
};

struct EuropeanDivPricer2 : public EuropeanDivPricer {

	EuropeanDivPricer2(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, 
	const std::vector<double>& ds);

	double tv() override;

	Result results() override;

	void initialize() override;

	double ustar_;
};

#endif