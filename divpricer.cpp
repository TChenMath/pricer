#include "divpricer.hpp"
#include <cmath>	// std::abs
#include <random>	// random
#include <iostream>

void test(const std::vector<double>& ts, const std::vector<double> &ds);

int  main() 
{

	std::cout << "\nsingle div:=====\n";
	{
		std::vector<double> ts{0.1};
		std::vector<double> ds{0.5};
		test(ts, ds);		
	}

	std::cout << "\nsmall divs:=====\n";
	{
		std::vector<double> ts{0.1, 0.2};
		std::vector<double> ds{0.0001, 0.0001};
		test(ts, ds);		
	}

	std::cout << "\nlarge divs:=====\n";
	{
		std::vector<double> ts{0.1, 0.2};
		std::vector<double> ds{10, 10};
		test(ts, ds);		
	}

	return 0;
}


void test(const std::vector<double>& ts, const std::vector<double> &ds)
{
	double S0 = 100;
	double K = 100;
	double vol = 0.2;
	double T  = 0.25;
	double r = 0.02;
	double netrate = 0.02;	// r - q
	size_t npath = 200000;

	EuropeanDivPricer1 pricer1(K, vol, S0, T, r, netrate, ts, ds);

	BlackPricer black(K, vol * std::sqrt(T), S0 * std::exp(netrate * T), std::exp(-r * T));

	MonteCarloPricer mc(K, vol, S0, T, r, netrate, ts, ds, npath);

	std::cout << "black tv " << black.tv() << " mc tv " << mc.tv() << " pricer1 tv " << pricer1.tv() << "\n";
}


int Math::newtonSolve(std::function< double(double) > f,
		std::function< double(double) > fp, const int maxIter,
		const double tol, double x0, double& xstar)
{
	double x = x0;
	int i = 0;
	for (; i < maxIter; i++) {
		double dx = f(x) / fp(x);
		x -= dx;
		if (std::abs(dx) < tol) {
			xstar = x;
			return i+1;
		}
	}
	return -1;
}


double BlackPricer::tv()
{
	double d1 = std::log(fwd_ / K_) / sigma_ + 0.5 * sigma_;
	double d2 = d1 - sigma_;
	return df_ * (fwd_ * Math::normalCDF(d1) - K_ * Math::normalCDF(d2));
}


double MonteCarloPricer::tv()
{
	std::default_random_engine generator(123);
	std::normal_distribution<double> dist(0, 1);

	const size_t ndiv = ts_.size();
	// W_T - W_{t_i}
	std::vector<double> ws(ndiv);
	double payoff = 0.0;
	double z;

	for (size_t n=0; n < npath_; n++) {
		// prepare Ws
		z = dist(generator);
		ws[ndiv-1] = z * std::sqrt(T_ - ts_[ndiv-1]);
		// backwards
		for (size_t i = ndiv - 1; i-- > 0;) {
			z = dist(generator);
			ws[i] = ws[i+1] + z * std::sqrt(ts_[i+1] - ts_[i]);
		}
		// W_T
		z = dist(generator);
		const double wt = ws[0] + z * std::sqrt(ts_[0]);
		// S_T
		double res = S0_ * std::exp((netrate_ - vol_ * vol_ / 2) * T_ + vol_ * wt);
		for (size_t i = 0; i < ndiv; i++) {
			double tau = T_ - ts_[i];
			res -= ds_[i] * std::exp((netrate_ - vol_ * vol_ / 2) * tau + vol_ * ws[i]);
		}
		payoff += std::max(0.0, res - K_);
	}

	double df = std::exp(-r_ * T_);
	return df * payoff / npath_;
}


EuropeanDivPricer::EuropeanDivPricer(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	K_(K), vol_(vol), S0_(S0), T_(T), r_(r), netrate_(netrate), ts_(ts), ds_(ds)
{}

void EuropeanDivPricer::initialize()
{
	const double sqrtt = std::sqrt(T_);

	// initial guess for zstar
	double z0 = (std::log(K_ / S0_) - (netrate_ - vol_ * vol_ / 2) * T_) / (vol_ * sqrtt);

	auto f = [&](double z) -> double {
		double res = S0_ * std::exp((netrate_ - vol_ * vol_ / 2) * T_ + vol_ * sqrtt * z);
		for (size_t i = 0; i < ts_.size(); i++) {
			double t = ts_[i];
			double tau = T_ - t;
			res -= ds_[i] * std::exp((netrate_ - vol_ * vol_ / 2) * tau + 
				0.5 * vol_ * vol_ * t  * tau / T_ + vol_ * tau / sqrtt * z);
		}
		return res - K_;
	};

	auto fp = [&](double z) -> double {
		double res = S0_ * std::exp((netrate_ - vol_ * vol_ / 2) * T_ + vol_ * sqrtt * z) * vol_ * sqrtt;
		for (size_t i = 0; i < ts_.size(); i++) {
			double t = ts_[i];
			double tau = T_ - t;
			res -= ds_[i] * std::exp((netrate_ - vol_ * vol_ / 2) * tau + 
				0.5 * vol_ * vol_ * t  * tau / T_ + vol_ * tau / sqrtt * z) * vol_ * tau / sqrtt;
		}
		return res;
	};

	constexpr int maxIter = 20;
	constexpr double tol = 1e-6;
	int niter = Math::newtonSolve(f, fp, maxIter, tol, z0, zstar_);
	if (niter < 0)
		std::cout << "unable to solve for zstar\n";

	//std::cout << "z0 " << z0 << " zstar " << zstar_ << " niter " << niter << "\n";

	fwd_ = S0_ * std::exp(netrate_ * T_);
	for (size_t i = 0; i < ts_.size(); i++) {
		double t = ts_[i];
		double tau = T_ - t;
		fwd_ -= ds_[i] * std::exp(netrate_ * tau);
	}
	if (fwd_ < 1e-10)
		std::cout << "negative fwd " << fwd_ << "\n";
}

EuropeanDivPricer1::EuropeanDivPricer1(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	EuropeanDivPricer(K, vol, S0, T, r, netrate, ts, ds)
{}

void EuropeanDivPricer1::initialize()
{
	EuropeanDivPricer::initialize();
	d1_ = vol_ * std::sqrt(T_) - zstar_;
}

double EuropeanDivPricer1::tv()
{
	initialize();

	const double sqrtt = std::sqrt(T_);

	tv0_ = S0_ * std::exp(netrate_ * T_) * Math::normalCDF(d1_);
	for (size_t i = 0; i < ts_.size(); i++) {
		double t = ts_[i];
		double tau = T_ - t;
		tv0_ -= ds_[i] * std::exp(netrate_ * tau) * 
			Math::normalCDF(d1_ - vol_ * t / sqrtt);
	}
	tv0_ -= K_ * Math::normalCDF(d1_ - vol_ * sqrtt);
	double df = std::exp(-r_ * T_);
	tv0_ *= df;
	// cap/floor with theo bounds
	double lbd = std::max(0.0, df * (fwd_ - K_));
	tv0_ = std::min(std::max(lbd, tv0_), S0_ * std::exp(netrate_ - r_) * T_);
	return tv0_;
}


EuropeanDivPricer2::EuropeanDivPricer2(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	EuropeanDivPricer(K, vol, S0, T, r, netrate, ts, ds)
{}

void EuropeanDivPricer2::initialize() 
{
	EuropeanDivPricer::initialize();
}

double EuropeanDivPricer2::tv()
{
	initialize();

	size_t ndiv = ts_.size();
	std::vector<double> dstars(ndiv);

	for (size_t i = 0; i < ts_.size(); i++) {
		double t = ts_[i];
		double tau = T_ - t;
	}

	return tv0_;
}
