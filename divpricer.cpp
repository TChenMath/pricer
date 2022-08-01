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
	size_t npath = 100000;

	EuropeanDivPricer1 pricer1(K, vol, S0, T, r, netrate, ts, ds);

	BlackPricer black(K, vol * std::sqrt(T), S0 * std::exp(netrate * T), std::exp(-r * T));

	MonteCarloPricer mc(K, vol, S0, T, r, netrate, ts, ds, npath);

	std::cout << "black tv " << black.tv() << " mc tv " << mc.tv() << " div tv " << pricer1.tv() << "\n";
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
	const double sqrtt = std::sqrt(T_);

	std::default_random_engine generator(123);
	std::normal_distribution<double> dist(0, 1);
	double payoff = 0.0;
	for (size_t n=0; n < npath_; n++) {
		double z = dist(generator);
		double res = S0_ * std::exp((netrate_ - vol_ * vol_ / 2) * T_ + vol_ * sqrtt * z);
		//std::cout << z << " " << vol_ << " " << sqrtt << " " << std::exp(vol_ * sqrtt * z) << " " << res << "\n";
		for (size_t i = 0; i < ts_.size(); i++) {
			double t = ts_[i];
			double tau = T_ - t;
			res -= ds_[i] * std::exp((netrate_ - vol_ * vol_ / 2) * tau + 
				0.5 * vol_ * vol_ * t  * tau / T_ + vol_ * tau / sqrtt * z);
		}
		payoff += std::max(0.0, res - K_);
	}
	double df = std::exp(-r_ * T_);
	return df * payoff / npath_;
}


double EuropeanDivPricer1::tv()
{
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


EuropeanDivPricer1::EuropeanDivPricer1(double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	K_(K), vol_(vol), S0_(S0), T_(T), r_(r), netrate_(netrate), ts_(ts), ds_(ds)
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
	double zstar;
	int niter = Math::newtonSolve(f, fp, maxIter, tol, z0, zstar);
	if (niter < 0)
		std::cout << "unable to solve for zstar\n";

	//std::cout << "z0 " << z0 << " zstar " << zstar << " niter " << niter << "\n";
	d1_ = vol_ * sqrtt - zstar;

	fwd_ = S0_ * std::exp(netrate_ * T_);
	for (size_t i = 0; i < ts_.size(); i++) {
		double t = ts_[i];
		double tau = T_ - t;
		fwd_ -= ds_[i] * std::exp(netrate_ * tau);
	}
	if (fwd_ < 1e-10)
		std::cout << "negative fwd " << fwd_ << "\n";

}