#include "divpricer.hpp"
#include <cmath>	// std::abs
#include <random>	// MC random
#include <iostream>

#define DEBUG 0

#define TEST_GREEKS 0


void test(const std::vector<double>& ts, const std::vector<double> &ds);

void testOne(const std::vector<double>& ts, const std::vector<double> &ds, 
			double K, double vol, double T);


int  main() 
{

	std::cout << "\nsingle div:=====\n";
	{
		const std::vector<double> ts{0.1};
		const std::vector<double> ds{1.0};
		test(ts, ds);		
	}

	std::cout << "\nsmall divs:=====\n";
	{
		const std::vector<double> ts{0.1, 0.2};
		const std::vector<double> ds{0.0001, 0.0001};
		test(ts, ds);		
	}

	std::cout << "\nlarge divs:=====\n";
	{
		const std::vector<double> ts{0.1, 0.2};
		const std::vector<double> ds{5, 10};
		test(ts, ds);		
	}

	std::cout << "\nhuge divs:=====\n";
	{
		const std::vector<double> ts{0.1, 0.2};
		const std::vector<double> ds{25, 25};
		test(ts, ds);		
	}

	return 0;
}


void test(const std::vector<double>& ts, const std::vector<double> &ds)
{
	const std::vector<double> Ks{80, 100, 120};
	const std::vector<double> vols{0.3, 0.6};
	const std::vector<double> Ts{0.25, 1.0};
#if TEST_GREEKS
	fprintf(stdout, "%10s%10s%10s%10s%10s\n",
		"tv", "delta", "gamma", "vega", "theta");
#else
	fprintf(stdout, "%10s%10s%10s%10s%10s%10s%10s%15s\n", 
		"K", "vol", "T", "Black76", "MC", "Pricer1", "Pricer2", "(MC-Pricer2)%");
#endif

	for (const auto K : Ks)
		for (const auto vol : vols)
			for (const auto T : Ts)
				testOne(ts, ds, K, vol, T);
}


void testOne(const std::vector<double>& ts, const std::vector<double> &ds, 
			double K, double vol, double T)
{	
	constexpr int cp = 1;
	constexpr double S0 = 100;
	constexpr double r = 0.02;
	constexpr double netrate = 0.01;	// r - q
	constexpr size_t npath = 200000;

	double fwd = S0 * std::exp(netrate * T);
	for (size_t i = 0; i < ts.size(); i++) {
		double tau = T - ts[i];
		fwd -= ds[i] * std::exp(netrate * tau);
	}
	BlackPricer black(cp, K, vol * std::sqrt(T), fwd, std::exp(-r * T));

	MonteCarloPricer mc(cp, K, vol, S0, T, r, netrate, ts, ds, npath);

	EuropeanDivPricer1 pricer1(cp, K, vol, S0, T, r, netrate, ts, ds);

	EuropeanDivPricer2 pricer2(cp, K, vol, S0, T, r, netrate, ts, ds);

#if TEST_GREEKS
	EuropeanDivPricer::Result res = pricer1.results();
	fprintf(stdout, "%10f%10f%10f%10f%10f\n", 
		res.tv_, res.delta_, res.gamma_, res.vega_, res.theta_);
#else
	fprintf(stdout, "%10.2f%10.2f%10.2f%10f%10f%10f%10f%15.1f\n",
		K, vol, T, black.tv(), mc.tv(), pricer1.tv(), pricer2.tv(), 
		(pricer2.tv() - mc.tv())/std::max(1e-6, pricer2.tv())*100);
#endif
}


int Math::newtonSolve(std::function< double(double) > f,
		std::function< double(double) > fp, const int maxIter,
		const double tol, double x0, double& xstar)
{
	double x = x0;
	int i = 0;
	for (; i < maxIter; i++) {
		double dx = f(x) / fp(x);
#if DEBUG		
		std::cout << "x " << x << " f " << f(x) << " fp " << fp(x) << " dx " << dx << "\n";
#endif
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
	const double d1 = std::log(fwd_ / K_) / sigma_ + 0.5 * sigma_;
	const double d2 = d1 - sigma_;
	double tv0 = df_ * (fwd_ * Math::normalCDF(d1) - K_ * Math::normalCDF(d2));
	if (cp_ == -1)
		tv0 = tv0 - df_ * (fwd_ - K_);
	return tv0;
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
		payoff += std::max(0.0, cp_ * (res - K_));
	}

	const double df = std::exp(-r_ * T_);
	return df * payoff / npath_;
}



EuropeanDivPricer::EuropeanDivPricer(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	cp_(cp), K_(K), vol_(vol), S0_(S0), T_(T), r_(r), netrate_(netrate), ts_(ts), ds_(ds)
{}


void EuropeanDivPricer::initialize()
{
	const double sqrtt = std::sqrt(T_);

	// initial guess for zstar
	double z0 = (std::log(K_ / S0_) - (netrate_ - vol_ * vol_ / 2) * T_) / (vol_ * sqrtt);

	// prepare for Newton solver
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
		std::cerr << "unable to solve for zstar\n";
#if DEBUG
	std::cout << "z0 " << z0 << " zstar " << zstar_ << " niter " << niter << "\n";
#endif
	fwd_ = S0_ * std::exp(netrate_ * T_);
	for (size_t i = 0; i < ts_.size(); i++) {
		double t = ts_[i];
		double tau = T_ - t;
		fwd_ -= ds_[i] * std::exp(netrate_ * tau);
	}
	if (fwd_ < 1e-10)
		std::cerr << "negative fwd " << fwd_ << "\n";
#if DEBUG
	std::cout << "fwd " << fwd_ << "\n";
#endif
	df_ = std::exp(-r_ * T_);
}


void EuropeanDivPricer::capFloorTV() {
	// cap/floor with theo bounds
	const double lbd = std::max(0.0, df_ * (fwd_ - K_));
	const double ubd = S0_ * std::exp((netrate_ - r_) * T_);
#if DEBUG	
	std::cout << "lbd " << lbd << " ubd " << ubd << " tv0 " << tv0_ << "\n";
#endif
	tv0_ = std::min(std::max(lbd, tv0_), ubd);
}

void EuropeanDivPricer::applyParityIfNeeded() {
	if (cp_ == -1)
		tv0_ = tv0_ - df_ * (fwd_ - K_);
}

EuropeanDivPricer1::EuropeanDivPricer1(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	EuropeanDivPricer(cp, K, vol, S0, T, r, netrate, ts, ds)
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
	tv0_ *= df_;
#if DEBUG
	std::cout << S0_ * std::exp(netrate_ * T_) * Math::normalCDF(d1_) << 
		" - " << K_ * Math::normalCDF(d1_ - vol_ * sqrtt) << " = " << tv0_ << "\n";
#endif
	capFloorTV();
	applyParityIfNeeded();	
	return tv0_;
}



EuropeanDivPricer2::EuropeanDivPricer2(int cp, double K, double vol, double S0, double T, 
	double r, double netrate, const std::vector<double>& ts, const std::vector<double>& ds) : 
	EuropeanDivPricer(cp, K, vol, S0, T, r, netrate, ts, ds)
{}


void EuropeanDivPricer2::initialize() 
{
	EuropeanDivPricer::initialize();
	const double sqrtt = std::sqrt(T_);

	const size_t ndiv = ts_.size();
	std::vector<double> dstars(ndiv);

	double rhs = 0.0;
	double dstarSum = 0.0;
	for (size_t i = 0; i < ndiv; i++) {
		double t = ts_[i];
		double tau = T_ - t;
		dstars[i] = ds_[i] * std::exp((netrate_ - vol_ * vol_ / 2) * tau + 
			0.5 * vol_ * vol_ * t  * tau / T_ + vol_ * tau / sqrtt * zstar_);
		// diagonal terms
		rhs += dstars[i] * dstars[i] * std::exp(vol_ * vol_ * t * tau / T_);
		dstarSum += dstars[i];
	}

#if DEBUG
	std::cout << "dstars: ";
	for (size_t i = 0; i < ndiv; i++)
		std::cout << dstars[i] << " ";
	std::cout << "\n";
#endif

	// cross terms
	for (size_t i = 0; i < ndiv; i++) {
		double ti = ts_[i];
		for (size_t j = i+1; j < ndiv; j++) {
			double tauj = T_ - ts_[j];
			rhs += 2 * dstars[i] * dstars[j] * std::exp(vol_ * vol_ * ti * tauj / T_);
		}
	}

	const double R = 1 + rhs / (dstarSum * dstarSum);
	if (R < 2.0)
		std::cout << "R = " << R << " < 2.0\n";

	ustar_ = 0.5 * (R + std::sqrt(R * R - 4.0));
#if DEBUG
	std::cout << "RHS = " << rhs << " R = " << R << " ustar = " << ustar_ << "\n";
#endif
}


double EuropeanDivPricer2::tv()
{
	initialize();

	const size_t ndiv = ts_.size();
	std::vector<double> ds(ndiv);

	for (size_t i = 0; i < ndiv; i++)
		ds[i] = ds_[i] * ustar_;
	EuropeanDivPricer1 p1(1, K_, vol_, S0_, T_, r_, netrate_, ts_, ds);
	const double tv1 = p1.tv();

	for (size_t i = 0; i < ndiv; i++)
		ds[i] = ds_[i] / ustar_;
	EuropeanDivPricer1 p2(1, K_, vol_, S0_, T_, r_, netrate_, ts_, ds);
	const double tv2 = p2.tv();

	tv0_ = tv1 / (1 + ustar_) + tv2 * ustar_ / (1 + ustar_);
#if DEBUG
	std::cout << "tv1 = " << tv1 << " tv2 = " << tv2 << " tv0 = " << tv0_ << "\n";
#endif
	capFloorTV();
	applyParityIfNeeded();
	return tv0_;
}


EuropeanDivPricer::Result EuropeanDivPricer2::results()
{
	Result res;
	res.tv_ = tv();

	{
		// delta, gamma
		const double ds = S0_ * 0.001;
		EuropeanDivPricer2 p1(cp_, K_, vol_, S0_ + ds, T_, r_, netrate_, ts_, ds_);
		const double tv1 = p1.tv();
		EuropeanDivPricer2 p2(cp_, K_, vol_, S0_ - ds, T_, r_, netrate_, ts_, ds_);
		const double tv2 = p2.tv();
		res.delta_ = (tv1 - tv2) / (2 * ds);
		res.gamma_ = (tv1 + tv2 - 2 * res.tv_) / (ds * ds);		
	}

	{
		// vega
		const double dvol = std::min(0.01 * vol_, 0.0025);
		EuropeanDivPricer2 p1(cp_, K_, vol_ + dvol, S0_, T_, r_, netrate_, ts_, ds_);
		const double tv1 = p1.tv();
		EuropeanDivPricer2 p2(cp_, K_, vol_ - dvol, S0_, T_, r_, netrate_, ts_, ds_);
		const double tv2 = p2.tv();
		res.delta_ = (tv1 - tv2) / (2 * dvol);
	}

	{
		// theta
		const double dt = std::min({0.004, 0.01 * T_, 0.9 * ts_[0]});
		std::vector<double> ts(ts_);
		for (size_t i = 0; i < ts.size(); i++)
			ts[i] -= dt;
		EuropeanDivPricer2 p1(cp_, K_, vol_, S0_, T_ - dt, r_, netrate_, ts, ds_);
		const double tv1 = p1.tv();
		res.theta_ = (tv1 - res.tv_) / dt;
	}
	return res;
}


EuropeanDivPricer::Result EuropeanDivPricer1::results()
{
	Result res;
	res.tv_ = tv();

	{
		// delta, gamma
		const double ds = S0_ * 0.001;
		EuropeanDivPricer1 p1(cp_, K_, vol_, S0_ + ds, T_, r_, netrate_, ts_, ds_);
		const double tv1 = p1.tv();
		EuropeanDivPricer1 p2(cp_, K_, vol_, S0_ - ds, T_, r_, netrate_, ts_, ds_);
		const double tv2 = p2.tv();
		res.delta_ = (tv1 - tv2) / (2 * ds);
		res.gamma_ = (tv1 + tv2 - 2 * res.tv_) / (ds * ds);		
	}

	{
		// one-vol vega
		const double dvol = std::min(0.01 * vol_, 0.0025);
		EuropeanDivPricer1 p1(cp_, K_, vol_ + dvol, S0_, T_, r_, netrate_, ts_, ds_);
		const double tv1 = p1.tv();
		EuropeanDivPricer1 p2(cp_, K_, vol_ - dvol, S0_, T_, r_, netrate_, ts_, ds_);
		const double tv2 = p2.tv();
		res.vega_ = (tv1 - tv2) / (2 * dvol) * 0.01;
	}

	{
		// one-day theta
		const double dt = std::min({0.004, 0.01 * T_, 0.9 * ts_[0]});
		std::vector<double> ts(ts_);
		for (size_t i = 0; i < ts.size(); i++)
			ts[i] -= dt;
		EuropeanDivPricer1 p1(cp_, K_, vol_, S0_, T_ - dt, r_, netrate_, ts, ds_);
		const double tv1 = p1.tv();
		res.theta_ = (tv1 - res.tv_) / dt / 365.25;
	}
	return res;
}


