#include <vector>
#include <algorithm>	// std::max
#include <iostream>
#include "bspricer.hpp"
#include "bsbarrierpricer.hpp"

//#define DEBUG 1

const double eps = 1e-10;


int main()
{
	//test_TridiagonalMatrix();

	BSVanillaPricer p;
	p.test();
	
	//test_FDGridHelper();
	//FDGridHelper::test();
	
	BSBarrierPricer bp;
	bp.test();

	BSWindowBarrierPricer wp;
	wp.test();
	
	return 0;
}


void BSVanillaPricer::test()
{
	double S0 = 100;
	double K = 100;
	double r = 0.05;
	double q = 0.02;
	double sigmas[] = {0.18, 0.25, 0.32, 0.44};
	double T = 1.5;
	int theta = -1;
	
	// convergence test
	std::cout << "\nconvergence test for vanilla\n";
	std::cout << "============================\n";
	//std::cout << "accuracy  npv  npv0  err\n";
	fprintf(stdout, "%10s %10s %10s %10s\n", "accuracy", "npv", "npv0", "err");
	double sigma = sigmas[0];
	double npv0 = bspricer(S0, K, r, q, sigma, T, theta);
	int accuracy = 1;
	for (size_t ii=0; ii<4; ii++) {
		accuracy *= 2;
		double npv = calculate(S0, K, r, q, sigma, T, theta, accuracy);
		//std::cout << accuracy << "  " << npv << " " << npv0 << " " << npv-npv0 << "\n";
		fprintf(stdout, "%10d %10.4f %10.4f %10.2e\n", accuracy, npv, npv0, npv - npv0);
	}
	
	// reproduce
	accuracy = 1;
	std::cout << "\npaper table for vanilla with accuracy " << accuracy << "\n";
	std::cout << "=======================================\n";
	//std::cout << "sigma  npv  npv0  err\n";
	fprintf(stdout, "%10s %10s %10s %10s\n", "sigma", "npv", "npv0", "err");
	for (size_t ii=0; ii<4; ii++) {
		double npv = calculate(S0, K, r, q, sigmas[ii], T, theta, accuracy);
		double npv0 = bspricer(S0, K, r, q, sigmas[ii], T, theta);
		//std::cout << sigmas[ii] << " " << npv << " " << npv0 << " " << npv-npv0 << "\n";
		fprintf(stdout, "%10.2f %10.4f %10.4f %10.2e\n", sigmas[ii], npv, npv0, npv - npv0);
	}
	
}


double BSVanillaPricer::calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,
		int accuracy)
{
	// finite difference pricer for vanilla options
	
	// space grids
	size_t ns = 100 * accuracy, s0_idx;
	std::vector<double> s(ns);
	set_space_grids(S0, K, r, q, sigma, T, s, s0_idx);

	// time grids
	size_t nt = 50 * accuracy;
	double dt = T / (nt - 1);
	
	// terminal conditions
	std::vector<double> vnext(ns, 0), v(ns, 0), rhs(ns, 0);
	for (size_t ii=0; ii<ns; ii++)
		vnext[ii] = std::max(0.0, theta * (s[ii] - K));
	
	// move backwards from nt-1 to nt-2 with FIM
	// helps vanilla greeks, e.g., gamma
	const size_t num_fim_steps = 4;
	TridiagonalMatrix fim_lhs_mat(ns), fim_rhs_mat(ns);
	set_matrix_rhs_fim(s, r, q, sigma, dt/num_fim_steps, fim_lhs_mat, fim_rhs_mat);
	
	for (size_t ii=0; ii<num_fim_steps; ii++) {
		fim_rhs_mat.product(vnext, rhs);
		fim_lhs_mat.solve(rhs, v);
		v.swap(vnext);
	}
#ifdef DEBUG
	std::cout << "vnext after FIM:\n";
	print(vnext);
#endif

	// move backwards from k+1 to k
	// CN tridiagonal matrix on lhs n rhs
	TridiagonalMatrix lhs_mat(ns), rhs_mat(ns);
	set_matrix_rhs_cn(s, r, q, sigma, dt, lhs_mat, rhs_mat);
	
	for (int k=nt-3; k>=0; k--) {
		rhs_mat.product(vnext, rhs);
		lhs_mat.solve(rhs, v);
		v.swap(vnext);
	}
	v.swap(vnext);	
	
	return v[s0_idx];
	
}


void BSVanillaPricer::set_space_grids(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		std::vector<double>& s,
		size_t& idx)
{
	double smin, smax;
	set_space_range(S0, K, r, q, sigma, T, smin, smax);
	set_concentrated_grids(smin, smax, K, s);
	idx = adjust_grids(S0, s);
	
#ifdef DEBUG
	std::cout << "ns " << s.size() << " smin  " << smin << " smax " << smax << "\n";
	print(s);
#endif

}


void BSVanillaPricer::set_space_range(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double& smin,
		double& smax)
{
	double F = S0 * exp((r - q) * T);
	smax = F * exp(5.0 * sigma * sqrt(T));
	smax = std::max(1.2 * K, std::max(1.2 * S0, smax));
	
	smin = F * exp(-5.0 * sigma * sqrt(T));
	smin = std::min(0.8 * K, std::min(0.8 * S0, smin));
	
}


void BSVanillaPricer::set_concentrated_grids(
		double smin,
		double smax,
		double sfixed,
		std::vector<double>& s)
{
	// bsv setup, so that grids denser around sfixed, typically K
	
	size_t ns = s.size();
	
	// asinh(1)
	double xf = 0.881373587;
	
	// set xi[]
	double xmin = asinh(smin / sfixed);
	double xmax = asinh(smax / sfixed);
	double alpha = 0.2 * (xmax - xmin);
	double beta = asinh((xmax - xf) / alpha);
	double gamma = asinh((xmin - xf) / alpha);
	
	double xif = -gamma / (beta - gamma);
	double dxi = 1.0 / (ns - 1);
	
	std::vector<double> xi(ns);
	xi[0] = 0;
	for (size_t ii=1; ii<ns-1; ii++)
		xi[ii] = xi[ii-1] + dxi;
	xi[ns-1] = 1;
	
	// space grids
	s[0] = smin;
	for (size_t ii=1; ii<ns-1; ii++)
		s[ii] = sfixed * sinh(xf + alpha * sinh(beta * xi[ii] + gamma * (1 - xi[ii])));
	s[ns-1] = smax;
	
}


size_t BSVanillaPricer::adjust_grids(
		double x0, 
		std::vector<double>& s) const
{
	// shift so that x0 is on grid
	
	size_t i0 = 0;
	for (size_t ii=0; ii<s.size(); ii++)
		if (s[ii] > x0) {
			i0 = ii;
			break;
		}
	size_t idx = (x0 - s[i0-1] > s[i0] - x0) ? i0 : i0-1;
	s[idx] = x0;
	return idx;
	
}


void BSVanillaPricer::set_aux(
		const std::vector<double>& s,
		std::vector<std::vector<double> >& aux)
{
	// coefficients for V_ss, V_s
	
	// step sizes
	size_t n = s.size();
	std::vector<double> d(n, 0);
	for (size_t ii=1; ii<n; ii++)
		d[ii] = s[ii] - s[ii-1];
		
	// inner points
	for (size_t ii=1; ii<n-1; ii++) {
		std::vector<double>& c = aux[ii];
		c.resize(4, 0);
		double d2 = d[ii] + d[ii+1];
		
		// coeff for ii s.t. sum is zero thus not specified here
		
		// V_ss
		// coeff for ii-1
		c[0] = 2 / d[ii] / d2;
		// coeff for ii+1
		c[1] = 2 / d[ii+1] / d2;
		
		// V_s
		// coeff for ii-1
		c[2] = -d[ii+1] / d[ii] / d2;
		// coeff for ii+1
		c[3] = d[ii] / d[ii+1] / d2;
	}
	
	// boundary points
	// V_ss == 0
	// V_s approximated with one-sided difference
	
	// ii = 0
	std::vector<double>& c0 = aux[0];
	c0.resize(4, 0);
	// coeff for ii+1
	c0[3] = 1 / d[1];
	
	// ii = n-1
	std::vector<double>& cn = aux[n-1];
	cn.resize(4, 0);
	// coeff for ii-1
	cn[2] = -1 / d[n-1];
	
}


void BSVanillaPricer::set_matrix_rhs_fim(
		const std::vector<double>& s,
		double r,
		double q,
		double sigma,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat)
{
	// set up linear system for FIM scheme
	
	size_t n = s.size();
	std::vector<std::vector<double> > aux(n);
	set_aux(s, aux);
	
	std::vector<double>& ll = lhs_mat.get_a();	ll.resize(n, 0);
	std::vector<double>& ld = lhs_mat.get_b();	ld.resize(n, 0);
	std::vector<double>& lu = lhs_mat.get_c();	lu.resize(n, 0);
	std::vector<double>& rl = rhs_mat.get_a();	rl.resize(n, 0);
	std::vector<double>& rd = rhs_mat.get_b();	rd.resize(n, 0);
	std::vector<double>& ru = rhs_mat.get_c();	ru.resize(n, 0);

	// inner + boundary points
	// boundary points share the same L thus formulas apply as well
	for (size_t ii=0; ii<n; ii++) {
		const std::vector<double>& c = aux[ii];
		
		// BS equation: V_t + L * V = 0
		// L := 0.5 * sigma  * sigma * S * S * V_ss + (r-q) * S * V_s - r * V
		// FIM scheme: (1/dt - L) * V^k = 1/ dt * V^{k+1}
		
		// rhs = 1/dt * V
		rd[ii] = 1 / dt;
		
		// lhs = (1/dt - L) * V
		// lhs lower diagonal
		ll[ii] = -(0.5 * sigma * sigma * s[ii] * s[ii] * c[0] + (r-q) * s[ii] * c[2]);
		// lhs upper diagonal
		lu[ii] = -(0.5 * sigma * sigma * s[ii] * s[ii] * c[1] + (r-q) * s[ii] * c[3]);
		// lhs diagonal
		ld[ii] = 1 / dt + r - ll[ii] - lu[ii];	
	}
	
}


void BSVanillaPricer::set_matrix_rhs_cn(
		const std::vector<double>& s,
		double r,
		double q,
		double sigma,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat)
{
	// set up linear system for CN scheme
	
	size_t n = s.size();
	std::vector<std::vector<double> > aux(n);
	set_aux(s, aux);
	
	std::vector<double>& ll = lhs_mat.get_a();	ll.resize(n, 0);
	std::vector<double>& ld = lhs_mat.get_b();	ld.resize(n, 0);
	std::vector<double>& lu = lhs_mat.get_c();	lu.resize(n, 0);
	std::vector<double>& rl = rhs_mat.get_a();	rl.resize(n, 0);
	std::vector<double>& rd = rhs_mat.get_b();	rd.resize(n, 0);
	std::vector<double>& ru = rhs_mat.get_c();	ru.resize(n, 0);

	// inner + boundary points
	// boundary points share the same L thus formulas apply as well
	for (size_t ii=0; ii<n; ii++) {
		const std::vector<double>& c = aux[ii];
		
		// BS equation: V_t + L * V = 0
		// L := 0.5 * sigma  * sigma * S * S * V_ss + (r-q) * S * V_s - r * V
		// CN scheme: (2/dt - L) * V^k = (2/ dt + L) * V^{k+1}
		
		// rhs = (2/dt + L) * V
		// rhs lower diagonal
		rl[ii] = (0.5 * sigma * sigma * s[ii] * s[ii] * c[0] + (r-q) * s[ii] * c[2]);
		// rhs upper diagonal
		ru[ii] = (0.5 * sigma * sigma * s[ii] * s[ii] * c[1] + (r-q) * s[ii] * c[3]);
		// rhs diagonal
		rd[ii] = 2 / dt - r - rl[ii] - ru[ii];
		
		// lhs = (2/dt - L) * V
		ll[ii] = -rl[ii];
		lu[ii] = -ru[ii];
		ld[ii] = 2 / dt + r - ll[ii] - lu[ii];	
	}
	
}


double BSVanillaPricer::bspricer(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta) const
{
	// classic bs formula
	
	double F = S0 * exp((r-q)*T);
	double stdev = sigma * sqrt(T);
	double d1 = log(F/K) / stdev + 0.5 * stdev;
	double d2 = d1 - stdev;
	double disc = exp(-r*T);
	
	double c = disc * (F * normalCDF(d1) - K * normalCDF(d2));
	if (theta == 1)
		return theta;
	double p = c - disc * (F - K);
	return p;
	
}


// =================================================================================================


void test_TridiagonalMatrix()
{
	std::cout << "testing TridiagonalMatrix ...\n";
	
	size_t n = 150;
	TridiagonalMatrix mat(n);
	
	std::vector<double>&  a = mat.get_a();	a.resize(n, 1.2);
	std::vector<double>&  b = mat.get_b();	b.resize(n, 2.5);
	std::vector<double>&  c = mat.get_c();	c.resize(n, 1.0);
	a[0] = 0;
	c[n-1] = 0;
	
	// A * x = rhs
	std::vector<double> rhs(n), rhs2(n), x(n);
	for (size_t ii=0; ii<n; ii++)
		rhs[ii] = ii;
	
	mat.solve(rhs, x);
	//std::cout << "x:\n";
	//print(x);
	mat.product(x, rhs2);

	// check residual	
	bool success = true;
	for (size_t ii=0; ii<n; ii++) {
		double res = rhs2[ii] - rhs[ii];
		if (std::abs(res) > eps) {
			success = false;
			std::cout << ii << " x " << x[ii] << " res " << res << "\n";
		}
	}
	
	if (success)
		std::cout << "lu test passed\n";
	
}


void TridiagonalMatrix::lu_decomposte() 
{	
	l.resize(n, 0);
	u.resize(n, 0);
	d.resize(n, 0);

	d[0] = b[0];
	u[0] = c[0];
	for (size_t ii=1; ii<n; ii++) {
		l[ii] = a[ii] / d[ii-1];
		u[ii] = c[ii];
		d[ii] = b[ii] - l[ii] * u[ii-1];
	}
	
}


void TridiagonalMatrix::lu_solve(
		const std::vector<double>& rhs,
		std::vector<double>& x) const
{
	// solve L * U * x = rhs for x
	
	// L * y = rhs
	std::vector<double> y(n);
	y[0] = rhs[0];
	for (size_t ii=0; ii<n; ii++)
		y[ii] = rhs[ii] - l[ii] * y[ii-1];
		
	// U * x = y
	x.resize(n, 0);
	x[n-1] = y[n-1] / d[n-1];
	for (int ii=n-2; ii>=0; ii--)
		x[ii] = (y[ii] - u[ii] * x[ii+1]) / d[ii];
}


void TridiagonalMatrix::solve(
		const std::vector<double>& rhs,
		std::vector<double>& x)
{
	// assume a, b, c already set n will not change
	if (l.size() != n)
		lu_decomposte();
	lu_solve(rhs, x);
	
}


void TridiagonalMatrix::product(
		const std::vector<double>& v,
		std::vector<double>& rhs) const
{
	rhs[0] = b[0] * v[0] + c[0] * v[1];
	for (size_t ii=1; ii<n-1; ii++)
		rhs[ii] = a[ii] * v[ii-1] + b[ii] * v[ii] + c[ii] * v[ii+1];
	rhs[n-1] = a[n-1] * v[n-2] + b[n-1] * v[n-1];
	
}


void print(const std::vector<double>& v) 
{
	for (size_t ii=0; ii<v.size(); ii++)
		std::cout << "(" << ii << ", " << v[ii] << "), ";
	std::cout << "\n";
	
}