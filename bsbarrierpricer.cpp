#include <iostream>
#include "bsbarrierpricer.hpp"

//#define DEBUG 1

const double eps = 1e-10;


void BSWindowBarrierPricer::test()
{
	// parameters per Table 1 from "Window Double Barrier Options" by Tristan
	// Table 1 wbarrier prices incorrect, thought formula produces values consistent with FD scheme

	double S0 = 3191.23;
	double K = 3125;
	double r = 0.02;
	double q = 0.0;
	double sigmas[] = {0.2, 0.25, 0.32, 0.44};
	double T = 0.175342;
	double Ls[] = {2850, 70, 60};
	double Us[] = {9120, 130, 140};
	double T1 = 0.0;
	double T2 = 0.02;
	int theta = -1;
	
	// convergence test
	std::cout << "\nconvergence test for window barrier\n";
	std::cout << "===================================\n";
	fprintf(stdout, "%10s %10s %10s %10s\n", "accuracy", "npv", "npv0", "err");
	double sigma = sigmas[0];
	double L = Ls[0];
	double U = Us[0];
	double npv0 = calculate(S0, K, r, q, sigma, T, theta, L, U, T1, T2, 16);
	
	int accuracy = 1;
#ifdef DEBUG
	for (size_t ii=0; ii<1; ii++) {
#else
	for (size_t ii=0; ii<4; ii++) {
#endif
		accuracy *= 2;
		double npv = calculate(S0, K, r, q, sigma, T, theta, L, U, T1, T2, accuracy);
		fprintf(stdout, "%10d %10.4f %10.4f %10.2e\n", accuracy, npv, npv0, npv - npv0);
	}
	
	// reproduce
	accuracy = 2;
	std::cout << "\npaper table for window barrier with accuracy " << accuracy << "\n";
	std::cout << "==============================================\n";
	fprintf(stdout, "%7s %10s %10s %10s %10s\n", "L/U", "sigma", "npv", "npv0", "err");
	for (size_t ii=0; ii<3; ii++) {
		double L = Ls[ii];
		double U = Us[ii];
		for (size_t jj=0; jj<4; jj++) {
			double npv = calculate(S0, K, r, q, sigmas[jj], T, theta, L, U, T1, T2, accuracy);
			double npv0 = calculate(S0, K, r, q, sigmas[jj], T, theta, L, U, T1, T2, 16);
			fprintf(stdout, "%3.0f/%3.0f %10.2f %10.4f %10.4f %10.2e\n", 
					L, U, sigmas[jj], npv, npv0, npv - npv0);
		}
	}
}


double BSWindowBarrierPricer::calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,
		double L,
		double U,
		double T1,
		double T2,
		int accuracy)
{
	// finite difference pricer for double window-barrier KO
	
	// space grids, concentrated around K, L, U
	// in particular, L, U are on grids
	size_t ns = accuracy * std::min(130*1.5, 130+20*sigma*sigma*T);
	std::vector<double> s(ns);
	// s0_idx is the index for S0
	size_t s0_idx;
	set_space_grids(S0, K, r, q, sigma, T, L, U, s, s0_idx);
	
	// time grids, mostly regular dt except around T1, T2
	// in particular, T1, T2 are on grids
	size_t nt = accuracy * std::min(50*1.5, 50+2*T);
	double dt = T / (nt - 1);
	std::vector<double> ts(nt);
	size_t t1_idx, t2_idx;
	set_time_grids(T1, T2, dt, ts, t1_idx, t2_idx);
	
	// terminal conditions
	std::vector<double> vnext(ns, 0), v(ns, 0), rhs(ns, 0);
	for (size_t ii=0; ii<ns; ii++)
		vnext[ii] = std::max(0.0, theta * (s[ii] - K));
	apply_barriers(L, U, s, T1, T2, T, vnext);	

	// CN tridiagonal matrix on lhs n rhs
	// without barriers
	TridiagonalMatrix lhs_mat(ns), rhs_mat(ns), blhs_mat(ns), brhs_mat(ns);
	set_matrix_rhs_cn(s, r, q, sigma, dt, lhs_mat, rhs_mat);
	// with barriers
	set_matrix_rhs_cn(s, r, q, sigma, dt, blhs_mat, brhs_mat);
	update_matrix_rhs(L, U, s, dt, blhs_mat, brhs_mat);

	for (int k=nt-2; k>=0; k--) {
		double t = ts[k];
		double tnext = ts[k+1];
		
		// apply FIM on entering n exiting barrier window
		// critical to achieve 2nd-order convergence
		if (k == nt-2 || k == t2_idx-1 || k == t1_idx-1)
			// discontinuous v
			fim_step(r, q, sigma, L, U, tnext-t, s, vnext, v, rhs, is_monitored(T1, T2, t));
		
		else if (k == t2_idx) {
			// not monitored except at T2
			cn_step(r, q, sigma, L, U, tnext-t, s, vnext, v, rhs, false);
			apply_barriers(L, U, s, T1, T2, t, vnext);	
				
		} else if (k  == t1_idx) {
			// irregular dt
			cn_step(r, q, sigma, L, U, tnext-t, s, vnext, v, rhs, true);
		
		} else if (is_monitored(T1, T2, t)) {
			brhs_mat.product(vnext, rhs);
			blhs_mat.solve(rhs, v);
			v.swap(vnext);
		
		} else {
			rhs_mat.product(vnext, rhs);
			lhs_mat.solve(rhs, v);
			v.swap(vnext);
		}
	}
	v.swap(vnext);	
	
	return v[s0_idx];

}


void BSWindowBarrierPricer::set_time_grids(
		double T1,
		double T2,
		double dt,
		std::vector<double>& ts,
		size_t& t1_idx,
		size_t& t2_idx) const
{
	// regular grids
	for (size_t ii=0; ii<ts.size(); ii++)
		ts[ii] = ii * dt;
	t1_idx = adjust_grids(T1, ts);
	t2_idx = adjust_grids(T2, ts);
	
#ifdef DEBUG
	std::cout << "t1 idx " << t1_idx << " t2 idx " << t2_idx << "\n";
	print(ts);
#endif
	
}


bool BSWindowBarrierPricer::is_monitored(
		double T1, 
		double T2, 
		double t)
{
	if (t < T1 || t > T2)
		return false;
	return true;
}


void BSWindowBarrierPricer::apply_barriers(
		double L,
		double U,
		const std::vector<double>& s,
		double T1,
		double T2,
		double t,
		std::vector<double>& v)
{
	if (!is_monitored(T1, T2, t))
		return;
	BSBarrierPricer::apply_barriers(L, U, s, v);
}


// =================================================================================================


void BSBarrierPricer::test()
{
	double S0 = 3190.93;
	double K = 3150;
	double r = 0.0;
	double q = 0.0;
	double sigmas[] = {0.2, 0.25, 0.32, 0.44};
	double T = 0.0328767;
	double Ls[] = {2850, 70, 60};
	double Us[] = {9120, 130, 140};
	int theta = -1;
	
	// convergence test
	std::cout << "\nconvergence test for barrier\n";
	std::cout << "============================\n";
	fprintf(stdout, "%10s %10s %10s %10s\n", "accuracy", "npv", "npv0", "err");
	double sigma = sigmas[0];
	double L = Ls[0];
	double U = Us[0];
	double npv0 = calculate(S0, K, r, q, sigma, T, theta, L, U, 16);
	
	int accuracy = 1;
#ifdef DEBUG
	for (size_t ii=0; ii<1; ii++) {
#else
	for (size_t ii=0; ii<4; ii++) {
#endif
		accuracy *= 2;
		double npv = calculate(S0, K, r, q, sigma, T, theta, L, U, accuracy);
		fprintf(stdout, "%10d %10.4f %10.4f %10.2e\n", accuracy, npv, npv0, npv - npv0);
	}
	
	// reproduce
	accuracy = 2;
	std::cout << "\npaper table for barrier with accuracy " << accuracy << "\n";
	std::cout << "=======================================\n";
	fprintf(stdout, "%7s %10s %10s %10s %10s\n", "L/U", "sigma", "npv", "npv0", "err");
	for (size_t ii=0; ii<3; ii++) {
		double L = Ls[ii];
		double U = Us[ii];
		for (size_t jj=0; jj<4; jj++) {
			double npv = calculate(S0, K, r, q, sigmas[jj], T, theta, L, U, accuracy);
			double npv0 = calculate(S0, K, r, q, sigmas[jj], T, theta, L, U, 16);
			fprintf(stdout, "%3.0f/%3.0f %10.2f %10.4f %10.4f %10.2e\n", 
					L, U, sigmas[jj], npv, npv0, npv - npv0);
		}
	}
	
}


double BSBarrierPricer::calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,
		double L,
		double U,
		int accuracy)
{
	// finite difference pricer for double barrier KO
	
	// space grids
	// either use L/U as bounds (better accuracy) or
	// similar to vanilla pricer but with concentration around L/U
	size_t ns = accuracy * std::min(130*1.5, 130+20*sigma*sigma*T);
	std::vector<double> s(ns);
	// idx is the index for S0
	size_t s0_idx;
	set_space_grids(S0, K, r, q, sigma, T, L, U, s, s0_idx);
	
	// time grids
	size_t nt = accuracy * std::min(50*1.5, 50+2*T);
	double dt = T / (nt - 1);
	
	// terminal conditions
	std::vector<double> vnext(ns, 0), v(ns, 0), rhs(ns, 0);
	for (size_t ii=0; ii<ns; ii++)
		vnext[ii] = std::max(0.0, theta * (s[ii] - K));
	apply_barriers(L, U, s, vnext);	
	
	// move backwards from nt-1 to nt-2 with FIM
	// helps reduce npv errors, e.g., from 0.01 to 0.001 n restore convergence rate from 1 to 2
	fim_step(r, q, sigma, L, U, dt, s, vnext, v, rhs, true);
	
#ifdef DEBUG
	std::cout << "vnext after FIM:\n";
	print(vnext);
#endif

	// move backwards from k+1 to k
	// CN tridiagonal matrix on lhs n rhs
	TridiagonalMatrix lhs_mat(ns), rhs_mat(ns);
	set_matrix_rhs_cn(s, r, q, sigma, dt, lhs_mat, rhs_mat);
	update_matrix_rhs(L, U, s, dt, lhs_mat, rhs_mat);
	
	for (int k=nt-3; k>=0; k--) {
		rhs_mat.product(vnext, rhs);
		lhs_mat.solve(rhs, v);
		v.swap(vnext);
	}
	v.swap(vnext);	
	
	return v[s0_idx];
	
}


void BSBarrierPricer::cn_step(
		double r,
		double q,
		double sigma,
		double L,
		double U,
		double dt,
		const std::vector<double>& s,
		std::vector<double>& vnext,
		std::vector<double>& v,
		std::vector<double>& rhs,
		bool is_monitored)
{
	// to be used for irregular dt
	size_t ns = s.size();
	TridiagonalMatrix lhs_mat(ns), rhs_mat(ns);
	set_matrix_rhs_cn(s, r, q, sigma, dt, lhs_mat, rhs_mat);
	
	if (is_monitored)
		update_matrix_rhs(L, U, s, dt, lhs_mat, rhs_mat);

	rhs_mat.product(vnext, rhs);
	lhs_mat.solve(rhs, v);
	v.swap(vnext);
	
}


void BSBarrierPricer::fim_step(
		double r,
		double q,
		double sigma,
		double L,
		double U,
		double dt,
		const std::vector<double>& s,
		std::vector<double>& vnext,
		std::vector<double>& v,
		std::vector<double>& rhs,
		bool is_monitored)
{
	const size_t num_fim_steps = 4;
	
	size_t ns = s.size();
	TridiagonalMatrix lhs_mat(ns), rhs_mat(ns);
	set_matrix_rhs_fim(s, r, q, sigma, dt/num_fim_steps, lhs_mat, rhs_mat);
	
	if (is_monitored)
		update_matrix_rhs(L, U, s, dt/num_fim_steps, lhs_mat, rhs_mat);

	for (size_t ii=0; ii<num_fim_steps; ii++) {
		rhs_mat.product(vnext, rhs);
		lhs_mat.solve(rhs, v);
		v.swap(vnext);
	}
	
}


void BSBarrierPricer::update_matrix_rhs(
		double L,
		double U,
		const std::vector<double>& s,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat)
{
	std::vector<double>& ll = lhs_mat.get_a();
	std::vector<double>& ld = lhs_mat.get_b();
	std::vector<double>& lu = lhs_mat.get_c();

	std::vector<double>& rl = rhs_mat.get_a();
	std::vector<double>& rd = rhs_mat.get_b();
	std::vector<double>& ru = rhs_mat.get_c();
	
	size_t ns = s.size();
	for (size_t ii=0; ii<ns; ii++) {
		if (s[ii] > L && s[ii] < U)
			continue;
		ll[ii] = 0;
		ld[ii] = 1 / dt;
		lu[ii] = 0;

		rl[ii] = 0;
		rd[ii] = 0;
		ru[ii] = 0;		
	}

}


void BSBarrierPricer::apply_barriers(
		double L,
		double U,
		const std::vector<double>& s,
		std::vector<double>& v)
{
	for (size_t ii=0; ii<s.size(); ii++) {
		if (s[ii] > L && s[ii] < U)
			continue;
		v[ii] = 0;
	}
}


void BSBarrierPricer::set_space_grids(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double L,
		double U,
		std::vector<double>& s,
		size_t& idx)
{
	double smin, smax;
	set_space_range(S0, K, r, q, sigma, T, smin, smax);
	smin = std::min(0.8*L, smin);
	smax = std::max(1.2*U, smax);
	
	std::vector<double> critical_points(3);
	critical_points[0] = K;
	critical_points[1] = L;
	critical_points[2] = U;
	
	FDGridHelper::set_concentrated_grids(smin, smax, critical_points, s);
	
	// shit so that S0, L, U are on grid
	// critical to achieve 2nd-order convergence
	idx = adjust_grids(S0, s);	
	size_t l_idx = adjust_grids(L, s);
	size_t u_idx = adjust_grids(U, s);
	
#ifdef DEBUG
	std::cout << "smin " << smin << " smax " << smax << " ns " << s.size() << "\n";
	print(s);
#endif
}


// =================================================================================================


void FDGridHelper::test()
{
	std::cout << "\ntesting FDGridHelper\n";
	std::cout << "====================\n";
	
	double smin = 50;
	double smax = 200;
	std::vector<double> critical_points(3), vec_space_grid(20);
	critical_points[0] = 100;
	critical_points[1] = 80;
	critical_points[2] = 120;
	
	set_concentrated_grids(smin, smax, critical_points, vec_space_grid);
	print(vec_space_grid);
	
}


void  FDGridHelper::set_concentrated_grids(
		double smin,
		double smax,
		const std::vector<double>& critical_points,
		std::vector<double>& vec_space_grid
		)
{
     /* critical points related */
     double alpha = 2.0;

     double A = 10.0;
	 double high;
     int rcode = upper_bound_for_ruger_kutta(
                 smin,
                 smax,
                 alpha,
                 critical_points,
                 vec_space_grid,
                 high
                 );
     if (rcode != 0){
         std::cout << "error from upper_bound_for_ruger_kutta\n";
         return;
     }

     /*
      * do the ruger-kutta integrator first to set up the initial space grids
      * then use bisection method to make the last space grid match smax
      */

     int count_loop = 0;
     double low = 0.0;
	 int is_success = 0;
	 const int max_num_loop = 50;

     const int num_space_grid = vec_space_grid.size();
     do{
         ruger_kutta(
                     0,
                     smin,
                     1,
                     A,
                     alpha,
                     critical_points,
                     vec_space_grid
                     );

         if (fabs(vec_space_grid[num_space_grid-1] - smax) <= eps ||
             fabs(vec_space_grid[num_space_grid-1] - smax)/smax <= eps){

             is_success = 1;

         } else {

             if ( vec_space_grid[num_space_grid-1] <= smax){
                 low = A;
                 A = (A + high) / 2;
             } else {
                 high = A;
                 A = (A + low) /2;
             }

             ++count_loop;
         }

     } while(!is_success && count_loop < max_num_loop);

     /* set up boundary space grids */
     vec_space_grid[0] = smin;
     vec_space_grid[num_space_grid-1] = smax;
}


double FDGridHelper::global_jacobian(
        const double y,
        const double alpha,
        const std::vector<double>& critical_points
        )
{
	const int num_points = critical_points.size();
    double result = 0.0;
    double temp;
    int ii;

    for(ii=0; ii<num_points; ++ii){
        temp = alpha * alpha + (y - critical_points[ii]) * (y - critical_points[ii]);
        result += 1 / temp;
    }

    result =  1.0 / sqrt( result );

    return result;
}


/**\brief Ruger-Kutta integrator
 *
 * Integrate dS(x)/dx = f(x,S(x))
 *
 * \param[in] (*f)              : derivative function to be integrated.
 * \param[in] x_start           : the start point of the integration.
 * \param[in] y_start           : the initial value of S(x0).
 * \param[in] xn                : the end point of the integration.
 * \param[in] A                 : constant to be determined.
 * \param[in] alpha             : parameter of the function, which indicates the concentration.
 * \param[in] critical points   : the critical points array,
 *                                local minimal points of the Jacobian function.
 * \param[in] num_points        : number of critical points array
 * \param[in] num               : number of integration steps.
 * \param[out] result           : output array of the integrated path S(x).
 *
 */
void FDGridHelper::ruger_kutta(
        const double x_start,
        const double y_start,
        const double xn,
        const double A,
        const double alpha,
        const std::vector<double>& critical_points,
        std::vector<double>& result
        )
{
	const int num = result.size();
	const int num_points = critical_points.size();
    std::vector<double> x(num);
    double h;
    double k1;
    double k2;
    double k3;
    double k4;
    int ii;

    x[0] = x_start;
    result[0] = y_start;
    h = (xn - x_start) / (num - 1);

    for (ii = 0; ii < num-1; ii++){
    
        k1 = A * global_jacobian(result[ii], alpha, critical_points);
        k2 = A * global_jacobian(result[ii] + h / 2 * k1, alpha, critical_points);
        k3 = A * global_jacobian(result[ii] + h / 2 * k2, alpha, critical_points);
        k4 = A * global_jacobian(result[ii] + h * k3, alpha, critical_points);

        result[ii+1] = result[ii] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        x[ii+1] = x[ii] + h;
    }
}


/**\brief get an upper bound for ruger-kutta integrator
*/
int FDGridHelper::upper_bound_for_ruger_kutta(
        double smin,
        double smax,
        double alpha,
        const std::vector<double>& critical_points,
        std::vector<double>& vec_space_grid,
        double& high
        )
{
    const int max_num_loop = 8;
    int count_loop = 0;
    int is_found = 0;
    double upper_bound = 20.0;
	const int num_space_grid = vec_space_grid.size();
	
    /* if the last space grid didn't exceed smax, then double high */
    while (!is_found && count_loop < max_num_loop){

        ruger_kutta(
                    0,
                    smin,
                    1,
                    upper_bound,
                    alpha,
                    critical_points,
                    vec_space_grid
                    );

        if (vec_space_grid[num_space_grid - 1] >= smax){
            is_found = 1;
        } else {
            upper_bound *= 2;
            ++count_loop;
        }
    }

    if (!is_found) {
        high = 0.0;
        return -1;
    }
    
    high = upper_bound;
    return 0;
}