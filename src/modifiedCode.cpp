#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins("cpp11")]]
//using namespace Rcpp;
//using namespace arma;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//#include <vector.h>
//#include <array.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

typedef Rcpp::NumericVector NumericVector;
typedef arma::vec vect;
typedef std::vector<vect> twoVect;
typedef arma::mat matrix;
typedef std::vector<arma::mat > array3d;

/*
void freemem(int ndim, ...);
void ifreemem(int ndim, ...);
int    fwhich_min(double *x, int n);
double inprod(double *v1, double *v2, int d);
double int_simp(double *f, double h, int n);
double int_romberg(double *f);
double int_trape(double *f);
double ldexgam(double x,double shape,double scale);
double norm2(double *vec, int len);
double pexgam(double x,double shape,double scale);
double qexgam(double p,double shape,double scale);
void quadsolve(double a, double b, double c, double *x, int get_both);
int    rdisc(int n,double *prob,int islog);
void r_simple_dirichlet(double *x, int k, double *shape);
double rexgam(double shape, double scale);
void   rperm(int size, int *samp, int n);
void   Rprintvec(char *a, double *x, int len);
void   Rprintveci(char *a, int *x, int len);
void   Rprintmat(char *a, double **m, int d1, int d2, int iscol);
void   Rprintmati(char *a, int **m, int d1, int d2, int iscol);
double sigmoid(double x);
double sigmoid2(double x);
double sum(double *x, int n);
double sumlog(double *x, int n);
void swap(double *x, int i, int j);
void transpose(double **A, int m, int n, double **B);
void Uprod(double **A, int n, double *x, double *b, int transpose);
double trace(double **a, int n, int is_log);
double fracf(double x);
double fold(double x);
double logadd(double lx, double ly);
*/


int order, *twopow, grid_len;
double grid_incr;

int nobs, nnode, ngrid, nsweep, nvar, *nrep, gAnchor;
double incr, lincr, nuL, muL, nuBp, muBp, nuBa, muBa;
double muBase, *betBase, sigmaBase, sigmaBaseSq, nuS = 2.0;

//classic functions :
double ss(const vect& myVec){ //sum of square
  return sum(myVec*myVec);
}

//--- fold a number to (0,1) by boundary reflection ---
double fold(double x){
  double ix;
	modf(x, &ix);
	x = fabs(x - ix);
	if((int)ix % 2)
		x = 1.0 - x;
	return x;
}

//Linear algebra functions
vect cprod(const matrix& A,const vect& b){
  return A.t()*b;
}
vect Usolve(const matrix& A,const vect& b, int transpose){
  if(transpose) return arma::solve(A.t(),b);
  return arma::solve(A,b);
}

//auxilliary functions
void get_proj(matrix& x, vect& proj, int nobs, int nvar, vect& px){
  
	int i;
	//cprod(x, proj, nvar, nobs, px);	
  px = cprod(x,proj);
//	double pmin = f_min(px, nobs);
//	double pmax = f_max(px, nobs);
	double pmin = px.min();
	double pmax = px.max();
	double pmid = 0.5 * (pmin + pmax);
	double pwid = (pmax - pmin);
	for(i = 0; i < nobs; i++)
		px[i] = (px[i] - pmid) / pwid;
	px[i++] = pmid;
	px[i++] = pwid;
}


double bigQ(double x){
  double a = R::qt(x, 1.0, 1, 0);
  if(a < -1.0e10) a = -1.0e10;
	if(a > 1.0e10) a = 1.0e10;
	return a;
}

void lpost(matrix& xMat, vect& y,
            vect& pnodes, vect& anodes,
    			 double mu0, vect& bet0,
           double sigma0, double sigma1,
					 vect& b2, vect& omega,
           vect& proj, vect& px, 
					 matrix& cmat_node, matrix& cmat_grid,
           double& logdetR, vect& xi,
           vect& rxi, matrix& q,
           vect& p0x, double & loglik,
           double & logpost, int doMat,
           int doFun, int doProj, 
           int doBet0, int toPrint,
           //added
           matrix& eta,vect& eta_cond,
           matrix& qd,vect& pgrid
					 ){
	int g, k, k2, n, j, skip;
	double q_curr, q_last;
	double vp, va0, va1;
	
	if(doBet0){
    p0x = cprod(xMat,bet0);
	}
			
			
	if(doMat){
		for(k = 0; k < nnode; k++){
			for(k2 = 0; k2 < k; k2++){
				cmat_node[k,k2] = exp(- b2[0] * (pnodes[k] - pnodes[k2])*(pnodes[k] - pnodes[k2])
															 - b2[1] * (anodes[k] - anodes[k2])*(anodes[k] - anodes[k2]));
				cmat_node[k2,k] = cmat_node[k,k2];
			}
			cmat_node[k,k] = 1.0 + 1.0e-10;
		}
    cmat_node = chol(cmat_node); //cholesky factorization
		for(logdetR = 0.0, k = 0; k < nnode; k++)
			logdetR += log(cmat_node[k,k]);
		
		for(k = 0; k < nnode; k++){
			va0 = b2[1] * anodes[k]*anodes[k];
			va1 = b2[1] * (anodes[k] - 1.0)*(anodes[k] - 1.0);
			for(g = 0; g < ngrid; g++){
				vp = b2[0] * (pnodes[k] - pgrid[g])*(pnodes[k] - pgrid[g]);
				cmat_grid[k,g] = exp(- vp - va0);
				cmat_grid[k,ngrid + g] = exp(- vp - va1);
			}
		}
	}
	

	if(doFun){
    rxi = Usolve(cmat_node,omega,1);
		xi = Usolve(cmat_node, rxi,0);
		
		for(g = 0; g < ngrid; g++){
			for(eta[0,g] = 0.0, k = 0; k < nnode; k++)
				eta[0,g] += xi[k] * cmat_grid[k,g];
			
			for(eta[1,g] = 0.0, k = 0; k < nnode; k++)
				eta[1,g] += xi[k] * cmat_grid[k,ngrid + g];
		}
		
		eta_cond[0] = max(eta.row(0));
		eta_cond[1] = max(eta.row(1));
		
		for(g = 0; g < ngrid; g++){
			qd[0,g] = exp(eta[0,g] - eta_cond[0]);
			qd[1,g] = exp(eta[1,g] - eta_cond[1]);
		}
		
		for(q[0,0] = 0.0, g = 1; g < ngrid; g++)
			q[0,g] = q[0,g - 1] + 0.5 * (qd[0,g - 1] + qd[0,g]);
		for(q[1,0] = 0.0, g = 1; g < ngrid; g++)
			q[1,g] = q[1,g - 1] + 0.5 * (qd[1,g - 1] + qd[1,g]);
		
		for(g = 0; g < ngrid; g++){
			q[0,g] /= q[0,ngrid - 1];
			q[1,g] /= q[1,ngrid - 1];
		}
		
	}
	
	if(doProj){
		get_proj(xMat, proj, nobs, nvar, px);
	}
			

	double alph, r;
	skip = 0;
	double bigQ_anchor[] = {0.0, 0.0};
	for(loglik = 0.0, n = 0; n < nobs; n++){
		g = 0;
		q_last = 0.0;
		alph = 0.5 - px[n];
		q_curr = (alph * sigma0 * (bigQ(q[0,g]) - bigQ_anchor[0]) 
							+ (1.0 - alph) * sigma1 * (bigQ(q[1,g]) - bigQ_anchor[1]));
		for(j = 0; j < nrep[n]; j++){
			r = y[skip + j] - mu0 - p0x[n];
			while(q_curr < r){
				g++;
				q_last = q_curr;
				q_curr = (alph * sigma0 * (bigQ(q[0,g]) - bigQ_anchor[0]) 
									+ (1.0 - alph) * sigma1 * (bigQ(q[1,g]) - bigQ_anchor[1]));
			}
			if(toPrint)
				Rprintf("x = %g, r = %g, g = %d, q_brack = [%g, %g]\n",
								px[n], r, g, q_last, q_curr);
			
			loglik += lincr - log(q_curr - q_last);
		}
		skip += nrep[n];
	}
	if(toPrint)
		Rprintf("skip = %d\n", skip);
	
	logpost = (loglik
								+ R::dgamma(b2[0], nuBp, muBp, 1)
								+ R::dgamma(b2[1], nuBa, muBa, 1)
								- logdetR
								- 0.5 * (nuL + nnode) * log(1.0 + ss(rxi) / muL)				  
								- 0.5 * (1.0 + nvar) * log(1.0 + ss(proj))
								);
}

void sweep(matrix& xMat, vect& y, vect& pnodes, vect& anodes,
					 array3d& cmat_node, array3d& cmat_grid, vect& tune,
					 vect& b2, vect& omega, vect& proj, vect& theta,
					 vect& st_b2, vect& st_xi, vect& st_omega, 
					 vect& st_proj, vect& st_px, vect& st_theta,
					 vect& st_loglik, vect& st_logpost, vect& st_logdetR, 
					 vect& acpt, int doRho, int doBase, vect& tuneBase, vect& acptBase,
           //added
           matrix& eta, vect& eta_cond,matrix& qd, vect& pgrid
					 ){
	
	int m, k, skipB = 0, skip = 0, skip_om = 0, skip_xi = 0, skip_proj = 0, skip_px = 0, skip_theta = 0;
  
  vect loglik(2);
	vect logpost(2);
	vect logdetR(2);
  
	twoVect xi(2,vect(nnode));
	twoVect rxi(2,vect(nnode));
	twoVect px(2,vect(nobs+2));
	twoVect p0x(2,vect(nobs));
  
	array3d q(2,matrix(2,ngrid)) ;

	double mu0 = theta[0];
  //double *bet0 = theta + 1;
  vect bet0(nvar);
  for(int j = 0; j < nvar; j++) bet0[j] = theta[j+1];
  double sigmaCom = theta[nvar + 1], rho = theta[nvar + 2];
	double sigma0 = sigmaCom * sqrt(2.0 * rho), sigma1 = sigmaCom * sqrt(2.0 * (1.0 - rho));
	
	int holdB = 0, propB = 1, holdBW = 0, propBW = 1, holdP = 0, propP = 1, hold = 0, prop = 1;
	double nacpt[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, b2_hold[2], omega_hold, u, proj_hold[nvar], rho_hold; 
	int hold_base = 0, prop_base = 1, move, bet0_pick; 
	double mu0_hold, sigma_hold[2], bet0_hold[nvar], bet0_prop, nacptBase[3] = {0.0, 0.0, 0.0};
	
	int k1, k2;
	double tot_norm, 	pi = atan(1.0) * 4.0;
	
	lpost(xMat, y, 
	      pnodes, anodes,
				//						mu0, bet0, sigma0, rho,
				mu0, bet0, sigma0, sigma1,
				b2, omega, proj,
				px[holdP], cmat_node[holdB], cmat_grid[holdB], logdetR[holdB],
				xi[holdBW], rxi[holdBW], q[holdBW], p0x[hold_base],
				loglik[hold], logpost[hold],
				1, 1, 1, 1, 0,eta,eta_cond,qd,pgrid);
	
	Rprintf("Initialized at: loglik = %g\tlogpost = %g\n", loglik[hold], logpost[hold]);
	double u0, u1;
	
	for(m = 0; m < nsweep; m++){
		
		
		// ++++++++++ Update omega ++++++++++++++ 
		
		k = floor(R::runif(0.0, (double) nnode));
		nacpt[0] += 1.0;
		propBW = !holdBW;
		prop = !hold;
		omega_hold = omega[k];
		omega[k] += tune[0] * static_cast<double>(R::rnorm(0.0, 1.0));
		
		lpost(xMat, y, 
					pnodes, anodes,
					mu0, bet0, sigma0, sigma1,
					b2, omega, proj, 
					px[holdP], cmat_node[holdB], cmat_grid[holdB], logdetR[holdB],
					xi[propBW], rxi[propBW], q[propBW], p0x[hold_base],
					loglik[prop], logpost[prop],
					0, 1, 0, 0, 0,eta,eta_cond,qd,pgrid);
		
		if(log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold])){
			holdBW = propBW; 
			hold = prop;
			acpt[0] += 1.0;
		} else {
			omega[k] = omega_hold;
		}
		
		
		// ++++++++++ Update b2 ++++++++++++++ 
		
		nacpt[1] += 1.0;
		propB = !holdB;
		propBW = !holdBW;
		prop  = !hold;		
		b2_hold[0] = b2[0];
		u = R::rnorm(0.0, tune[1]);
		b2[0] *= exp(u);
		
		lpost(xMat, y, 
					pnodes, anodes,
					//						mu0, bet0, sigma0, rho,
					mu0, bet0, sigma0, sigma1,
					b2, omega, proj,
					px[holdP], cmat_node[propB], cmat_grid[propB], logdetR[propB],
					xi[propBW], rxi[propBW], q[propBW], p0x[hold_base],
					loglik[prop], logpost[prop],
					1, 1, 0, 0, 0,eta,eta_cond,qd,pgrid);
		
		if(log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold] + u)){
			holdB = propB;
			holdBW = propBW; 
			hold  = prop;
			acpt[1] += 1.0;
		} else {
			b2[0] = b2_hold[0];
		}
		
		nacpt[2] += 1.0;
		propB = !holdB;
		propBW = !holdBW;
		prop  = !hold;		
		b2_hold[1] = b2[1];
		u = R::rnorm(0.0, tune[2]);
		b2[1] *= exp(u);
		
		lpost(xMat, y, 
					pnodes, anodes,
					mu0, bet0, sigma0, sigma1,
					b2, omega, proj,
					px[holdP], cmat_node[propB], cmat_grid[propB], logdetR[propB],
					xi[propBW], rxi[propBW], q[propBW], p0x[hold_base],
					loglik[prop], logpost[prop],
					1, 1, 0, 0, 0,eta,eta_cond,qd,pgrid);
		
		if(log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold] + u)){
			holdB = propB;
			holdBW = propBW; 
			hold  = prop;
			acpt[2] += 1.0;
		} else {
			b2[1] = b2_hold[1];
		}
		
		
		// ++++++++++ Update proj ++++++++++++++ 
		
		if(nvar > 1){
			nacpt[3] += 1.0;
			propP = !holdP;
			prop  = !hold;
			for(k = 0; k < nvar; k++){
				proj_hold[k] = proj[k];
				proj[k] += tune[3] * R::rnorm(0.0, 1.0) * R::rbeta(0.5, 1.0);	
			}
			
			lpost(xMat, y, 
						pnodes, anodes,
						mu0, bet0, sigma0, sigma1,
						b2, omega, proj,
						px[propP], cmat_node[holdB], cmat_grid[holdB], logdetR[holdB],
						xi[holdBW], rxi[holdBW], q[holdBW], p0x[hold_base],
						loglik[prop], logpost[prop],
						0, 0, 1, 0, 0,eta,eta_cond,qd,pgrid);
			
			if(log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold])){
				holdP = propP;
				hold  = prop;
				acpt[3] += 1.0;
			} else {
				for(k = 0; k < nvar; k++)
					proj[k] = proj_hold[k];
			}
		}		
		
		
		// ++++++++++ Update rho +++++++++++
		
		if(doRho){
			nacpt[4] += 1.0;
			prop  = !hold;		
			rho_hold = rho;
			rho = fold(rho + R::rnorm(0.0, tune[4]));
			u = R::dbeta(rho, nuS, nuS, 1) - R::dbeta(rho_hold, nuS, nuS, 1);
			lpost(xMat, y, 
						pnodes, anodes,
						//						mu0, bet0, sigma0, rho,
						mu0, bet0, sigma0, sigma1,
						b2, omega, proj,
						px[holdP], cmat_node[holdB], cmat_grid[holdB], logdetR[holdB],
						xi[holdBW], rxi[holdBW], q[holdBW], p0x[hold_base],
						loglik[prop], logpost[prop],
						0, 0, 0, 0, 0,eta,eta_cond,qd,pgrid);
			
			if( log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold] + u) ){
				hold  = prop;
				acpt[4] += 1.0;
			} else {
				rho = rho_hold;
			}
		}
		

		// ++++++++ Update base parameters mu0, bet0, sigma  ++++ 
		
		if(doBase){
			move = floor(R::runif(0.0, 3.0));			
			
			if(move == 0){
				mu0_hold = mu0;
				mu0 += R::rnorm(0.0, tuneBase[0]);
				u = R::dt(3.0 * (mu0 - muBase) / sigmaBase, 1.0, 1) - R::dt(3.0 * (mu0_hold - muBase) / sigmaBase, 1.0, 1);
			} else if(move == 1){
				u = 0.0;
				for(k = 0; k < nvar; k++){				
					bet0_hold[k] = bet0[k]; 
					bet0[k] += tuneBase[1] * R::rnorm(0.0, 1.0) * R::rbeta(1.0, 2.0);			
					u += R::dt(3.0 * (bet0[k] - betBase[k]) / sigmaBase, 1.0, 1) - R::dt(3.0 * (bet0_hold[k] - betBase[k]) / sigmaBase, 1.0, 1);
				}
			} else {
				sigma_hold[0] = sigma0;
				sigma_hold[1] = sigma1;
				sigma0 *= exp(R::rnorm(0.0, tuneBase[2]));
				sigma1 *= exp(R::rnorm(0.0, tuneBase[2]));
											
				u = 0.0;
				u += (R::dgamma(sigma0*sigma0, nuS, sigmaBaseSq / nuS, 1) + 2.0 * log(sigma0)
							- R::dgamma(sigma_hold[0]*sigma_hold[0], nuS, sigmaBaseSq / nuS, 1) - 2.0 * log(sigma_hold[0])
							+ R::dgamma(sigma1*sigma1, nuS, sigmaBaseSq / nuS, 1) + 2.0 * log(sigma1)
							- R::dgamma(sigma_hold[1]*sigma_hold[1], nuS, sigmaBaseSq / nuS, 1) - 2.0 * log(sigma_hold[1])
							);
			}
			prop = !hold;
			nacptBase[move] += 1.0;
			prop_base = (move == 1) ? !hold_base : hold_base;
			
			lpost(xMat, y, 
						pnodes, anodes,
						mu0, bet0, sigma0, sigma1,
						b2, omega, proj,
						px[holdP], cmat_node[holdB], cmat_grid[holdB], logdetR[holdB],
						xi[holdBW], rxi[holdBW], q[holdBW], p0x[prop_base],
						loglik[prop], logpost[prop],
						0, 0, 0, (move == 1), 0,eta,eta_cond,qd,pgrid);
			
			if(log(R::runif(0.0, 1.0)) < (logpost[prop] - logpost[hold] + u)){
				hold = prop;
				hold_base = prop_base;
				acptBase[move] += 1.0;
			} else {
				if(move == 0){
					mu0 = mu0_hold;
				} else if(move == 1){
					for(k = 0; k < nvar; k++){
						bet0[k] = bet0_hold[k];
					}
				} else {
					sigma0 = sigma_hold[0];
					sigma1 = sigma_hold[1];
				}
			}
		}
		
		// +++++++++++++ Store ++++++++++++++
		
		st_b2[skipB++] = b2[0];
		st_b2[skipB++] = b2[1];
		for(k = 0; k < nnode; k++)
			st_xi[skip_xi++] = xi[holdBW][k];	
		for(k = 0; k < nnode; k++)
			st_omega[skip_om++] = omega[k];
		for(k = 0; k < nvar; k++)
			st_proj[skip_proj++] = proj[k];
		for(k = 0; k < 2; k++)
			st_px[skip_px++] = px[holdP][nobs + k];		
		st_theta[skip_theta++] = mu0;
		for(k = 0; k < nvar; k++)
			st_theta[skip_theta++] = bet0[k];
		st_theta[skip_theta++] = sqrt(0.5 * (sigma0 * sigma0 + sigma1 * sigma1));
		st_theta[skip_theta++] = (sigma0 * sigma0) / (sigma0 * sigma0 + sigma1 * sigma1);
		st_loglik[m] = loglik[hold];
		st_logpost[m] = logpost[hold];	
		st_logdetR[m] = logdetR[holdB];
		
		
		// +++++++++++ Print progress +++++++++++++
		
		if((10 * (m + 1)) % nsweep == 0)
			Rprintf("iter = %d, loglik = %5g logpost = %5g\n", 
							m + 1, loglik[hold], logpost[hold]);
		
	}
  
	acpt[0] /= nacpt[0] + 1.0e-10;
	acpt[1] /= nacpt[1] + 1.0e-10;
	acpt[2] /= nacpt[2] + 1.0e-10;
	acpt[3] /= nacpt[3] + 1.0e-10;
	acpt[4] /= nacpt[4] + 1.0e-10;
	acptBase[0] /= nacptBase[0] + 1.0e-10;
	acptBase[1] /= nacptBase[1] + 1.0e-10;
	acptBase[2] /= nacptBase[2] + 1.0e-10;
}

void mcmc(vect& x, vect& y, vect& pnodes, vect& anodes,
					vect& tune, vect& b2, vect& omega, vect& proj, vect& theta,
					vect& st_b2, vect& st_xi, vect& st_omega, 
					vect& st_proj, vect& st_px, vect& st_theta,
					vect& st_loglik, vect& st_logpost, vect& st_logdetR,
					vect& acpt, int doRho, int doBase, vect& tuneBase, vect& acptBase
					){

	int k, k2, g, i;
	double vp, va0, va1;
	
	//double **xMat = setmem(2, nvar, nobs);
  matrix xMat(2,nvar);
	int skip_x = 0;
	for(k = 0; k < nvar; k++){
		for(i = 0; i < nobs; i++){
			xMat[k,i] = x[skip_x];
			skip_x += nrep[i];
		}
	}

	//pgrid = setmem(1, ngrid);
  vect pgrid(ngrid,0.0);
	for(g = 1; g < ngrid; g++)
		pgrid[g] = pgrid[g - 1] + incr;
	
  array3d cmat_grid(2,matrix(nnode,2*ngrid));
  array3d cmat_node(2,matrix(nnode,nnode));
  matrix eta(2,ngrid);
  vect eta_cond(2);
  matrix qd(2,ngrid);

	GetRNGstate();
	sweep(xMat, y, 
	      pnodes, anodes,
	      cmat_node, cmat_grid, 
				tune, 
				b2, omega, proj, theta,
				st_b2, st_xi, st_omega, st_proj, st_px, st_theta,
				st_loglik, st_logpost, st_logdetR,
				acpt, doRho, doBase, tuneBase, acptBase,
        eta,eta_cond,qd,pgrid);
	PutRNGstate();
}

void initialize_ljqrf(int * dims,
                      NumericVector& hpar
  										){

	nobs = dims[0];
	nnode = dims[1];
	ngrid = dims[2];
	nsweep = dims[3];
	nvar = dims[4];
	nrep = dims + 5;
	gAnchor = ngrid / 2;
	
	incr = 1.0 / ((double)ngrid - 1.0);
	lincr = log(incr);
	
	nuL = hpar[0];
	muL = hpar[0] * hpar[1];
	nuBp = hpar[2];
	muBp = hpar[3];
	nuBa = hpar[4];
	muBa = hpar[5];
	muBase = hpar[6];
	betBase = hpar + 7;
	sigmaBase = hpar[7 + nvar];
	sigmaBaseSq = sigmaBase * sigmaBase;
	
}

// [[Rcpp::export]]
void slqr(NumericVector x, NumericVector y, int * dims, NumericVector pnodes, NumericVector anodes,
      		NumericVector hpar, NumericVector tune, NumericVector theta,
					NumericVector b2, NumericVector omega, NumericVector proj, 
					NumericVector st_b2, NumericVector st_xi, NumericVector st_omega, 
					NumericVector st_proj, NumericVector st_px, NumericVector st_theta,
					NumericVector st_loglik, NumericVector st_logpost, NumericVector st_logdetR,
					NumericVector acpt, int *doRho, int *doBase, NumericVector tuneBase, NumericVector acptBase
				 ){
	initialize_ljqrf(dims, hpar);
	Rprintf("Parameters: nuL = %g, muL = %g, nuBp = %g, muBp = %g, nuBa = %g, muBa = %g\n",
					nuL, muL, nuBp, muBp, nuBa, muBa);
	mcmc(x, y, pnodes, anodes, tune, b2, omega, proj, theta,
			 st_b2, st_xi, st_omega, st_proj, st_px, st_theta,
			 st_loglik, st_logpost, st_logdetR,
			 acpt, doRho[0], doBase[0], tuneBase, acptBase);
	
}

/*
void getq(vect& pnodes,
					vect& anodes,
					vect& b2,
					vect& xi,
					vect& q0,
					vect& q1,
          matrix& eta,vect& eta_cond,matrix& q,matrix& qd,vect& pgrid,matrix& cmat_grid
					){
	
	int g, k, n;
	double vp, va0, va1;
	
	for(k = 0; k < nnode; k++){
		va0 = b2[1] * anodes[k]*anodes[k];
		va1 = b2[1] * (anodes[k] - 1.0)*(anodes[k] - 1.0);
		for(g = 0; g < ngrid; g++){
			vp = b2[0] * (pnodes[k] - pgrid[g])*(pnodes[k] - pgrid[g]);
			cmat_grid[k,g] = exp(- vp - va0);
			cmat_grid[k,ngrid + g] = exp(- vp - va1);
		}
	}
	
	
	for(g = 0; g < ngrid; g++){
		for(eta[0,g] = 0.0, k = 0; k < nnode; k++)
			eta[0,g] += xi[k] * cmat_grid[k,g];
		
		for(eta[1,g] = 0.0, k = 0; k < nnode; k++)
			eta[1,g] += xi[k] * cmat_grid[k,ngrid + g];
	}
	
	eta_cond[0] = max(eta.row(0));
	eta_cond[1] = max(eta.row(1));
	
	for(g = 0; g < ngrid; g++){
		qd[0,g] = exp(eta[0,g] - eta_cond[0]);
		qd[1,g] = exp(eta[1,g] - eta_cond[1]);
	}
	
	for(q[0,0] = 0.0, g = 1; g < ngrid; g++)
		q[0,g] = q[0,g - 1] + 0.5 * (qd[0,g - 1] + qd[0,g]);
	for(q[1,0] = 0.0, g = 1; g < ngrid; g++)
		q[1,g] = q[1,g - 1] + 0.5 * (qd[1,g - 1] + qd[1,g]);
	
	for(g = 0; g < ngrid; g++){
		q0[g] = q[0,g] / q[0,ngrid - 1];
		q1[g] = q[1,g] / q[1,ngrid - 1];
	}
	
}

void fit(int *dims,
				 vect& pnodes,
				 vect& anodes,
				 vect& b2,
				 vect& xi,
				 vect& st_q,
         //
         matrix& cmat_grid,vect& pgrid
				 ){
	
	nnode = dims[0];
	ngrid = dims[1];
	nsweep = dims[2];
		
	incr = 1.0 / ((double)ngrid - 1.0);
	
	int k, k2, g, m;
	
	//pgrid = (double *)calloc(ngrid, sizeof(double));
	//for(pgrid[0] = 0.0, g = 1; g < ngrid; g++)
		//pgrid[g] = pgrid[g - 1] + incr;
	//int twongrid = 2 * ngrid;
	
	//cmat_grid = (double **)calloc(nnode, sizeof(double *));
	//for(k = 0; k < nnode; k++)
		//cmat_grid[k] = (double *)calloc(twongrid, sizeof(double));

   matrix eta(2,ngrid);
   vect eta_cond(2);
   matrix qd(2,ngrid);
   matrix q(2,ngrid);
   
	int b2_skip = 0, xi_skip = 0, q0_skip = 0, q1_skip = ngrid;
	for(m = 0; m < nsweep; m++){
		getq(pnodes, anodes, b2 + b2_skip, xi + xi_skip, st_q + q0_skip, st_q + q1_skip,eta,eta_cond,q,qd,cmat_grid);
		b2_skip += 2;
		xi_skip += nnode;
		q0_skip += twongrid;
		q1_skip += twongrid;
	}
}
*/