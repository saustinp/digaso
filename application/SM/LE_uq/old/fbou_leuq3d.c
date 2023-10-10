void fbou_leuq3d(double *fh, double *fh_udg, double *fh_uh, double *pg, double *udg, double *uh, 
        double *nl, double *ui, double *param, double time, int ib, int ng, int nc, int ncu, int nd, int ncd)
{
	double mu = param[0];
	double lambda = param[1];
	double tau = param[2];
    double uinf1 = ui[0];
    double uinf2 = ui[1];
    double uinf3 = ui[2];
    double t2, t3, t4, t5;
    
	for (int i = 0; i <ng; i++) {
		double x1 = pg[0*ng+i];
		double x2 = pg[1*ng+i];
		double x3 = pg[2*ng+i];
		double u1 = udg[0*ng+i];
		double u2 = udg[1*ng+i];
		double u3 = udg[2*ng+i];
		double q11 = udg[3*ng+i];
		double q21 = udg[4*ng+i];
		double q31 = udg[5*ng+i];
		double q12 = udg[6*ng+i];
		double q22 = udg[7*ng+i];
		double q32 = udg[8*ng+i];
		double q13 = udg[9*ng+i];
		double q23 = udg[10*ng+i];
		double q33 = udg[11*ng+i];
		double uh1 = uh[0*ng+i];
		double uh2 = uh[1*ng+i];
		double uh3 = uh[2*ng+i];
		double nl1 = nl[0*ng+i];
		double nl2 = nl[1*ng+i];
		double nl3 = nl[2*ng+i];
        
        
    }    
}

from numpy import *

def fbou(ib,uinf,nl,pg,udg,uh,param,time):

	ng = pg.shape[0]
	nd = pg.shape[1]
	nc = udg.shape[1]
	ncu = 3
	fh = zeros((ng,ncu))
	fh_udg = zeros((ng,ncu,nc))
	fh_uh = zeros((ng,ncu,ncu))

	x1 = pg[:,0]
	x2 = pg[:,1]
	x3 = pg[:,2]
	nl1 = nl[:,0]
	nl2 = nl[:,1]
	nl3 = nl[:,2]
	uinf1 = uinf[:,0]
	uinf2 = uinf[:,1]
	uinf3 = uinf[:,2]
	u1 = udg[:,0]
	u2 = udg[:,1]
	u3 = udg[:,2]
	u4 = udg[:,3]
	u5 = udg[:,4]
	u6 = udg[:,5]
	u7 = udg[:,6]
	u8 = udg[:,7]
	u9 = udg[:,8]
	u10 = udg[:,9]
	u11 = udg[:,10]
	u12 = udg[:,11]
	uh1 = uh[:,0]
	uh2 = uh[:,1]
	uh3 = uh[:,2]
	param1 = param[0]
	param2 = param[1]
	param3 = param[2]

	if ib ==1:
		  fh[:,0,0] = -uh1+uinf1+x1
		  fh[:,1,0] = -uh2+uinf2+x2
		  fh[:,2,0] = -uh3+uinf3+x3


		  fh_uh[:,0,0] = -1.0
		  fh_uh[:,1,1] = -1.0
		  fh_uh[:,2,2] = -1.0

	elif ib ==2:
		  t2 = q11+q22+q33+3.0
		  t3 = lambda*t2
		  t4 = q12+q21
		  t5 = q13+q31
		  t6 = q23+q32
		  fh[:,0,0] = -uinf1+nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*t4+mu*nl3*t5
		  fh[:,1,0] = -uinf2+nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*t4+mu*nl3*t6
		  fh[:,2,0] = -uinf3+nl3*(t3+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl1*t5+mu*nl2*t6

		  t2 = mu*nl2
		  t3 = mu*nl3
		  t4 = lambda*nl1
		  t5 = mu*nl1
		  t6 = mu*2.0
		  t7 = lambda+t6
		  t8 = lambda*nl2
		  t9 = lambda*nl3
		  fh_udg[:,0,0] = tau
		  fh_udg[:,0,3] = nl1*t7
		  fh_udg[:,0,4] = t2
		  fh_udg[:,0,5] = t3
		  fh_udg[:,0,6] = t2
		  fh_udg[:,0,7] = t4
		  fh_udg[:,0,9] = t3
		  fh_udg[:,0,11] = t4
		  fh_udg[:,1,1] = tau
		  fh_udg[:,1,3] = t8
		  fh_udg[:,1,4] = t5
		  fh_udg[:,1,6] = t5
		  fh_udg[:,1,7] = nl2*t7
		  fh_udg[:,1,8] = t3
		  fh_udg[:,1,10] = t3
		  fh_udg[:,1,11] = t8
		  fh_udg[:,2,2] = tau
		  fh_udg[:,2,3] = t9
		  fh_udg[:,2,5] = t5
		  fh_udg[:,2,7] = t9
		  fh_udg[:,2,8] = t2
		  fh_udg[:,2,9] = t5
		  fh_udg[:,2,10] = t2
		  fh_udg[:,2,11] = nl3*t7

		  fh_uh[:,0,0] = -tau
		  fh_uh[:,1,1] = -tau
		  fh_uh[:,2,2] = -tau

	elif ib ==3:
		  t2 = q11+q22+q33+3.0
		  t3 = lambda*t2
		  t4 = q23+q32
		  fh[:,0,0] = -uh1+uinf1+x1
		  fh[:,1,0] = -uinf2+nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*(q12+q21)+mu*nl3*t4
		  fh[:,2,0] = -uinf3+nl3*(t3+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl1*(q13+q31)+mu*nl2*t4

		  t2 = mu*nl1
		  t3 = mu*nl3
		  t4 = lambda*nl2
		  t5 = lambda*nl3
		  t6 = mu*nl2
		  t7 = mu*2.0
		  t8 = lambda+t7
		  fh_udg[:,1,1] = tau
		  fh_udg[:,1,3] = t4
		  fh_udg[:,1,4] = t2
		  fh_udg[:,1,6] = t2
		  fh_udg[:,1,7] = nl2*t8
		  fh_udg[:,1,8] = t3
		  fh_udg[:,1,10] = t3
		  fh_udg[:,1,11] = t4
		  fh_udg[:,2,2] = tau
		  fh_udg[:,2,3] = t5
		  fh_udg[:,2,5] = t2
		  fh_udg[:,2,7] = t5
		  fh_udg[:,2,8] = t6
		  fh_udg[:,2,9] = t2
		  fh_udg[:,2,10] = t6
		  fh_udg[:,2,11] = nl3*t8

		  fh_uh[:,0,0] = -1.0
		  fh_uh[:,1,1] = -tau
		  fh_uh[:,2,2] = -tau

	elif ib ==4:
		  t2 = q11+q22+q33+3.0
		  t3 = lambda*t2
		  t4 = q13+q31
		  fh[:,0,0] = -uinf1+nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*(q12+q21)+mu*nl3*t4
		  fh[:,1,0] = -uh2+uinf2+x2
		  fh[:,2,0] = -uinf3+nl3*(t3+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl2*(q23+q32)+mu*nl1*t4

		  t2 = mu*nl2
		  t3 = mu*nl3
		  t4 = lambda*nl1
		  t5 = lambda*nl3
		  t6 = mu*nl1
		  t7 = mu*2.0
		  t8 = lambda+t7
		  fh_udg[:,0,0] = tau
		  fh_udg[:,0,3] = nl1*t8
		  fh_udg[:,0,4] = t2
		  fh_udg[:,0,5] = t3
		  fh_udg[:,0,6] = t2
		  fh_udg[:,0,7] = t4
		  fh_udg[:,0,9] = t3
		  fh_udg[:,0,11] = t4
		  fh_udg[:,2,2] = tau
		  fh_udg[:,2,3] = t5
		  fh_udg[:,2,5] = t6
		  fh_udg[:,2,7] = t5
		  fh_udg[:,2,8] = t2
		  fh_udg[:,2,9] = t6
		  fh_udg[:,2,10] = t2
		  fh_udg[:,2,11] = nl3*t8

		  fh_uh[:,0,0] = -tau
		  fh_uh[:,1,1] = -1.0
		  fh_uh[:,2,2] = -tau

	elif ib ==5:
		  t2 = q11+q22+q33+3.0
		  t3 = lambda*t2
		  t4 = q12+q21
		  fh[:,0,0] = -uinf1+nl1*(t3+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl3*(q13+q31)+mu*nl2*t4
		  fh[:,1,0] = -uinf2+nl2*(t3+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl3*(q23+q32)+mu*nl1*t4
		  fh[:,2,0] = -uh3+uinf3+x3

		  t2 = mu*nl2
		  t3 = mu*nl3
		  t4 = lambda*nl1
		  t5 = mu*nl1
		  t6 = mu*2.0
		  t7 = lambda+t6
		  t8 = lambda*nl2
		  fh_udg[:,0,0] = tau
		  fh_udg[:,0,3] = nl1*t7
		  fh_udg[:,0,4] = t2
		  fh_udg[:,0,5] = t3
		  fh_udg[:,0,6] = t2
		  fh_udg[:,0,7] = t4
		  fh_udg[:,0,9] = t3
		  fh_udg[:,0,11] = t4
		  fh_udg[:,1,1] = tau
		  fh_udg[:,1,3] = t8
		  fh_udg[:,1,4] = t5
		  fh_udg[:,1,6] = t5
		  fh_udg[:,1,7] = nl2*t7
		  fh_udg[:,1,8] = t3
		  fh_udg[:,1,10] = t3
		  fh_udg[:,1,11] = t8

		  fh_uh[:,0,0] = -tau
		  fh_uh[:,1,1] = -tau
		  fh_uh[:,2,2] = -1.0

	elif ib ==6:
		  fh[:,0,0] = -uh1+uinf1+x1
		  fh[:,1,0] = -uh2+uinf2+x2
		  fh[:,2,0] = -uinf3+nl3*(lambda*(q11+q22+q33+3.0)+mu*(q33+1.0)*2.0)+tau*(u3-uh3)+mu*nl1*(q13+q31)+mu*nl2*(q23+q32)

		  t2 = lambda*nl3
		  t3 = mu*nl1
		  t4 = mu*nl2
		  fh_udg[:,2,2] = tau
		  fh_udg[:,2,3] = t2
		  fh_udg[:,2,5] = t3
		  fh_udg[:,2,7] = t2
		  fh_udg[:,2,8] = t4
		  fh_udg[:,2,9] = t3
		  fh_udg[:,2,10] = t4
		  fh_udg[:,2,11] = nl3*(lambda+mu*2.0)

		  fh_uh[:,0,0] = -1.0
		  fh_uh[:,1,1] = -1.0
		  fh_uh[:,2,2] = -tau

	elif ib ==7:
		  fh[:,0,0] = -uinf1+nl1*(lambda*(q11+q22+q33+3.0)+mu*(q11+1.0)*2.0)+tau*(u1-uh1)+mu*nl2*(q12+q21)+mu*nl3*(q13+q31)
		  fh[:,1,0] = -uh2+uinf2+x2
		  fh[:,2,0] = -uh3+uinf3+x3

		  t2 = mu*nl2
		  t3 = mu*nl3
		  t4 = lambda*nl1
		  fh_udg[:,0,0] = tau
		  fh_udg[:,0,3] = nl1*(lambda+mu*2.0)
		  fh_udg[:,0,4] = t2
		  fh_udg[:,0,5] = t3
		  fh_udg[:,0,6] = t2
		  fh_udg[:,0,7] = t4
		  fh_udg[:,0,9] = t3
		  fh_udg[:,0,11] = t4

		  fh_uh[:,0,0] = -tau
		  fh_uh[:,1,1] = -1.0
		  fh_uh[:,2,2] = -1.0

	elif ib ==8:
		  fh[:,0,0] = -uh1+uinf1+x1
		  fh[:,1,0] = -uinf2+nl2*(lambda*(q11+q22+q33+3.0)+mu*(q22+1.0)*2.0)+tau*(u2-uh2)+mu*nl1*(q12+q21)+mu*nl3*(q23+q32)
		  fh[:,2,0] = -uh3+uinf3+x3

		  t2 = mu*nl1
		  t3 = mu*nl3
		  t4 = lambda*nl2
		  fh_udg[:,1,1] = tau
		  fh_udg[:,1,3] = t4
		  fh_udg[:,1,4] = t2
		  fh_udg[:,1,6] = t2
		  fh_udg[:,1,7] = nl2*(lambda+mu*2.0)
		  fh_udg[:,1,8] = t3
		  fh_udg[:,1,10] = t3
		  fh_udg[:,1,11] = t4

		  fh_uh[:,0,0] = -1.0
		  fh_uh[:,1,1] = -tau
		  fh_uh[:,2,2] = -1.0

	return fh, fh_udg, fh_uh
