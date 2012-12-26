#!/usr/bin/awk -f 

# returns an analytical solution for temperature distribution 
# for infinit slab (see Cleary1999, eq. (35) 
# to test erfc function 
# in gnuplot
# plot "<awk -v test_erf=1 -v x=1.0 -f infslab.awk ", erfc(x)
# Usage
# seq 0.0 0.01 1.0 | awk -v k_l=1.0 -v rho_l=1e3 -v cv_l=1.0 -v Tl=0.0 -v Tr=1.0 -v xm=0.5 -v t=0.1 -f infslab.awk 
# seq 0.0 0.001 1.0 | awk -v k_l=1.0 -v rho_l=1e3 -v cv_l=1.0 -v Tl=0.0 -v Tr=1.0 -v xm=0.5 -v t=0.1 -f infslab.awk

function abs(x) {
    if (x>0) {
	return x
    } else {
	return -x
    }
}

# Handbook of Mathematical Functions: with Formulas, Graphs, and Mathematical Tables , 
# formula 7.1.26, error<1.5e-7
function erf(x,        sign, a1, a2, a3, a4, a5, t, y) {
    # save the sign of x
    sign = 1
    if (x < 0) {
        sign = -1
    }
    x = abs(x)
    if (x<5.0) {
	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911
	
	# A&S formula 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
    } else {
	y = 1.0
    }
    return sign*y
}

function erfc(x) {
    return 1.0 - erf(x)
}

BEGIN {
    if (test_erf) {
	dx=0.01
	for (i=-400; i<400; i++) {
	    print i*dx, erfc(i*dx)
	}
	exit
    }

    if (length(k_l)==0) {
	printf("k_l should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(rho_l)==0) {
	printf("rho_l should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(cv_l)==0) {
	printf("cv_l should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(k_r)==0) {
	# use "left" conductivity if not given
	k_r = k_l
    }
    if (length(rho_r)==0) {
	# use "left" density if not given
	rho_r = rho_l
    }
    if (length(cv_r)==0) {
	# use "left" heat capacity if not given
	cv_r = cv_l
    }
    if (length(Tl)==0) {
	printf("Tl should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(Tr)==0) {
	printf("Tr should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(xm)==0) {
	printf("xm should be given\n") > "/dev/stderr"
	exit(2)
    }
    if (length(t)==0) {
	printf("t should be given\n") > "/dev/stderr"
	exit(2)
    }
    alpha_l=k_l/(rho_l*cv_l)
    alpha_r=k_r/(rho_r*cv_r)
    printf("alpha_l = %e\n", alpha_l) > "/dev/stderr"
    printf("alpha_r = %e\n", alpha_r) > "/dev/stderr"
    Tc = (Tr-Tl) * sqrt(alpha_l) / ( sqrt(alpha_l) + sqrt(alpha_r) )
    printf("Tc = %e\n", Tc) > "/dev/stderr"
    printf("xm = %e\n", xm) > "/dev/stderr"
}

NF{
    x=$1
    if (x<xm) {
	relT= erfc( (xm - x)/ sqrt(alpha_l) / sqrt(t) / 2  )
    } else {
	relT= 1 + sqrt(alpha_l/alpha_r) * erf( (x - xm)/ sqrt(alpha_r)  / sqrt(t) / 2 )
    }
    print x, relT*Tc + Tl
}