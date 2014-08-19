#!/usr/bin/awk -f

# Get an analytical solution for gravity driven Poiseulle flow
# Usage:
# plot "<seq 0 1e-5 1e-3 | awk -v t=0.1 -v verbose=1 -f ./poisuille.awk"

function f(y, n,   aux) {
    return 4*F*L^2 / (nu*pi^3*(2*n+1)^3) * sin(pi*y/L * (2*n+1)) * exp(- (2*n+1)^2 * pi^2 * nu / L^2 * t)
}

function abs(x) {
    if (x>0) {
	return x
    }
    else {
	return -x
    }
}

BEGIN {
    if (!F) {F=1e-4}
    if (!nu) {nu=1e-6}
    if (!L) {L=1e-3}
    eps=1e-15
    if (verbose) {
	printf("F=%6.2e\n", F) > "/dev/stderr"
	printf("nu=%6.2e\n", nu) > "/dev/stderr"
	printf("L=%6.2e\n", L) > "/dev/stderr"
	printf("t=%6.2e\n", t) > "/dev/stderr"
	printf("eps=%6.2e\n", eps) > "/dev/stderr"
    }
    pi=3.14159265358979
}

{
    y = $1
    $1 = ""
    rest = $0
    if (length(t)==0) {
	fsum = -F/(2*nu) * y * (y - L)
    } else {
	n=0; fsum=0
	do {
	    fn = f(y, n)
	    fsum += fn
	    n++
	} while (abs(fn) > eps)
    }
    print y, rest, -F/(2*nu)*y*(y - L) - fsum
}