# build initial conditions
# eq[1]: 0 << x << 1;
# eq[2]: 0 << y << 1/2;
# eq[3]: y << x + 0.25;
# eq[4]: y << -x + 0.75;

function phase_dispatch(x, y) {
    if (x>0.5) {
	if (y > -x + 0.75) {
	    return 2
	} else {
	    return 3
	}
    } else {
	if (y > x - 0.25) {
	    return 1
	} else {
	    return 3
	}
    }
}


!in_atom && $0~/^Atoms / {
    in_atom = 1
    print
    getline
    print
    next
}

in_atom && NF==0 {
    in_atom = 0
    print
    next
}

!in_atom {
    print
    next
}

in_atom {
    x = $6
    y = $7
    $2 = phase_dispatch(x, y)
    print

    next
}
