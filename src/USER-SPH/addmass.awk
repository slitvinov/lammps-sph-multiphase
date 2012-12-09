BEGIN {
    FS="//"
}

$1~/[^a-z]e\[/ {
    aux=$0
    gsub(/e\[/, "rmass[")
    print
    print aux
    next
}

{
    print
}

