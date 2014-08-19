# transform to c++

BEGIN {
    FS=","
}

{
    gsub("^\\[", "")
    gsub("\\]$", "")
    for (i=1; i<=NF; i++) {
	printf("b%i[%i] = %s\n", NR, i, $i)
    }
}