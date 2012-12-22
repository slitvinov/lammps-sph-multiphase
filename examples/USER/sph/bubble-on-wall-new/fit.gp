last(n) = sprintf("%s/rg.dat", system(flast(n)))
flast(n) = sprintf("ls -1td data* | awk -v n=%i 'NR==n'", n)

order(n) = sprintf("%s/rg.dat", system(forder(n)))
forder(n) = sprintf("ls -1d data* | sort -g | awk -v n=%i 'NR==n'", n)

fg(n) = sprintf("ls -1d data* | sort -g | awk -v n=%i 'NR==n' | awk -vRS='-' -vFS='gy' '/gy/{print $2}' ", n)
g(n) = system(fg(n))+0.0

# get all gy as one string
gidx = system("ls  -d data*  | awk -vRS='-' -vFS='gy' '/gy/{print $2}' | xargs")
g(n) = word(gidx,n)


set style data line
set macro
dx=1.0/40.0
x1='(sqrt($2)*g(i)**0.50)'
#y1='($7/sqrt($2))'
y1='($7/sqrt($2))'
set term x11 1
plot [][:] for [i=1:2] order(i) u @x1:@y1 t g(i)

#set term x11 2
#plot [][:1] for [i=20:54:5] order(i) u (sqrt($2)):@y1 t g(i)
