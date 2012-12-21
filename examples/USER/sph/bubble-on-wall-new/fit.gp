last(n) = sprintf("%s/rg.dat", system(flast(n)))
flast(n) = sprintf("ls -1td data* | awk -v n=%i 'NR==n'", n)

order(n) = sprintf("%s/rg.dat", system(forder(n)))
forder(n) = sprintf("ls -1d data* | sort -g | awk -v n=%i 'NR==n'", n)

fg(n) = sprintf("ls -1d data* | sort -g | awk -v n=%i 'NR==n' | awk -vRS='-' -vFS='gy' '/gy/{print $2}' ", n)
g(n) = system(fg(n))+0.0

q1=g(1); q2=g(2); q3=g(3)
q4=g(4); q5=g(5); q6=g(6)
q7=g(7); q8=g(8)

set style data line
set macro
dx=1.0/40.0
y1='($7/sqrt($2))'
set term x11 1
plot order(1) u (sqrt($2)*sqrt(q1)):@y1,\
     order(2) u (sqrt($2)*sqrt(q2)):@y1,\
     order(3) u (sqrt($2)*sqrt(q3)):@y1,\
     order(4) u (sqrt($2)*sqrt(q4)):@y1,\
     order(5) u (sqrt($2)*sqrt(q5)):@y1,\
     order(6) u (sqrt($2)*sqrt(q6)):@y1,\
     order(7) u (sqrt($2)*sqrt(q7)):@y1,\
     order(8) u (sqrt($2)*sqrt(q8)):@y1

set term x11 2
plot order(1) u (sqrt($2)):@y1,\
     order(2) u (sqrt($2)):@y1,\
     order(3) u (sqrt($2)):@y1,\
     order(4) u (sqrt($2)):@y1,\
     order(5) u (sqrt($2)):@y1,\
     order(6) u (sqrt($2)):@y1,\
     order(7) u (sqrt($2)):@y1,\
     order(8) u (sqrt($2)):@y1
     
