reset 
# set term pdfcairo size 10,10 font "Times-New-Roman,40"
# set output
set terminal qt size 1500,1500
unset output

set key on box lc rgb "grey10" left top Left vertical reverse spacing 1.5 samplen 1 height -0.1 width -0.1 

set format x "%.2f"

set xrange[-1:10]
set yrange[-1:10]

steps=5000
xlim=10
ylim=20
set object rectangle from 0,0 to xlim,ylim fs empty border 3

while(1){
    do for[i=0:steps-1:10]{
        set title sprintf("Iteration %d",i)
        plot sprintf("../build/positions_step_%d.txt",i)  with points pt 7 ps 1 notitle
        pause 0.08
    }
}

exit
