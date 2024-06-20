#!/bin/bash
# Generate figures for ionarbitrary.tex

export FIGDIRECTORY=inputfiles/

./chiofv ${FIGDIRECTORY}onegaus.dat -t -d
ps2png '-r300 plot0001.ps' maxwellian.png

./chiofv ${FIGDIRECTORY}onegaus.dat -d
ps2png '-r300 plot0001.ps' maxwellion.png

mv plot0002.ps maxwelok.ps
epstopdf maxwelok.ps
rm maxwelok.ps

./chiofv ${FIGDIRECTORY}flattail.dat -d
ps2png '-r300 plot0001.ps' flattailion.png
mv plot0002.ps flattailok.ps
epstopdf flattailok.ps
rm flattailok.ps

./chiofv ${FIGDIRECTORY}twogaussLowTe.dat -d
ps2png '-r300 plot0001.ps' twogaussLowTe.png

./chiofv ${FIGDIRECTORY}twonarrowgauss.dat -d
ps2png '-r300 plot0001.ps'  twonarrowgauss.png
mv plot0002.ps twonarrowok.ps
epstopdf twonarrowok.ps
rm twonarrowok.ps

./chiofv  -T20 -v1,-10,1,1 -Va-4 -d
ps2png '-r300 plot0001.ps' ionacoustic1.png 

./chiofv  -T10 -v1,-40,1,1 -Va-40 -d 
ps2png '-r300 plot0001.ps' ionacoustic2.png 
ps2png '-r300 plot0002.ps' ovk2.png 

./chiofv  -T100 -v1,-470,1,1 -Va-450 -i.8 -m20 -d
ps2png '-r300 plot0001.ps' ionacoustic3.png 

./chiofv  -T100 -v1,-900,1,1 -Va-880 -i1.5 -m20 -d
ps2png '-r300 plot0001.ps' ionacoustic4.png 
ps2png '-r300 plot0002.ps' ovk4.png 

./chiofv  -v1.,1.,.05,.05 -r1. -x.6 -Vi-.4 -Va1.4 -CP -d
ps2png '-r300 plot0001.ps' eamode.png 

./chiofv  -T99 -v.7,-1.3,1,1 -v.3,1.4,.5,.5 -i.9 -d
ps2png '-r300 plot0001.ps' langmuir.png 

./chiofv  -T99 -v.21,9,8,8 -v.19,3,8,8 -v.6,0,1,1 -x.6 -d
ps2png '-r300 plot0001.ps' eamode2.png 
