c Switch whether we do pathdocmentation or not: iacpsw
c iacppt is the pointer in the arrays. 
c iacpcon is the pointers to starts of contours which are 
c xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
c iacpcon(i)+1 is the start index of each contour segment.
c iacpcon(i+1)-iacpcon(i) is the length of contour i.
c iacpcp is the current contour number, and when completed
c is 1+ the total number of contours. iacpcon(iacpcp) points to 
c the end of the last contour.
      integer iacpsw,iacppt,iacpcp,iacpmax,iacpconmax
      parameter (iacpmax=10000,iacpconmax=50)
      real xacp(iacpmax),yacp(iacpmax)
      integer iacpcon(iacpconmax)
      common /acpathcom/iacpsw,iacppt,iacpcp,xacp,yacp,iacpcon
