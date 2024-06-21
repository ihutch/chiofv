!    Switch whether we do pathdocmentation or not: iacpsw
! iacppt is the pointer in the arrays. 
! iacpcon is the pointers to starts of contours which are 
! xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
! iacpcon(i)+1 is the start index of each contour segment.
! iacpcon(i+1)-iacpcon(i) is the length of contour i.
! iacpcp is the current contour number, and when completed
! is 1+ the total number of contours. iacpcon(iacpcp) points to 
! the end of the last contour.
      integer iacpsw,iacppt,iacpcp,iacpmax,iacpconmax
      parameter (iacpmax=10000,iacpconmax=50)
      real xacp(iacpmax),yacp(iacpmax)
      integer iacpcon(iacpconmax)
      common /acpathcom/iacpsw,iacppt,iacpcp,xacp,yacp,iacpcon
