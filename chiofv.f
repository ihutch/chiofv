c New version of chicomplex to normalize to T0 and allow different Te to
c be explored graphically. T0 is set to T_e on input and never
c changes. So the velocity scale labellings are relative to that T0 and
c in sound-speed units sqrt(T0/mi).  However, T_e can be considered to
c vary within the interpretation of the contours, keeping T0 fixed.

c Velocities are input and output in the same sqrt(T0/m_i) units.
c The range of velocities plotted is 2vw+vdmax-vdmin, where vw is 
c propto sqrt(max(gpar(3,1),gpar(4,1)). So first Gaussian sets plot width.
      program chiofv
      integer nv,npts
      real vw
c Number of velocities in integrals nv, number of Re(vp) points npts
c These are two different velocity arrays.
c Number of Im(vp) points: nim. 
      parameter (nv=400,npts=200,nimax=60,nim2=2*nimax+1,nmmax=-nimax)
c Input to the analysis routines on nv-mesh distrib, velocity, f^prime.
      real f(0:nv+2),v(0:nv+2),ww(0:nv+2)
c Intermediate analysis output on npts mesh
      real vp(npts),fv(npts),fvex(npts)
c ,fvex(npts)
c Final chi(vp) for complex grid.
      real cvreal(npts,-nimax:nimax),cvimag(npts,-nimax:nimax)
      real ceimag(npts,-nimax:nimax)
c Contouring work array:
      character cworka(npts,-nimax:nimax)
c The velocity step and start.
      real dv,v0,v1
      real Te
      real rmitome
c
c Things describing the input distribution as sum of Gaussians. 
      integer npar
      parameter (npar=4,ngmax=6)
c Parameters; density, velocity shift, vtperp, vtparallel.      
      real gpar(npar,ngmax)
      real fstat(ngmax)
      real theta
      complex cvint,vpec,chihat,chie
      external cvint,chihat
      integer nim,nmm,nvin
      logical laspect,lcolor,lextra,lgrowth,lcoptic,ltwotone,lthresh
      real vpimax,vpi(-nimax:nimax)
      character*50 cfilename

      data nim/nimax/nmm/nmmax/laspect/.false./lcolor/.true./
      data lextra/.false./lgrowth/.true./ltwotone/.true./
      data vpimax/.4/vpi/nim2*0/

c Default: Plot to screen
      call pfset(3)
      call dcharsize(.02)
      fflat=0.
c Silence warnings.
      cworka(1,1)=char(0)
c Extent of integration factor
      vwfac=6.
c Default gaussian data:
      ng=1
      gpar(1,1)=1.
      gpar(2,1)=0.
      gpar(3,1)=1.
      gpar(4,1)=1.
      theta=0.
      Te=1.
      rmitome=1836
      vdmin=0.
      vdmax=0.
      vw=vwfac*sqrt(max(gpar(3,1),gpar(4,1)))
      vrgfac=1.
c Default electron landau damping weight is what a maxwellian of
c temperature Te will give.
c Needs to account for the mass ratio 1836 by default.
      eldweight=sqrt(3.141593/(2.*rmitome))
c Smoothing of analytic continuation default setting.
      smfac=1.
c Noise amplitude:
      amp=.00
c Default no coptic data
      lcoptic=.false.
c---------------------------------------------------------------
c Parse command line for changes to defaults
      call parsecmdline(gpar,npar,ngmax,ng,vdmin,vdmax,vw,Te,theta
     $     ,vwfac,vpimax,nim,nmm,nimax,smfac,laspect,lcolor,eldweight
     $     ,vflat,wflat,fflat,lextra,amp,nw,lgrowth,lcoptic,cfilename
     $     ,ltwotone,vrgfac,lthresh,rmitome)
c---------------------------------------------------------------
      if(.not.lcoptic)then
c Construct the test distribution function as a sum of Gaussians.
         vw=vwfac*sqrt(max(gpar(3,1),gpar(4,1)))
         dv=(vdmax-vdmin+2.*vw)/nv
         v0=vdmin-vw
c Scale to sound-speed units assuming inputs relative to Ti=1:
         v0=v0/sqrt(Te)
         dv=dv/sqrt(Te)
         do i=1,ng
            gpar(2,i)=gpar(2,i)/sqrt(Te)
            gpar(3,i)=gpar(3,i)/(Te)
            gpar(4,i)=gpar(4,i)/(Te)
         enddo
c
         pne=0.
c Calculate total ion density
         do i=1,ng
            pne=pne+gpar(1,i)
         enddo
c Set the distribution function to be sum of maxwellians.
c      write(*,*)'dv,v0',dv,v0
         call fvinit(nv,dv,v0,f,ng,gpar,theta)
c Normalize it and initialize velocity.
         do i=0,nv+2
            f(i)=f(i)/pne
            v(i)=(v0+i*dv)
         enddo
         
         if(fflat.ne.0.)then
            iv=int((vflat-v0)/dv)
            iw=int(wflat/dv)
            call flatspot(nv,dv,f,iv,iw,fflat)
            write(*,'(a,3f6.2)')'Flatspot vflat,wflat,fflat'
     $           ,vflat,wflat,fflat
         endif
         call COPTICverif(nv,dv,v0,f)
c Distribution created.
      else
c External data
         open(13,file=cfilename, form='formatted',status='old',err=101)
         read(13,*) nvin
         read(13,*) dv
         read(13,*) v1
         v0=v1-dv
         if(nvin.gt.nv)then
           write(*,*) 'Input f longer than expected.','nvin=',nvin
           call exit(1)
         endif
c There's a big problem here. nv ought to be in the file.
c Also the velocity has to be a uniform grid. Really only v0 and dv
c ought to be specified.
         read(13,'(1g14.6)',end=102)(f(i),i=1,nvin)
c Also there must not be non-zero values at the ends. So we ought not
c to read the ends we should just put them to zero.
c Also it ought to be normalized. ...
         f(0)=0.
         f(1)=0.
c Set the end of the distribution function to zero.
         do i=nvin+1,nv+2
            f(i)=0.
         enddo
         do i=0,nv+2
            v(i)=(v0+i*dv)
         enddo
      endif
c Increment Smoothing of analytic continuation based on ranges.
      smfac=smfac+ 10.*(npts/(3.*nv))*dv/vpimax

c      write(*,*)'smfac=',smfac,' nv=',nv,' dv=',dv
c Correct the electron susceptibility for propagation angle theta
      if(cos(theta).gt.0.)then 
         eldweight=eldweight/cos(theta)
      else
         write(*,*)'CosTheta is non-positive',theta,cos(theta)
     $        ,' no angle factors.'
      endif

c Perhaps add noise?
      call addnoise(nv,f,amp)
c Perhaps smooth:
c      write(*,*)'nw',nw
      if(nw.gt.0)call trismooth(nv-3,nw,f(3),ww,f(3))
c-----------------------------------------------------------------
c Calculate the chi along just the real axis for thresholds
      if(lthresh)then
         call stationary(nv,f,ngmax,fstat,nstat)
         if(nstat.eq.3)then
c Find fractional depth assuming that the middle one is the local minimum.
            fmax=max(fstat(1),fstat(3))
            fmid=min(fstat(1),fstat(3))
            depth=(fmid-fstat(2))/fmid
            write(*,*)nstat,(fstat(i),i=1,nstat),depth
         else
            write(*,*)'Incorrect number of stationaries',nstat
         endif
         call chiion(v0,dv,f,nv,   fv,vp,npts,fmax
     $     ,nimax,0,0,vpi,vpimax,smfac,cvreal,cvimag)
!         write(*,'(i5,3f8.4)')(i,vp(i),cvreal(i,0),cvimag(i,0),i=1,npts)
         do i=npts/10,npts-1
            if(cvimag(i,0)*cvimag(i+1,0).le.0.and.cvreal(i,0).le.0)then
               ri=cvimag(i,0)/(cvimag(i,0)-cvimag(i+1,0))
               cvr0=ri*cvreal(i+1,0)+(1-ri)*cvreal(i,0)
               cvi0=ri*cvimag(i+1,0)+(1-ri)*cvimag(i,0)
               vp0=ri*vp(i+1)+(1-ri)*vp(i)
!               write(*,*)i,ri,vp0,cvr0
               write(*,*)depth,abs(gpar(2,1)-gpar(2,2))/2.,-Te/cvr0,
     $              gpar(1,2)/(gpar(1,1)+gpar(1,2)),gpar(3,2)/gpar(3,1)
               stop
            endif
         enddo
      endif
c-----------------------------------------------------------------
c Calculate chi-ion from it.
      call chiion(v0,dv,f,nv,   fv,vp,npts,fmax
     $     ,nimax,nim,nmm,vpi,vpimax,smfac,cvreal,cvimag)
c Get chi electron and add on real part, and imag for original plot.
      do j=1,npts
         do k=nmm,nim
            vefactor=1./sqrt(2.*rmitome)
            ct=cos(theta)
            if(ct.gt..01)then
               vefactor=vefactor/ct
            else
               write(*,*)'Dangerously low cos theta',ct,theta
               stop
            endif
            vpec=cmplx(vp(j),vpi(k))*vefactor
            chie=chihat(vpec)
            if(.not.ltwotone)then
               cvimag(j,k)=cvimag(j,k)+imag(chie)
               cvreal(j,k)=cvreal(j,k)+real(chie)
            endif
            ceimag(j,k)=imag(chie)
         enddo
      enddo
      if(.false.)then !Obsolete but saving for theta tests.
         call chicalc(v0,dv,f,nv,  fv,vp,npts,fmax,theta,rmitome
     $        ,nimax,nim,nmm,vpi,vpimax,eldweight,smfac,cvreal,cvimag)
      endif
c---------------------------------------------------------------

c Usually unused plots:
      if(lextra)then
         call realvpplot(vp,cvreal(1,0),cvimag(1,0),npts,nimax+1)
         call argdiagplot(cvreal(1,0),cvimag(1,0),npts,Te)
c         call vertslice(vpi(-nim),npts,nim, cvreal(1,-nim))
      endif

c Contour the susceptibility for stability analysis
      call multiframe(2,1,0)
      call ticnumset(7)
      vrgfac=.6
      koff=int(min(max(1.-vrgfac,0.)*npts/2.,npts-1.))
      nrmin=1+koff
      nrmax=npts-koff

      call plotfofv(vp(nrmin),fv(nrmin),nrmax-nrmin+1,fmax)
      if(.true.)then
c Plot the electron distribution function effective shape.
         vte2=2.*rmitome
         do i=1,npts
            fvex(i)=rmitome/sqrt(3.141593*vte2)
     $           *(exp(-vp(i)**2/vte2)-1.)
         enddo
         call dashset(2)
         call color(3)
         call polyline(vp(nrmin),fvex(nrmin),nrmax-nrmin+1)
         call legendline(.1,.1,0,' (f!de!d(v)-f!de!d(0))m!di!d/m!de!d')
         call color(15)
         call dashset(0)
      endif
      call chicontour(vp,vpi(nmm),npts,nrmin,nrmax,nim,nmm,cvreal(1,nmm)
     $     ,cvimag(1,nmm),ceimag(1,nmm),cworka,laspect,lcolor,ltwotone)
      call multiframe(0,0,0)
c Make the cvreal the total chi if it is not already (for acpplot).
      if(ltwotone)cvreal=cvreal+1.
c Calculate parameters along contours and plot:
      if(lgrowth)call acpplot(vp,vpi(nmm),npts,nim,nmm,cvreal(1,nmm)
     $     ,cvimag(1,nmm),cworka,laspect,lcolor)
      call exit(0)
 101  write(*,*)'Could not open coptic file',cfilename
      call exit(1)
 102  write(*,*)'Premature coptic file end'
      end
c End of main program.
c*****************************************************************
c*****************************************************************
c Calculate the real and imaginary parts of chi for a specified
c ion distribution function. Ion chi only, here.
      subroutine chiion(v0,dv,f,nv, fv,vp,npts,fmax
     $     ,nimax,nim,nmm,vpi,vpimax,smfac,cvreal,cvimag)

c On entry: v0, dv are the origin and step of the velocity for f(v),
c           f(v) is the distribution of length 0:nv+2.
c           nimax defines the imaginary velocity allocated length
c           nim,nmm define its actual array length.
c           vpimax defines the imag-velocity range,
c           smfac the smoothing factor.
c
c On exit:
c          Arrays of real velocity of length npts are in vp, 
c          fv is the distribution, 
c          vpi contains the imaginary velocity array. 
c          cvreal, cvimag contain the real and imaginary part of
c          chi.(k\lambda_{De})^2.
c
c Input to the analysis routines on nv-mesh distrib, velocity, f^prime.
      real f(0:nv+2)
c      real v(0:nv+2) is not used, uniform grid v0,dv instead
c Intermediate analysis output on npts mesh
      real vp(npts),fv(npts)
      real cvreal(npts,-nimax:nimax),cvimag(npts,-nimax:nimax)

c Fourier transform variables
      parameter (nptsmax=200)
      parameter (efold2=10.)
      integer lensav,nfftwork
      complex ftfp(nptsmax),ftbk(nptsmax)
c Has to be at least 2*N + INT(LOG(REAL(N))/LOG(2.)) + 4. , (N=npts)
      parameter (lensav=3*nptsmax,nfftwork=2*nptsmax)
      real wsave(lensav),fftwork(nfftwork)
      real vpimax,vpi(-nimax:nimax)
      complex cvint,cvi,vpc
      complex chihat
c      logical lchaprox,lchieinc
      external cvint,chihat

c      write(*,*)nv,npts,nim,nmm,nimax,v0

      if(npts.gt.nptsmax)stop 'npts .gt. nptsmax'
c-----------------------------------------------------------------
      fmax=0.
      do j=1,npts
         vp(j)=v0+dv+dv*(nv-1)*(j-1.)/(npts-1.)
         vi=1.+(nv-1)*(j-1.)/(npts-1.)
         i=int(vi)
         xi=vi-i
         if(i.gt.nv)write(*,*)'Too large i',i,j
c Interpolated values of f(v) and f'(v)
         fv(j)=(f(i)*(1.-xi)+f(i+1)*xi)
         if(fv(j).gt.fmax)fmax=fv(j)
c Derivative indices
         ixd=nint(xi)
         id=i+ixd
         xid=xi+0.5-ixd
c--------------------
c Real values of vp.
c For a range of vp's, calculate the pvintegral and imaginary part.
         cvimag(j,0)=-3.1415926*(
     $        (f(id)-f(id-1))/(dv)*(1.-xid)
     $        +(f(id+1)-f(id))/(dv)*xid)
c The ion principal value integral
         cvreal(j,0)=-pvint(v0,dv,f,nv,vp(j))
c--------------------
c Positive values of imaginary part of vp.
c Integration along the real axis gives chi integral.
         do k=1,nim
            vpi(k)=vpimax*k/nim
c Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
c Integral (which is chi.(k\lambda_{De})^2
            cvi=-cvint(v0,dv,f,nv,vpc)
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
c--------------------
c Also get chi for negative imaginary part.
         do k=-1,nmm,-1
            vpi(k)=vpimax*k/abs(nmm)
c Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
c Integral (which is chi.(k\lambda_{De})^2 plus the electron terms.
            cvi=-cvint(v0,dv,f,nv,vpc)
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
      enddo
c---------------------------------------------------------------
c Analytically continue fp into the lower half plane to give pole 
c contribution proportional to df/dv.
c Copy f(v) data. 
c      write(*,*)'Analytic Continuation'
      do j=1,npts
         ftfp(j)=fv(j)
      enddo
c      write(*,*)'Calling cfft1i',npts,lensav
      call cfft1i(npts,wsave,lensav,ier)
      if(ier.ne.0)stop 'cfft1i error'
c      write(*,*)'Calling cfft1f',nv+3,lensav,ier
      call cfft1f(npts,1,ftfp,npts,wsave,lensav,fftwork,nfftwork,ier)
      if(ier.ne.0)stop 'cfft1f error'      
c Now ftfp contains the fourier transform of f(). v<->t p<->omega
      dt=(nv*dv)/npts    ! Spacing of the (npts) v-values: delta v in tex.
c 1/dt is the maximum p 
      dodt=2.*3.1415926/npts   ! delta omega * delta v ->  dp*dt
      do i=0,nmm,-1        ! For each imaginary vp 
c Bandwidth limit
         do j=1,npts       ! Along the real direction
            kmod=mod(j-1+npts/2,npts)-npts/2    ! Reordered real index
            p=kmod*dodt    ! p*dt
c  Differentiate (once only)
            if(i.eq.0)ftfp(j)=ftfp(j)*cmplx(0.,1./dt)*sin(p)
            pv=p*vpi(i)/dt
c smfac is approximately 1. The max pzi should be about 5 (e-foldings).
c The max of the following expression is  efold2/(2*smfac) which is
c ~5 when efold2=10.
            pzi=pv/(1.+(pv*smfac/efold2)**2)
            ftbk(j)=ftfp(j)     *exp(pzi)
         enddo
c For testing only:
c            ftbk(j)=ftfp(j)
         call cfft1b(npts,1,ftbk,npts,wsave,lensav,fftwork,nfftwork,ier)
c then ftbk contains the analytic continuation of f' to the
c corresponding value of vpi, (as a fn of vpr) which we use as the pole
c contribution to add (times 2\pi i) to the negative vpi cases.
c Check of (last) transform. 
c       write(*,'(5f10.4)')(fv(j),ftbk(j),ftfp(j),j=1,npts)
         do j=1,npts
            if(i.eq.0)then 
            else
c Add appropriately weighted to cv
               cvimag(j,i)=cvimag(j,i)-2.*3.1415926*real(ftbk(j))
               cvreal(j,i)=cvreal(j,i)-2.*3.1415926*imag(ftbk(j))
            endif
         enddo
      enddo
c      write(*,*)'Analytic Continuation Completed'
      end
c**********************************************************************
c*****************************************************************
c Calculate the real and imaginary parts of chi for a specified
c ion distribution function and electron temperature.
      subroutine chicalc(v0,dv,f,nv,  fv,vp,npts,fmax,theta,rmitome
     $     ,nimax,nim,nmm,vpi,vpimax,eldweight,smfac,cvreal,cvimag)

c On entry: v0, dv are the origin and step of the velocity for f(v),
c           f(v) is the distribution of length 0:nv+2.
c           nimax defines the imaginary velocity allocated length
c           nim,nmm define its actual array length.
c           vpimax defines the imag-velocity range,
c           eldweight the electron landau damping weight
c           theta is the angle of propagation to the B-field
c           rmitome is the ratio of ion to electron mass
c           smfac the smoothing factor.
c
c On exit: fp contains derivative of f.  
c          Arrays of real velocity of length npts are in vp, 
c          fv is the distribution, 
c          vpi contains the imaginary velocity array. 
c          cvreal, cvimag contain the real and imaginary part of
c          chi.(k\lambda_{De})^2.
c
c Input to the analysis routines on nv-mesh distrib, velocity, f^prime.
      real f(0:nv+2)
c      real v(0:nv+2) is not used, uniform grid v0,dv instead
c Intermediate analysis output on npts mesh
      real vp(npts),fv(npts)
      real cvreal(npts,-nimax:nimax),cvimag(npts,-nimax:nimax)

c Fourier transform variables
      parameter (nptsmax=200)
      parameter (efold2=10.)
      integer lensav,nfftwork
      complex ftfp(nptsmax),ftbk(nptsmax)
c Has to be at least 2*N + INT(LOG(REAL(N))/LOG(2.)) + 4. , (N=npts)
      parameter (lensav=3*nptsmax,nfftwork=2*nptsmax)
      real wsave(lensav),fftwork(nfftwork)
      real vpimax,vpi(-nimax:nimax)
      complex cvint,cvi,vpc,vpec
      complex chihat,chie
      logical lchaprox
      external cvint,chihat

c      write(*,*)nv,npts,nim,nmm,nimax,v0

c      lchaprox=.true.
      lchaprox=.false.
      if(npts.gt.nptsmax)stop 'npts .gt. nptsmax'
c-----------------------------------------------------------------
      vefactor=1./sqrt(2.*rmitome)
      ct=cos(theta)
      if(ct.gt..01)then
         vefactor=vefactor/ct
         write(*,'(a,2f8.5)')'vefactor=',vefactor,ct
      else
         write(*,*)'Dangerously low cos theta',ct,theta
      endif
      fmax=0.
      do j=1,npts
         vp(j)=v0+dv+dv*(nv-1)*(j-1.)/(npts-1.)
         vi=1.+(nv-1)*(j-1.)/(npts-1.)
         i=int(vi)
         xi=vi-i
         if(i.gt.nv)write(*,*)'Too large i',i,j
c Interpolated values of f(v) and f'(v)
         fv(j)=(f(i)*(1.-xi)+f(i+1)*xi)
         if(fv(j).gt.fmax)fmax=fv(j)
c Derivative indices
         ixd=nint(xi)
         id=i+ixd
         xid=xi+0.5-ixd
c--------------------
c Real values of vp.
c For a range of vp's, calculate the pvintegral and imaginary part.
         if(lchaprox)then
c f' also includes the electron contribution. 
            cvimag(j,0)=-3.1415926*(
     $           (f(id)-f(id-1))/(dv)*(1.-xid)
     $           +(f(id+1)-f(id))/(dv)*xid)
     $           +vp(j)*eldweight
c The ion principal value integral plus the electron contribution
            cvreal(j,0)=-pvint(v0,dv,f,nv,vp(j))+1.
         else
            vpec=complex(vp(j),0.)*vefactor
            chie=chihat(vpec)
c            write(*,'(a,4f12.6)')'chie',chie,1.,vp(j)*eldweight
            cvimag(j,0)=-3.1415926*(
     $           (f(id)-f(id-1))/(dv)*(1.-xid)
     $           +(f(id+1)-f(id))/(dv)*xid)
     $           + imag(chie)
            cvreal(j,0)=-pvint(v0,dv,f,nv,vp(j))+real(chie)
         endif
c--------------------
c Positive values of imaginary part of vp.
c Integration along the real axis gives chi integral.
         do k=1,nim
            vpi(k)=vpimax*k/nim
c Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
c Integral (which is chi.(k\lambda_{De})^2 plus the electron terms.
            cvi=-cvint(v0,dv,f,nv,vpc)
            if(lchaprox)then
c Plus the approximated electron susceptibility:
               cvi=cvi+1.+cmplx(0.,vp(j))*eldweight
            else
               vpec=vpc*vefactor
               chie=chihat(vpec)
               cvi=cvi+chie
            endif
c (An alternative would be to use the full dispersion function for 
c oblique propagation. The above does not add on real part of i vpc.)
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
c--------------------
c Also get chi for negative imaginary part.
         do k=-1,nmm,-1
            vpi(k)=vpimax*k/abs(nmm)
c Complex v_p.
            vpc=cmplx(vp(j),vpi(k))
c Integral (which is chi.(k\lambda_{De})^2 plus the electron terms.
            cvi=-cvint(v0,dv,f,nv,vpc)
            if(lchaprox)then
c Plus the approximated electron susceptibility:
               cvi=cvi+1.+cmplx(0.,vp(j))*eldweight
            else
               vpec=vpc*vefactor
               chie=chihat(vpec)
               cvi=cvi+chie
            endif
            cvreal(j,k)=real(cvi)
            cvimag(j,k)=imag(cvi)
         enddo
      enddo
c      write(*,*)nim,nmm,vpimax
c      write(*,*)vpi
c---------------------------------------------------------------
c Analytically continue fp into the lower half plane to give chi.
c Copy f(v) data. 
c      write(*,*)'Analytic Continuation'
      do j=1,npts
         ftfp(j)=fv(j)
      enddo
c      write(*,*)'Calling cfft1i',npts,lensav
      call cfft1i(npts,wsave,lensav,ier)
      if(ier.ne.0)stop 'cfft1i error'
c      write(*,*)'Calling cfft1f',nv+3,lensav,ier
      call cfft1f(npts,1,ftfp,npts,wsave,lensav,fftwork,nfftwork,ier)
      if(ier.ne.0)stop 'cfft1f error'      
c Now ftfp contains the fourier transform of f().
      dt=(nv*dv)/npts
      dodt=2.*3.1415926/npts
      do i=0,nmm,-1        ! For each imaginary vp 
c Bandwidth limit
         do j=1,npts       ! Along the real direction
            kmod=mod(j-1+npts/2,npts)-npts/2    ! Reordered real index
            p=kmod*dodt    ! p*dt
c  Differentiate (once only)
            if(i.eq.0)ftfp(j)=ftfp(j)*cmplx(0.,1./dt)*sin(p)
            pv=p*vpi(i)/dt
c smfac is approximately 1. The max pzi should be about 5 (e-foldings).
c The max of the following expression is  efold2/(2*smfac) which is
c ~5 when efold2=10.
            pzi=pv/(1.+(pv*smfac/efold2)**2)
            ftbk(j)=ftfp(j)     *exp(pzi)
         enddo
c For testing only:
c            ftbk(j)=ftfp(j)
         call cfft1b(npts,1,ftbk,npts,wsave,lensav,fftwork,nfftwork,ier)
c then ftbk contains the analytic continuation of f' to the
c corresponding value of vpi, (as a fn of vpr) which we use as the pole
c contribution to add (times 2\pi i) to the negative vpi cases.
c Check of (last) transform. 
c       write(*,'(5f10.4)')(fv(j),ftbk(j),ftfp(j),j=1,npts)
         do j=1,npts
            if(i.eq.0)then 
            else
c Add appropriately weighted to cv
               cvimag(j,i)=cvimag(j,i)-2.*3.1415926*real(ftbk(j))
               cvreal(j,i)=cvreal(j,i)-2.*3.1415926*imag(ftbk(j))
            endif
         enddo
      enddo

c      write(*,*)'Analytic Continuation Completed'
      end
c**********************************************************************
c*********************************************************************
c      real function pvint(v0,dv,f,nv,vp)
c Integrate a distribution function f to obtain the susceptibility
c principal value integral defined as I(v_p)= \int {df(v)\over
c dv}/(v-v_p) dv. Here v_p is the phase velocity omega/k. f is the 1-d
c distribution function along v, and the resultant must be multiplied by
c (\omega_{pi}/k)^2/n_i to yield susceptibility.
c 
c The principal value integral is implemented by integrating by parts
c twice to obtain the form I = \int (v-v_p)(ln|v-v_p|-1)
c f'''(v)dv. (Which requires that f', f'' are zero at the extrema of the
c integration.  Then the value of f''' is evaluated by natural finite
c differences on a uniform mesh:
c
c  f'''_{i+1/2} = ([f_{i+2}+f_i-2f_{i+1}]- [f_{i+1}+f_{i-1}-2f_i])/dv^3
c 
c This is then integrated with I = \sum dv x value at i+1/2.  The
c limitations arise from the high order difference and rounding.  But
c the velocity function is non-singular once integrated, so no big
c difficulty arises in evaluating the integral through the former pole.
c f(v) ought to be zero for at least two values at the ends of its
c range.  The integral is done only over the non-zero values. The
c velocity grid must be uniform, and only v_0,dv is passed. 2:nv is the
c non-zero range. Convergence is second order in dv, i.e. in 1/nv.
c
c If velocities are in units of sqrt(2T_i/m_i), and the total i-th species
c density is given by n_i = \int f dv, then the i-susceptibility is
c \chi_i(\omega,k)= [Z_i^2 T_e/T_i] [-I(v_p)/2]/(n_e k^2\lambda_{De}^2)
c
c The imaginary part of the susceptibility in these units comes from the
c above equation with I(v_p) replaced by \pi df/dv giving half the residue.
c 
c Alternatively, if velocities are in units of sqrt(T_e/m_i), then
c \chi_i(\omega,k)= [Z_i^2] [-I(v_p)]/(n_e k^2\lambda_{De}^2)
c
c In either case, if f_i is normalized, one must multiply by n_i.
c
c Return pvint(=I) for the distribution function f, which is an array of
c values on a velocity grid starting at v0 (at index 0), spaced by dv,
c and extending to index nv+2. The first two and last two values would
c normally be zero, so as to satisfy boundary conditions.
c 
      real function pvint(v0,dv,f,nv,vp)
      integer nv
      real v0,vp,dv
      real f(0:nv+2)
      parameter (small=1.e-20)
      pvint=0.
      do i=1,nv
         v=v0+(i+0.5)*dv
         vd=v-vp
         if(abs(vd).lt.small)then 
            vfn=0.
         else
            vfn=vd*(alog(abs(vd))-1)
         endif
         pvint=pvint+vfn*((f(i+2)-f(i-1))-3.*(f(i+1)-f(i)))
      enddo
      pvint=pvint/dv**2
      end
c*********************************************************************
c      complex function cvint(v0,dv,f,nv,vp)
c Perform the real and imaginary susceptibility integrations for a complex
c \omega/k =v_p with non-zero imaginary part, along the real v-axis.
c Denoting the real and imaginary parts of vp by vr, and vi respectively,
c the integral is
c J = \int (df/dv) (v-vr)/[(v-vr)^2+vi^2] df
c         + i \int (df/dv) vi /[(v-vr)^2 + vi^2] df
c   = -\int (d^2f/dv^2)[ln[{(v-vr)/vi}^2+1]/2 + i arctan{(v-vr)/vi}] dv.
c
c Interpolation to f''(i) = (f(i-1)-2*f(i)+f(i+1))/dv^2 
c So f''(i+1/2)=(f(i-1)-f(i)-f(i+1)+f(i+2))/2.dv^2
c
c f must go to zero at the ends of the integration to make this form valid.

      complex function cvint(v0,dv,f,nv,vp)
      complex vp
      integer nv
      real v0,dv
      real f(0:nv+2)
      parameter (small=1.e-6)
      real cr,ci

      vr=real(vp)
      vi=imag(vp)
      if(abs(vi).lt.small)then
         write(*,*)'cvint error: vi too small',vi
         cvint=cmplx(0.,0.)
         return
      endif
      cr=0.
      ci=0.
      do i=1,nv
c         v=v0+(i+0.5)*dv
c Four-point formula f''= this combination /2dv^2
c         fpp=(f(i-1)-f(i)-f(i+1)+f(i+2))

         v=v0+(i)*dv
cc Three-point formula. f''= this combination/dv^2
         fpp=(f(i-1)-2.*f(i)+f(i+1))
         va=(v-vr)/vi
         cr=cr-alog(va**2+1.)*fpp
         ci=ci-atan(va)*fpp
      enddo
c Four-point:
c      cr=cr*0.25/dv
c      ci=ci*0.5/dv
c Three-point:
      cr=cr*0.5/dv
      ci=ci*1./dv
      
      cvint=cmplx(cr,ci)
      end

c********************************************************************
c Return the 1-D distribution function in array f(0:nv+2) for nv+3
c velocity points from v0 to v0+(nv+2)*dv being the sum of ng gaussians
c whose parameters are in gpar(*,ig)=(den, vel, vtperp, vtparallel) 
c At an angle of velocity component to drift velocity, theta.

c gpar(3,4) here have already been scaled dividing by Te and velocities
c dv, v0, gpar(2) by sqrt(T0).  Consequently the maxwellian is
c exp(-v^2/(2*vt2)), vt2=Ti/T0 and the velocities v0 dv etc can be
c considered to be in units of sqrt(T0/mi).
      subroutine fvinit(nv,dv,v0,f,ng,gpar,theta)
      integer nv
      real dv,v0
      real f(0:nv+2)
      integer npar
      parameter (npar=4)
c Parameters; den, vel, vtperp, vtparallel.      
      real gpar(npar,ng)
      real theta

      do i=0,nv+2
         f(i)=0.
      enddo
      ct=cos(theta)
      st=sin(theta)
      do k=1,ng
         vt2=gpar(3,k)*st**2+gpar(4,k)*ct**2
         fnorm=1/(sqrt(2.*3.141593*vt2))
         do i=0,nv+2
            v=v0+i*dv
            f(i)=f(i)+gpar(1,k)*fnorm*exp(-(v-gpar(2,k)*ct)**2
     $           /(2.*vt2))
         enddo
      enddo
      end
c*********************************************************************
      subroutine realvpplot(vp,cvreal,cvimag,npts,nimax)
      integer npts,nimax
      real vp(npts),cvreal(npts,0:nimax),cvimag(npts,0:nimax)
      call autoplot(vp,cvreal(1,0),npts)
      call polymark(vp,cvreal(1,0),npts,1)
      call vecw(vp(1),0.,0)
      call vecw(vp(npts),0.,1)
c      call axlabels('v!dp!d/(2Ti/m_i)!u1/2!u',
c     $     '!Ax!@ .T!di!d(k!Al!@!dD!d)!u2!u/T!de!d')
      call axlabels('v!dp!d/(T!de!d/m!di!d)!u1/2!u',
     $     '!Ax!@ .(k!Al!@!dD!d)!u2!u')
      call legendline(.7,.9,0,'!AR!@(!Ax!@!di!d)')
      call boxtitle
     $     ('Susceptibility along the real !Aw!@/k=v!dp!d axis.')
      call winset(.true.)
      call dashset(2)
      call color(2)
      call polyline(vp,cvimag(1,0),npts)
      call legendline(.7,.85,0,'!AI!@(!Ax!@!di!d)')
c cv lines:
      call color(6)
      call dashset(3)
      call polyline(vp,cvreal(1,1),npts)
      call dashset(4)
      call polyline(vp,cvimag(1,1),npts)
      call color(15)
      call dashset(0)
      call pltend()
      end
c******************************************************************
      subroutine argdiagplot(pv,fvp,npts,Te)
      integer npts
      real pv(npts),fvp(npts),Te
c Argand diagram plot.
      call autoplot(pv,fvp,npts)
      call axlabels('Real','Imaginary')
      call winset(.true.)
c Set the space between labels to .05 normal units and line drawing off.
      call labeline(.05,0.,0,'',-99)
      call labeline(pv,fvp,npts,'>',1) 
      call labeline(.005,0.,0,'',-99)
      call labeline(pv,fvp,npts,'"',1)
      call labeline(0.,0.,0,'',-99)
c Draw axes through origin
      call vecw(-Te,0.,0)
      call vecw(Te,0.,1)
      call vecw(0.,-Te,0)
      call vecw(0.,Te,1)
      call pltend()
      end
c*******************************************************************
c Contour the real and imaginary parts of chi on 2D grid of vp
c vp(1:npts) vpi(nmm:nim) cworka is a work array at least as big.
c laspect if true says plot conserving aspect ratio.  
      subroutine chicontour(vp,vpi,npts,nrmin,nrmax,nim,nmm,cvreal
     $     ,cvimag,ceimag,cworka,laspect,lcolor,ltwotone)
      integer npts,nim,nmm,nrmin,nrmax
      real vp(npts),vpi(nmm:nim),cvreal(npts,nmm:nim),cvimag(npts,
     $     nmm:nim),ceimag(npts,nmm:nim)
      character cworka(npts,nmm:nim)
      logical laspect,lcolor
      logical ltwotone

      include 'acpathcom.f'

      parameter (npixel=240)
      integer ired(npixel),igreen(npixel),iblue(npixel)
      parameter (ncl=9)
      real zcl(-ncl:ncl)
c Cunning contours so as to be able to do coloring
      data zcl/-120.,-50.,-20.,-10.,-5.,-2.,-1.,-.5
     $     ,-.2,0.,.2,.5,1.,2.,5.,10.,20.,50.,120./
      if(laspect)then
         call pltinaspect(vp(nrmin),vp(nrmax),vpi(nmm),vpi(nim))
      else
         call pltinit(vp(nrmin),vp(nrmax),vpi(nmm),vpi(nim))
      endif

      if(lcolor)then
         do k=1,npixel 
            if(k.le.npixel/2-1.or.(k.le.npixel/2.and..not.ltwotone))then
               ired(k)=65535
               iblue(k)=65535
               igreen(k)=65535
            elseif(k.le.npixel/2)then
               ired(k)=58000
               iblue(k)=58000
               igreen(k)=58000
            else
               ired(k)=50000
               iblue(k)=50000
               igreen(k)=50000
            endif
         enddo
         call accisgradset(ired,igreen,iblue,npixel)
         icl=2*ncl+1
c Vector, Color, Omit lines.
         icsw=1+16+32+64
         call contourl(cvreal(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim
     $        -nmm+1,zcl(-ncl),icl,vp(nrmin),vpi,icsw)
      endif

c Draw zero line
      call color(12)
      call vecw(vp(nrmin),0.,0)
      call vecw(vp(nrmax),0.,1)
      call color(3)
c Vector, negative real part contours.
      icsw=1
      icl=ncl+1
      call contourl(cvreal(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim-nmm
     $     +1,zcl(-ncl),icl,vp(nrmin),vpi,icsw)
      call color(15)
      call axis()
      call axis2()
c Correct the cvreal
!      if(ltwotone)cvreal=cvreal+1.
c Zero imaginary part contour,
      icsw=1
      icl=1
      if(.not.ltwotone)iacpsw=1    ! Path documentation
      call contourl(cvimag(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim-nmm
     $     +1,zcl(0),icl,vp(nrmin),vpi,icsw)
      if(ltwotone)then
         call axlabels('!AR!@(v!dp!d) /(T!d0!d/m!di!d)!u1/2!u',
     $        '  !AI!@(v!dp!d) /(T!d0!d/m!di!d)!u1/2!u')
         call legendline(.6,.9,258,
     $        '(k!Al!@!dD!d)!u2!u!Ax!@!di!d contours')
         do k=nmm,nim
            do j=1,npts
               cvimag(j,k)=cvimag(j,k)+ceimag(j,k)
            enddo
         enddo
c Second imaginary contour
         call color(1)
         icl=-1
c Turn on path documentation.
         iacpsw=1
         iacpcp=1
         iacpcon(iacpcp)=0
         iacppt=0
         call contourl(cvimag(nrmin,nmm),cworka,npts,nrmax-nrmin+1,nim
     $        -nmm+1,zcl(0),icl,vp(nrmin),vpi,icsw)
         call legendline(.77,.65,0,'imag(!Ax!@)')
         iacpsw=0
         call color(15)
      else
         call axlabels('!AR!@(v!dp!d) /(T!de!d/m!di!d)!u1/2!u',
     $        '  !AI!@(v!dp!d) /(T!de!d/m!di!d)!u1/2!u')
         call legendline(.6,.9,258,'(k!al!@!dDe!d)!u2!u!Ax!@ contours')
      endif
      call legendline(.77,.75,0,'imag')
      call color(3)
      call legendline(.77,.85,0,'real')
      call color(15)
      call pltend()
      end
c*******************************************************************
c*******************************************************************
      subroutine plotfofv(vp,fv,npts,fmax)
      integer npts
      real vp(npts),fv(npts),fmax
      yt=1.1*fmax
      yb=-0.18*yt
      call pltinit(vp(1),vp(npts),yb,yt)
      call axptset(0.,abs(yb)/(yt-yb))
      call axis()
      call axis2()
      call polyline(vp,fv,npts)
      call axlabels('','f(v)')
      end
c********************************************************************
c Calculate the result of the path/contour following.
      subroutine acpplot(vp,vpi,npts,nim,nmm,cvreal,cvimag,cworka
     $     ,laspect,lcolor)
      integer npts,nim,nmm
      real vp(npts),vpi(nmm:nim),cvreal(npts,nmm:nim),cvimag(npts,
     $     nmm:nim)
      character cworka(npts,nmm:nim)
      logical laspect,lcolor

c iacpcon is the pointers to starts of contours which are 
c xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
c iacpcon(i)+1 is the start index of each contour segment.
c iacpcon(i+1)-iacpcon(i) is the length of contour i.
c iacpcp is the current contour number, and when completed
c is 1+ the total number of contours. iacpcon(iacpcp) points to 
c the end of the last contour.
      include 'acpathcom.f'
      real curk(iacpmax),curfr(iacpmax),curfi(iacpmax)
      character*20 label

      wimfac=.1
      ymin=-1.5
      ymax=2.
      xmax=3.2
      xmin=0.
      call pltinit(xmin,xmax,ymin,ymax)
      call axis()
      call axlabels('k!Al!@!dD!d','!AR!@(!Aw!@)/!Aw!@!dpi!d')
      call axptset(1.,0.)
      call ticrev()
      call altyaxis(wimfac,wimfac)
      call axlabels('','!AI!@(!Aw!@)/!Aw!@!dpi!d')
      call ticrev()
      call axptset(0.,0.)
      call winset(.true.)
      call vecw(0.,0.,0)
      call vecw(5.,0.,1)
      call legendline(.02,.95,0,'!AR!@(!aw!@)')
      call dashset(2)
      call legendline(.7,.95,0,'!AI!@(!Aw!@)')
      call dashset(0)
c Do over contour (segments)
      curfimx=-1.e6
      ncur=0
      do i=1,iacpcp-1
         icur=0
         id=0
         vpimx=-1.e6
         do j=1,iacpcon(i+1)-iacpcon(i)
c Convert x,y values to fractional indexes then interpolate for chi
            xf=(1.+(xacp(iacpcon(i)+j)-vp(1))/(vp(npts)-vp(1))*(npts-1))
            yf=(1.+(yacp(iacpcon(i)+j)-vpi(nmm))/(vpi(nim)-vpi(nmm))
     $           *(nim-nmm))
            chirept=amultilin2(cvreal(1,nmm),npts,npts,nim-nmm+1,
     $           xf,yf,ierr)
c Add only if in the negative chi-real region
            if(chirept.lt.0.)then
c Add point to curve
               icur=icur+1
               curk(icur)=sqrt(-chirept)
c Real part of this phase velocity
c               vprm=(vp(1) +(xacp(iacpcon(i)+j)-1)*(vp(npts)-vp(1))
c     $              /(npts-1))
               vprm=xacp(iacpcon(i)+j)
c Hence real part of frequency k.v_p:
               curfr(icur)=curk(icur)*vprm
c Imaginary part of phase velocity
               vpim=yacp(iacpcon(i)+j)
c Imaginary part of frequency divided by scaling factor.
               curfi(icur)=curk(icur)*vpim/wimfac
               if(vpim.gt.vpimx)vpimx=vpim
c This is scaled wrong because of wimfac.
               if(curfi(icur).gt.curfimx)curfimx=curfi(icur)*wimfac
c Decide the position of the line's label:
               xpos=xmax*(.6-.05*i)
               if(curk(icur).lt.xpos .and. id.eq.0)id=icur
c               write(*,'(2i4,7f9.4)')i,j,xacp(iacpcon(i)+j)
c     $              ,yacp(iacpcon(i)+j),curk(icur),curfr(icur)
c     $              ,curfi(icur),vprm,vpim
            endif
         enddo
c Segment completed
         if(id.gt.0.and.icur.gt.1 .and. vpimx.gt.-0.05)then
            ncur=ncur+1
c There is a curve to plot. Plot it.
            call dashset(0)
            call color(ncur)
            call iwrite(ncur,iwdth,label)
c            call labeline(curk,curfr,icur,label,10)
c            id=0.7*icur
            call polyline(curk,curfr,icur)
            call jdrwstr(wx2nx(curk(id)),wy2ny(curfr(id)),label,0.)
            call dashset(2)
            call polyline(curk,curfi,icur)
            call jdrwstr(wx2nx(curk(id)),wy2ny(curfi(id)),label,0.)
c            write(*,*)'   k       Re(o)       Im(o)  Seg',i
c     $           ,j
c            do k=1,icur
c               write(*,'(3f10.4)')curk(k),curfr(k),curfi(k)
c            enddo
         endif
      enddo
      call pltend()
      end
c******************************************************************
c Use interp at two columns to interpolate values; then interpolate
c those two values. Hence giving a multilinear interpolation in 2D.
      real function amultilin2(data,ld,nx,ny,xf,yf,ierr)
      real data(ld,ny)
      real xf,yf
      ierr=0
      ix=int(xf)
      iy=int(yf)
      yff=yf-iy
      xff=xf-ix
c Allow values to the upper limit
      if(ix.eq.nx)then
         ix=ix-1
         xff=1.
      endif
      if(iy.eq.ny)then
         iy=iy-1
         yff=1.
      endif
c But not beyond.
      if(ix.lt.1.or.ix.ge.nx.or.iy.lt.1.or.iy.ge.ny)then
         ierr=1
         write(*,*)'amultlin2 index error',xf,nx,yf,ny
      endif
      f11=data(ix,iy)
      f12=data(ix,iy+1)
      f21=data(ix+1,iy)
      f22=data(ix+1,iy+1)
      d1=f11*(1-xff)+f21*xff
      d2=f12*(1-xff)+f22*xff
      amultilin2=d1*(1.-yff)+d2*yff
      end
c********************************************************************
c******************************************************************
      subroutine parsecmdline(gpar,npar,ngmax,ng,vdmin,vdmax,vw,Te,theta
     $     ,vwfac,vpimax,nim,nmm,nimax,smfac,laspect,lcolor,eldweight
     $     ,vflat,wflat,fflat,lextra,amp,nw,lgrowth,lcoptic,cfilename
     $     ,ltwotone,vrgfac,lthresh,rmitome)
      integer npar,ngmax
      real gpar(npar,ngmax)
      real vdmin,vdmax,vw,Te,theta,vwfac,amp
      logical laspect,lcolor,lextra,lgrowth,lcoptic,ltwotone,lthresh
      character*50 cfilename

      character*50 argument
      logical lgset
      data lgset/.false./
      argument=' '
      eld=1.
      nw=0
      ng=1

c Deal with arguments.
      do ia=1,iargc()
         call getarg(ia,argument)
c Read the input file.
c         write(*,*)argument(1:10)
         if(argument(1:1).ne.'-')then
            if(.not.lcoptic)then
               open(12,file=argument,form='formatted',status='old',err
     $              =101)
               read(12,*,end=102)theta,Te
               write(*,'(a,f8.3,a,f8.3)')'theta=',theta,'  Te=',Te
               ct=cos(theta)
               if(ct.le.0)then
                  write(*,*)'CosTheta is negative. Not allowed.',theta
                  stop
               endif
               read(12,*,end=102)ng
               if(ng.gt.ngmax)then
                  write(*,*)'Too many gaussians'
                  stop
               endif
               vdmin=0.
               vdmax=0.
               write(*,'(a)')'          Den   Velocity   Tperp   Tpara'
               do i=1,ng
                  read(12,*,end=102)(gpar(k,i),k=1,4)
                  write(*,'(a,i1,a,4f8.3)')'gpar(',i,')=',(gpar(k,i),k=1
     $                 ,4)
                  if(gpar(2,i).lt.vdmin)vdmin=ct*gpar(2,i)
                  if(gpar(2,i).gt.vdmax)vdmax=ct*gpar(2,i)
                  vwhere=vwfac*sqrt(sqrt(1-ct**2)*gpar(3,i)
     $                 +ct**2*gpar(4,i))
                  if(vwhere.gt.vw)vw=vwhere
               enddo
            else
c File name of coptic data.
               cfilename=argument
            endif
         else
c Switch arguments.
            if(argument(1:2).eq.'-i')read(argument(3:),*,err=103)vpimax
            if(argument(1:2).eq.'-n')read(argument(3:),*,err=103)nim
            if(argument(1:2).eq.'-m')read(argument(3:),*,err=103)nmm
            if(argument(1:2).eq.'-s')read(argument(3:),*,err=103)smfac
            if(argument(1:2).eq.'-w')read(argument(3:),*,err=103)vwfac
            if(argument(1:2).eq.'-u')read(argument(3:),*,err=103)amp
            if(argument(1:2).eq.'-x')read(argument(3:),*,err=103)vrgfac
            if(argument(1:2).eq.'-p')read(argument(3:),*,err=103)nw
            if(argument(1:2).eq.'-T')read(argument(3:),*,err=103)Te
            if(argument(1:2).eq.'-r')read(argument(3:),*,err=103)rmitome
            if(argument(1:2).eq.'-v')then
               if(lgset)ng=ng+1
               read(argument(3:),*,err=103,end=105)(gpar(k,ng),k=1,4)
               if(gpar(2,ng).gt.vdmax)vdmax=gpar(2,ng)
               if(gpar(2,ng).lt.vdmin)vdmin=gpar(2,ng)
               lgset=.true.
               endif
            if(argument(1:2).eq.'-h')goto 104
            if(argument(1:2).eq.'-?')goto 104
            if(argument(1:2).eq.'-a')laspect=.not.laspect
            if(argument(1:2).eq.'-c')lcolor=.not.lcolor
            if(argument(1:2).eq.'-b')lextra=.not.lextra
            if(argument(1:2).eq.'-g')lgrowth=.not.lgrowth
            if(argument(1:2).eq.'-t')ltwotone=.not.ltwotone
            if(argument(1:2).eq.'-S')lthresh=.not.lthresh
            if(argument(1:2).eq.'-e')read(argument(3:),*,err=103)eld
            if(argument(1:2).eq.'-d')call pfset(-3)
            if(argument(1:7).eq.'-COPTIC')lcoptic=.not.lcoptic
            if(argument(1:2).eq.'-f')then
               if(fflat.eq.0.)then
                  fflat=1.
                  vflat=.5
                  wflat=.1
               endif
               if(argument(1:3).eq.'-fv')read(argument(4:),*,err
     $              =103)vflat
               if(argument(1:3).eq.'-fw')read(argument(4:),*,err
     $              =103)wflat
               if(argument(1:3).eq.'-ff')read(argument(4:),*,err
     $              =103)fflat
c               write(*,*)fflat,vflat,wflat
            endif
         endif
 105     continue
      enddo

      eldweight=eldweight*eld
      if(iargc().le.0)then
         write(*,*)'No input file specified. Using defaults'
      endif
      if(nim.gt.nimax)nim=nimax
      if(abs(nmm).gt.nimax)nmm=nimax
      if(nmm.gt.0)nmm=-nmm

      return

 101  write(*,*)'Failed to open file: ',argument
 103  write(*,*)'Error parsing argument: ',argument
 104  write(*,*)'Usage: chiofv [switches] [filename]'
      write(*,'(2a,2i4,a)')' -n -m imag mesh numbers'
     $     ,' +ve,-ve [',nim,nmm,']'
      write(*,'(a,f7.4,a)')' -i max imaginary vp        [',vpimax,']'
      write(*,'(a,f7.4,a)')' -u uncertainty (noise)     [',amp,']'
      write(*,'(a,f7.4,a)')' -T Te                      [',Te,']'
      write(*,'(a,f7.1,a)')' -r mi/me                   [',rmitome,']'
      write(*,'(a,f7.4,a)')' -x vrange fraction plot    [',vrgfac,']'
      write(*,*)'-s smoothing factor [1]   -w integration width factor'
     $     ,' [4]'
      write(*,*)'-a toggle aspect ratio      [',laspect
     $     ,'] -c toggle coloring          [',lcolor,']'
      write(*,*)'-b toggle extra detail plots[',laspect,']',
     $     ' -g toggle growth rate plot  [',lgrowth,']'
      write(*,*)'-t toggle two tone plot     [',ltwotone,']',
     $     ' -S stability threshold      [',lthresh,']'
      write(*,*)'-d disable display of plots'
      write(*,*)'-e ELD weight factor'
      write(*,*)'-p integer width of smoothing triangle box'
      write(*,'(a,4f7.3,a)')' -v<n>,<v>,<Tpp>,<Tpl> Specify Gaussian ['
     $     ,(gpar(i,ng),i=1,4),']'
      write(*,*)'   First -v overwrites default, succeeding adds new'
      write(*,*)'-f flatspot: -fv center-velocity -fw width'
     $     ,' -ff slope-factor'
      write(*,*)'-COPTIC toggle external distribution reading lcoptic'
      write(*,*)'   subsequent filename refers to external data.'
      call exit
 102  write(*,*)'File ended prematurely',ng,k
      
      end

c*******************************************************************
c Create a flatspot by adding a zero-integral function at 
c position iv, of width iw
      subroutine flatspot(nv,dv,f,iv,iw,fflat)
      real f(0:nv+2)
      integer nv,iv,iw

      if(iw.le.1)return
      if((iv+2*iw).gt.(nv+2).or.((iv-2*iw).lt.0))then
         write(*,*)'Flat spot overlaps ends',iw,iv,nv
         return
      endif

      slope=(f(iv+1)-f(iv-1))/(2.*dv)

c Triangle shaped add-on has analyticity problems.
      if(.false.)then
      do i=iv-2*iw,iv+2*iw
         ir=i-iv
         if(abs(ir).gt.iw)ir=sign(2*iw-abs(ir),ir)
c         write(*,*)i,ir
         f(i)=f(i)-fflat*slope*dv*ir
      enddo
      endif

c Small gaussian derivative.
      do i=0,nv+2
         ir=i-iv
         f(i)=f(i)-fflat*slope*dv*ir*exp(-(float(ir)/iw)**2)
      enddo

      end
c******************************************************************
c This replaces a dummy routine in the accis library. 
c Stores the just completed contour of length imax and value cv.
!      subroutine acpathdoc(z,cv,L,imax,xc,yc,i)
      subroutine acpathdoc(imax,xc,yc,cv)
      integer imax
      real cv,xc(imax),yc(imax)
c common block for communications includes iacp...
      include 'acpathcom.f'
      if(iacpsw.ne.0)then
c Switched on:
         if(imax.ne.0.)then
c Store onto the contour stack.
            do i=1,imax
               if(iacppt.lt.iacpmax)then
                  iacppt=iacppt+1
                  xacp(iacppt)=xc(i)
                  yacp(iacppt)=yc(i)
c                  if(i.eq.1)write(*,'(a,3i4,3f9.3)')'acpathdoc',
c     $        iacppt,iacpcp,iacpcon(iacpcp),cv,xacp(iacppt),yacp(iacppt)
               else
                  write(*,*)'Exhausted iacpmax length',iacpmax
                  return
               endif
            enddo
c Terminate a contour and store its end.
            if(iacpcp.lt.iacpconmax)then
               iacpcp=iacpcp+1
               iacpcon(iacpcp)=iacppt
            else
               write(*,*)'Exhausted iacpconmax length',iacpconmax
               return
            endif
         endif
c         write(*,'(a,3i4,3f9.3)')'acpathdoc',
c     $        iacppt,iacpcp,iacpcon(iacpcp),cv,xacp(iacppt),yacp(iacppt)
      else ! Reset call
         iacpcp=1
         iacppt=0
      endif
c iacpcon is the pointer to start of contour which is 
c xacp(iacpcon(i)+1),yacp(iacpcon(i)+1) 
c iacpcon(i)+1 is the start index of each contour segment.
c iacpcon(i+1)-iacpcon(id) is the length of contour i.
c iacpcp is the current contour number, and when completed
c is 1+ the total number of contours. iacpcon(iacpcp) points to 
c the end of the last.
      end
c**********************************************************************
      subroutine addnoise(nv,f,amp)
      integer nv
      real f(0:nv+2)

c Initialize random number generator.
      df=ran1(-1)
      do i=3,nv-1
         df=ran1(0)
c         write(*,*)i,' ran1=',df
         f(i)=f(i)+df*amp
      enddo

      end
c**********************************************************************
      FUNCTION ran1(IDUM)
c Returns a uniform random deviate between 0 and 1.
c Reentrant version of ran1 makes the state visible in common.
c Get explicit 
      implicit none
      real ran1
      integer idum
c Common here requires separation of processes for reentrancy.
c So this is not thread safe.
c      include 'ran1com.f'
c Internal state of the fortran random number generator ran1 is made
c visible through this common block. It can then be saved by just
c writing ranstate(1:103) and restored by reading it back.
c This requires integers to be no longer than reals.
      real ranstate(103)
c Internal state of ran1:
      real rrnd(97)
      integer irx1,irx2,irx3,jrc
      equivalence (rrnd,ranstate)
      equivalence (irx1,ranstate(98)),(irx2,ranstate(99)),
     $     (irx3,ranstate(100)),(jrc,ranstate(101))
c Internal state of the gaussian random generator gasdev:
      integer gd_iset
      real gd_gset
      equivalence (gd_iset,ranstate(102)),(gd_gset,ranstate(103))
c The whole thing:
      common /ran1com/ranstate

      real RM1,RM2
      integer M1,M2,M3,ic1,ic2,ic3,ia1,ia2,ia3,j
      save
c      DIMENSION Rrnd(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
c Can't do auto initialize. Have to do it by hand.
c      DATA IFF /0/
c      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
      IF (IDUM.LT.0) THEN
c        IFF=1
        IrX1=MOD(IC1-IDUM,M1)
        IrX1=MOD(IA1*IrX1+IC1,M1)
        IrX2=MOD(IrX1,M2)
        IrX1=MOD(IA1*IrX1+IC1,M1)
        IrX3=MOD(IrX1,M3)
        do J=1,97
           IrX1=MOD(IA1*IrX1+IC1,M1)
           IrX2=MOD(IA2*IrX2+IC2,M2)
           Rrnd(J)=(FLOAT(IrX1)+FLOAT(IrX2)*RM2)*RM1
        enddo
c Tell gasdev to reset.
        gd_iset=0
      ENDIF
      IrX1=MOD(IA1*IrX1+IC1,M1)
      IrX2=MOD(IA2*IrX2+IC2,M2)
      IrX3=MOD(IA3*IrX3+IC3,M3)
      JRC=1+(97*IrX3)/M3
      if(JRC.GT.97.OR.JRC.LT.1)then
         write(*,*)'RAN1 Error!'
      endif
      RAN1=Rrnd(JRC)
      Rrnd(JRC)=(FLOAT(IrX1)+FLOAT(IrX2)*RM2)*RM1
      RETURN
      END
c*****************************************************************
      subroutine trismooth(n,nw,v,ww,vs)
c Smooth a vector using a triangular box.
c So vs= 1/(nw+1) sum_{-nw}^{+nw} v_j*(1-j/(nw+1))
c On entry
c   n is the length of v
c   nw is the smoothing length in samples (0 => do nothing)
c   v is the input vector
c   ww is a work space of length n
c On exit
c   vs is the smoothed result. 
c The arrays v and vs or the arrays v and ww, but not ww and vs, 
c can be the same storage locations.
c If the vector is on a non-uniform grid, then this smoothing is not
c area conserving, because it knows nothing about the grid spacing. 
      integer n,nw
      real v(n),vs(n),ww(n)
      do i=1,n
         ww(i)=v(i)
      enddo
      p=nw+1.
      do i=1,n
         vi=0.
         do j=-nw,nw
c Deal with the ends...
            k=min(n,max(1,i+j))
            vi=vi+ww(k)*(1.-abs(float(j))/p)
         enddo
         vs(i)=vi/p
      enddo
      end
c*******************************************************************
c Test the trismooth routine on a step function.
      subroutine trismoothtest(nw)
      integer nw,nt
      parameter (nt=200)
      real v(nt),vs(nt),ww(nt),x(nt)

      do i=1,nt
         x(i)=i
         if(i.lt.nt/2)then
            v(i)=1.
         else
            v(i)=0.
         endif
      enddo
      call trismooth(nt,nw,v,ww,vs)
c      write(*,*)v,vs
      call autoplot(x,v,nt)
      call color(1)
      call polyline(x,vs,nt)
      call pltend()
      end
c**********************************************************************
c Test the COPTIC input data reading
      subroutine COPTICverif(nv,dv,v0,f)
      integer nv
      real dv,v0,f(0:nv+2)
      open(14,file='COPTICverif.dat',status='unknown',err=101)
      close(14,status='delete')
      open(14,file='COPTICverif.dat',status='new',err=101)
          write(14,*) nv
          write(14,*) dv
          write(14,*) v0
          write(14,'(1g14.6)')(f(i),i=0,nv+2)
      close(14)
      goto 102
 101  write(*,*)'Error opening file'
      close(14,status='delete')
      call exit(1)
 102  continue
      end
c*********************************************************************
      subroutine stationary(nv,f,ngmax,fstat,nstat)
      real f(nv),fstat(ngmax)
! Finds stationary values of f and return nstat values in fstat
      nstat=0
      do i=2,nv-1
      if((f(i+1)-f(i))*(f(i)-f(i-1)).lt.0..or.f(i+1)-f(i).eq.0)then
            nstat=nstat+1
            fstat(nstat)=f(i)
         endif
         if(nstat.eq.ngmax)return
      enddo
      end
