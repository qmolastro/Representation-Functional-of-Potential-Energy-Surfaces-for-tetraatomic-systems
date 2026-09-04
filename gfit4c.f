********************************************************************************
C     PROGRAM TO FIT A FOUR-BODY(6D) G(lobal)P(otential)E(energy)S(urface)
C
C     REFERENCES:
C     A. Aguado, M. Paniagua, J. Chem. Phys. 96(1992) 1265-1275.
C
C     A. Aguado, C. Suarez, M. Paniagua, J. Chem. Phys. 101(1994) 4004-4010.
C
C     A. Aguado, C. Tablero, M. Paniagua, Comput. Phys. Commun. 108(1998) 
C                                         259-266.
C
C     A. Aguado, C. Tablero, M. Paniagua, Comput. Phys. Commun. 134(2001)
C                                         97-109. (One of the paper corresponding to
C                                         the documentation of this program).
C
C     A. Aguado, C. Tablero, M. Paniagua, Comput. Phys. Commun. 140(2001)
C                                         412-417. (One of the paper corresponding to
C                                         the documentation of this program).
********************************************************************************
      implicit real*8 (a-h,o-z)
      external cal1,cal2,cal3,cal4,cal6

      include 'dimensions.inc'
      
      parameter (bfact=627.51d0)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      common/ipot41/ipa41(nexp41,7),nfas41(igrad4)
      common/ipot42/ipa42(nexp42,7),nfas42(igrad4)
      common/ipot43/ipa43(nexp43,7),nfas43(igrad4)
      common/ipot44/ipa44(nexp44,7),nfas44(igrad4)
      common/ipot45/ipa45(nexp45,6),nfas45(igrad4)
      dimension vex12(2),vex13(2),vex14(2),vex23(2),vex24(2),vex34(2),
     ,          vex123(3),vex124(3),vex134(3),vex234(3),
c    ,          vex1,vex2(2),vex3(3),vex4(4),vex6(6)
     ,          vex2(2),vex3(3),vex4(4),vex6(6)
      dimension c123(nexp3),c124(nexp3),c134(nexp3),
     ,          c234(nexp3),c12(nexp2),c13(nexp2),c14(nexp2),
     ,          c23(nexp2),c24(nexp2),c34(nexp2)
      dimension aux1(1),aux2(2),aux3(3),aux4(4),aux6(6)
      character*6 nombre
*************************************************************************
c     the input is read from the file gfit4c.inp
c     for a detailed description see the documentation
C
c     molecula = 2(diatomic molecule), 3(triatomic m.), 4(tetraatomic m.)
*************************************************************************
      ifor1=0
 1015 format(/,1x,'tetratomic fit ',a6)
 1010 format(1x,'iteration',6x,'exponent ',8x,'rms(u.a.)',
     ,       3x,'rms(kcal)',3x,'e max (kcal)'/)
 1020 format(1x,'rms=dsqrt(sum(del**2)/n)= ',f12.6,' u.a. ',1x,f12.4,
     - ' kcal/mol',1x,' emax = ', f12.6,' kcal/mol')
  417 format(1x,'             x               v-inp        v-fit   ',
     *          '     diff(u.a.)  diff(Kcal)',/)
*************************************************************************
      open (unit=7,status='unknown',file='fitout.f')
      read(5,*) indice,ifor,numiter,lim
      if (lim.le.0) lim=500
      call cabeza(indice)
      if (indice.eq.1) then
         it123=1
         it124=1
         it134=1
         it234=1
      endif
      if (indice.eq.2) then
         it123=2
         it124=2
         it134=2
         it234=1
      endif
      if (indice.eq.3) then
         it123=2
         it124=2
         it134=2
         it234=2
      endif
      if (indice.eq.4) then
         it123=3
         it124=3
         it134=2
         it234=2
      endif
      if (indice.eq.5) then
         it123=3
         it124=3
         it134=3
         it234=3
      endif
*************************************************************************
c     DIATOMICS
*************************************************************************
      call ajus2('diat12',mt12,vex12,c12)
*************************************************************************
      if (indice.eq.3.or.indice.eq.4.or.indice.eq.5) then
         call ajus2('diat13',mt13,vex13,c13)
        else
         mt13=mt12
         vex13(1)=vex12(1)
         vex13(2)=vex12(2)
         call digual(c13,nexp2,c12,nexp2,mt13)
      endif
*************************************************************************
      if (indice.eq.5) then
         call ajus2('diat14',mt14,vex14,c14)
        else
         mt14=mt13
         vex14(1)=vex13(1)
         vex14(2)=vex13(2)
         call digual(c14,nexp2,c13,nexp2,mt14)
      endif
*************************************************************************
      if (indice.eq.2.or.indice.eq.4.or.indice.eq.5) then
         call ajus2('diat23',mt23,vex23,c23)
        else
         mt23=mt14
         vex23(1)=vex14(1)
         vex23(2)=vex14(2)
         call digual(c23,nexp2,c14,nexp2,mt23)
      endif
*************************************************************************
      if (indice.eq.5) then
         call ajus2('diat24',mt24,vex24,c24)
        else
         mt24=mt23
         vex24(1)=vex23(1)
         vex24(2)=vex23(2)
         call digual(c24,nexp2,c23,nexp2,mt24)
      endif
*************************************************************************
      if (indice.eq.3.or.indice.eq.4.or.indice.eq.5) then
         call ajus2('diat34',mt34,vex34,c34)
        else
         mt34=mt24
         vex34(1)=vex24(1)
         vex34(2)=vex24(2)
         call digual(c34,nexp2,c24,nexp2,mt34)
      endif
*************************************************************************
c     TRIATOMICS
*************************************************************************
      call ajus3('tri123',it123,mt12,mt13,mt23,vex12,vex13,vex23,
     ,           c12,c13,c23,mt123,vex123,c123)
*************************************************************************
      if (indice.eq.5) then
         call ajus3('tri124',it124,mt12,mt14,mt24,vex12,vex14,vex24,
     ,              c12,c14,c24,mt124,vex124,c124)
       else
         mt124=mt123
         vex124(1)=vex123(1)
         vex124(2)=vex123(2)
         vex124(3)=vex123(3)
         call digual(c124,nexp3,c123,nexp3,nexp3)
      endif
*************************************************************************
      if (indice.eq.3.or.indice.eq.4.or.indice.eq.5) then
         call ajus3('tri134',it134,mt13,mt14,mt34,vex13,vex14,vex34,
     ,              c13,c14,c34,mt134,vex134,c134)
       else
         mt134=mt124
         vex134(1)=vex124(1)
         vex134(2)=vex124(2)
         vex134(3)=vex124(3)
         call digual(c134,nexp3,c124,nexp3,nexp3)
      endif
*************************************************************************
      if (indice.eq.2.or.indice.eq.4.or.indice.eq.5) then
         call ajus3('tri234',it234,mt23,mt24,mt34,vex23,vex24,vex34,
     ,              c23,c24,c34,mt234,vex234,c234)
       else
         mt234=mt134
         vex234(1)=vex134(1)
         vex234(2)=vex134(2)
         vex234(3)=vex134(3)
         call digual(c234,nexp3,c134,nexp3,nexp3)
      endif
*************************************************************************
c     TETRAATOMIC
*************************************************************************
      molecula=4
      npar=indice
      write(6,1015) 'tt1234'
      if (npar.eq.1) then
         read(5,*) n,mt,vex1,e01234
         print*,'input parameters ',n,mt,vex1,e01234
         if(mt.gt.12) stop 'maximum degree = 12'
         nombre='tfaaaa'
      endif
      if (npar.eq.2) then
         read(5,*) n,mt,vex2(1),vex2(2),e01234
         print*,'input parameters ', n,mt,vex2(1),vex2(2),e01234
         if(mt.gt.9) stop 'maximum degree = 9'
         nombre='tfabbb'
      endif
      if (npar.eq.3) then
         read(5,*) n,mt,vex3(1),vex3(2),vex3(3),e01234
         print*,'input parameters ', n,mt,vex3(1),vex3(2),vex3(3),e01234
         if(mt.gt.8) stop 'maximum degree = 8'
         nombre='tfaabb'
      endif
      if (npar.eq.4) then
         read(5,*) n,mt,vex4(1),vex4(2),vex4(3),vex4(4),e01234
         print*,'input parameters ', n,mt,vex4(1),vex4(2),vex4(3),
     *                                    vex4(4),e01234
         nombre='tfabcc'
         if(mt.gt.7) stop 'maximum degree = 7'
      endif
      if (npar.eq.5) then
         read(5,*) n,mt,vex6(1),vex6(2),vex6(3),vex6(4),vex6(5),vex6(6),
     ,             e01234
         print*, 'input parameters ',n,mt,vex6(1),vex6(2),vex6(3),
     ,           vex6(4),vex6(5),vex6(6), e01234
         if(mt.gt.6) stop 'maximum degree = 6'
         nombre='tfabcd'
      endif
      if (n.gt.nmax) stop 100
      if (mt.gt.n) stop 200
      do 11 i=1,n
         read(5,*) r12(i),r13(i),r14(i),r23(i),r24(i),r34(i),xener(i)
         call diat(mt12,vex12,r12(i),c12,0.d0,e12)
         call diat(mt13,vex13,r13(i),c13,0.d0,e13)
         call diat(mt14,vex14,r14(i),c14,0.d0,e14)
         call diat(mt23,vex23,r23(i),c23,0.d0,e23)
         call diat(mt24,vex24,r24(i),c24,0.d0,e24)
         call diat(mt34,vex34,r34(i),c34,0.d0,e34)
         if (indice.ne.3) then
           call triat(it123,mt123,r12(i),r13(i),r23(i),vex123,c123,e123)
           else
           call triat(it123,mt123,r13(i),r23(i),r12(i),vex123,c123,e123)
         endif
         if (indice.ne.3) then
           call triat(it124,mt124,r12(i),r14(i),r24(i),vex124,c124,e124)
           else
           call triat(it124,mt124,r14(i),r24(i),r12(i),vex124,c124,e124)
         endif
         call triat(it134,mt134,r13(i),r14(i),r34(i),vex134,c134,e134)
         call triat(it234,mt234,r23(i),r24(i),r34(i),vex234,c234,e234)
         v123=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e01234

         xener(i)=xener(i)-(e12+e13+e14+e23+e24+e34+e123+e124+e134+e234)
     ,            -e01234
11    continue
      iter = 0
      if (ifor.eq.1.and.numiter.le.2) write(6,1010)
      if (npar.eq.1) then
         nfas=nfas41(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.2) call minim(cal1,vex1,rmsf,1,lim,aux1,ier)
         vex1=dabs(vex1)
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.2) bux=cal1(vex1)
         ifor1=0
         write(6,*) vex1
         call saaaa('aaaa',n,mt,vex1,rms,emax)
      endif
      if (npar.eq.2) then
         nfas=nfas42(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.2) call minim(cal2,vex2,rmsf,2,lim,aux2,ier)
         vex2(1)=dabs(vex2(1))
         vex2(2)=dabs(vex2(2))
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.2) bux=cal2(vex2)
         ifor1=0
         write(6,*) vex2(1),vex2(2)
         call sabbb('abbb',n,mt,vex2,rms,emax)
      endif
      if (npar.eq.3) then
         nfas=nfas43(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.2) call minim(cal3,vex3,rmsf,3,lim,aux3,ier)
         vex3(1)=dabs(vex3(1))
         vex3(2)=dabs(vex3(2))
         vex3(3)=dabs(vex3(3))
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.2) bux=cal3(vex3)
         ifor1=0
         write(6,*) vex3(1),vex3(2),vex3(3)
         call saabb('aabb',n,mt,vex3,rms,emax)
      endif
      if (npar.eq.4) then
         nfas=nfas44(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.2) call minim(cal4,vex4,rmsf,4,lim,aux4,ier)
         vex4(1)=dabs(vex4(1))
         vex4(2)=dabs(vex4(2))
         vex4(3)=dabs(vex4(3))
         vex4(4)=dabs(vex4(4))
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.2) bux=cal4(vex4)
         ifor1=0
         write(6,*) vex4(1),vex4(2),vex4(3),vex4(4)
         call sabcc('abcc',n,mt,vex4,rms,emax)
      endif
      if (npar.eq.5) then
         nfas=nfas45(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.2) call minim(cal6,vex6,rmsf,6,lim,aux6,ier)
         vex6(1)=dabs(vex6(1))
         vex6(2)=dabs(vex6(2))
         vex6(3)=dabs(vex6(3))
         vex6(4)=dabs(vex6(4))
         vex6(5)=dabs(vex6(5))
         vex6(6)=dabs(vex6(6))
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.2) bux=cal6(vex6)
         ifor1=0
         write(6,*) vex6(1),vex6(2),vex6(3),vex6(4),vex6(5),vex6(6)
         call sabcd('abcd',n,mt,vex6,rms,emax)
      endif
      write(6,1020) rms,rms*bfact,emax*bfact
      close(7)
      end
*************************************************************************
************************************************************************
*                                                                      *
*     general case                                                     *
*                                                                      *
************************************************************************
************************************************************************
      function cal1(vex1)
************************************************************************
      implicit  real * 8 (a-h,o-z)
      parameter (bfact=627.51d0)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      vex1=dabs(vex1)
      iter = iter + 1
      if (molecula.eq.3) call gfaaa(n,mt,vex1)
      if (molecula.eq.4) call gfaaaa(n,mt,vex1)
      call calmat
      if (ifor.eq.1) write(6,772) iter, vex1,rms, rms*bfact, emax*bfact
  772 format(1x,i5,5x,4(f10.5,3x))
      cal1=rms
      return
      end
************************************************************************
      function cal2(vex2)
************************************************************************
      implicit  real * 8 (a-h,o-z)
      parameter (bfact=627.51d0)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      dimension vex2(2)
      vex2(1)=dabs(vex2(1))
      vex2(2)=dabs(vex2(2))
      iter = iter + 1
      if (molecula.eq.2) call gfab(n,mt,vex2)
      if (molecula.eq.3) call gfabb(n,mt,vex2)
      if (molecula.eq.4) call gfabbb(n,mt,vex2)
      call calmat
      if (ifor.eq.1) write(6,772) iter, vex2(1),vex2(2),rms, rms*bfact,
     ,emax*bfact
  772 format(1x,i5,5x,5(f10.5,3x))
      cal2=rms
      return
      end
************************************************************************
      function cal3(vex3)
************************************************************************
      implicit  real * 8 (a-h,o-z)
      parameter (bfact=627.51d0)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      dimension vex3(3)
      vex3(1)=dabs(vex3(1))
      vex3(2)=dabs(vex3(2))
      vex3(3)=dabs(vex3(3))
      iter = iter + 1
      if (molecula.eq.3) call gfabc(n,mt,vex3)
      if (molecula.eq.4) call gfaabb(n,mt,vex3)
      call calmat
      if (ifor.eq.1) write(6,772) iter,vex3(1),vex3(2),vex3(3), rms,
     , rms*bfact,emax*bfact
  772 format(1x,i5,5x,6(f10.5,3x))
      cal3=rms
      return
      end
************************************************************************
      function cal4(vex4)
************************************************************************
      implicit  real * 8 (a-h,o-z)
      parameter (bfact=627.51d0)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      dimension vex4(4)
      vex4(1)=dabs(vex4(1))
      vex4(2)=dabs(vex4(2))
      vex4(3)=dabs(vex4(3))
      vex4(4)=dabs(vex4(4))
      iter = iter + 1
      if (molecula.eq.4) call gfabcc(n,mt,vex4)
      call calmat
      if (ifor.eq.1) write(6,772) iter,vex4(1),vex4(2),vex4(3),
     , vex4(4),rms,rms*bfact,emax*bfact
  772 format(1x,i5,5x,7(f10.5,3x))
      cal4=rms
      return
      end
************************************************************************
      function cal6(vex6)
************************************************************************
      implicit  real * 8 (a-h,o-z)
      parameter (bfact=627.51d0)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      dimension vex6(6)
      vex6(1)=dabs(vex6(1))
      vex6(2)=dabs(vex6(2))
      vex6(3)=dabs(vex6(3))
      vex6(4)=dabs(vex6(4))
      vex6(5)=dabs(vex6(5))
      vex6(6)=dabs(vex6(6))
      iter = iter + 1
      if (molecula.eq.4) call gfabcd(n,mt,vex6)
      call calmat
      if (ifor.eq.1) write(6,772) iter,vex6(1),vex6(2),vex6(3),
     , vex6(4),vex6(5),vex6(6),rms,rms*bfact,emax*bfact
  772 format(1x,i5,5x,9(f10.5,3x))
      cal6=rms
      return
      end
************************************************************************
      subroutine calmat
************************************************************************
c     computes the needed matrices
************************************************************************
      implicit  real * 8 (a-h,o-z)
      include 'dimensions.inc'
      parameter (bfact=627.51d0)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      do 20 k=1,nfas
         sm(k) = 0.0d0
         do 20 i=1,n
            sm(k) = sm(k) + xener(i)*xx(k,i)
   20 continue
      ind = 0
      do 50 k=1,nfas
      do 50 l=1,nfas
         ind = ind + 1
         sg(ind) = 0.0d0
         do 60 i=1,n
            sg(ind) = sg(ind) + xx(k,i)*xx(l,i)
   60    continue
   50 continue
      call gauss(sg,sm,nfas,ks)

c     if(ks.ne.0) stop 400
      if(ks.ne.0) stop 'There are redundancies in the data or functions'
      sum = 0.0d0
      emax =0.d0
      do 7 i=1,n
         v = 0.0d0
         do 8 k=1,nfas
    8       v = v + sm(k)*xx(k,i)
         vdi = v - xener(i)
         if(ifor1.eq.1) write(6,777) i,r12(i),r13(i),r23(i),
     ,          r14(i),r24(i),r34(i),xener(i),v,vdi,vdi*bfact
         if (dabs(vdi).gt.emax) emax = dabs(vdi)
         sum = sum + vdi**2
    7 continue
  777 format(i3,1x,6(1x,f7.3),4f13.8)
      rms = dsqrt(sum/n)
      sdv = dsqrt(sum/(n-nfas))
      return
      end
*********************************************************************
      subroutine cabeza(i)
      implicit  real * 8 (a-h,o-z)
      if (i.eq.1) write(7,1)
      if (i.eq.1) write(7,6)
      if (i.eq.2) write(7,2)
      if (i.eq.2) write(7,6)
      if (i.eq.3) write(7,3)
      if (i.eq.3) write(7,13)
      if (i.eq.4) write(7,4)
      if (i.eq.4) write(7,6)
      if (i.eq.5) write(7,5)
      if (i.eq.5) write(7,6)

1     format(
     -6x,'subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -6x,'implicit  real * 8 (a-h,o-z)',/,
     -6x,'dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)',/,
     -6x,'call diat12(r12,e12,d12,iop)',/,
     -6x,'call diat12(r13,e13,d13,iop)',/,
     -6x,'call diat12(r14,e14,d14,iop)',/,
     -6x,'call diat12(r23,e23,d23,iop)',/,
     -6x,'call diat12(r24,e24,d24,iop)',/,
     -6x,'call diat12(r34,e34,d34,iop)',/,
     -6x,'call tri123(r12,r13,r23,e123,d123,iop)',/,
     -6x,'call tri123(r12,r14,r24,e124,d124,iop)',/,
     -6x,'call tri123(r13,r14,r34,e134,d134,iop)',/,
     -6x,'call tri123(r23,r24,r34,e234,d234,iop)',/,
     -6x,'call aaaa(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)',/,
     -6x,'ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234')
2     format(
     -6x,'subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -6x,'implicit  real * 8 (a-h,o-z)',/,
     -6x,'dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)',/,
     -6x,'call diat12(r12,e12,d12,iop)',/,
     -6x,'call diat12(r13,e13,d13,iop)',/,
     -6x,'call diat12(r14,e14,d14,iop)',/,
     -6x,'call diat23(r23,e23,d23,iop)',/,
     -6x,'call diat23(r24,e24,d24,iop)',/,
     -6x,'call diat23(r34,e34,d34,iop)',/,
     -6x,'call tri123(r12,r13,r23,e123,d123,iop)',/,
     -6x,'call tri123(r12,r14,r24,e124,d124,iop)',/,
     -6x,'call tri123(r13,r14,r34,e134,d134,iop)',/,
     -6x,'call tri234(r23,r24,r34,e234,d234,iop)',/,
     -6x,'call abbb(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)',/,
     -6x,'ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234')
3     format(
     -6x,'subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -6x,'implicit  real * 8 (a-h,o-z)',/,
     -6x,'dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)',/,
     -6x,'call diat12(r12,e12,d12,iop)',/,
     -6x,'call diat13(r13,e13,d13,iop)',/,
     -6x,'call diat13(r14,e14,d14,iop)',/,
     -6x,'call diat13(r23,e23,d23,iop)',/,
     -6x,'call diat13(r24,e24,d24,iop)',/,
     -6x,'call diat34(r34,e34,d34,iop)',/,
     -6x,'call tri123(r13,r23,r12,e123,d123,iop)',/,
     -6x,'call tri123(r14,r24,r12,e124,d124,iop)',/,
     -6x,'call tri134(r13,r14,r34,e134,d134,iop)',/,
     -6x,'call tri134(r23,r24,r34,e234,d234,iop)',/,
     -6x,'call aabb(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)',/,
     -6x,'ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234')
4     format(
     -6x,'subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -6x,'implicit  real * 8 (a-h,o-z)',/,
     -6x,'dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)',/,
     -6x,'call diat12(r12,e12,d12,iop)',/,
     -6x,'call diat13(r13,e13,d13,iop)',/,
     -6x,'call diat13(r14,e14,d14,iop)',/,
     -6x,'call diat23(r23,e23,d23,iop)',/,
     -6x,'call diat23(r24,e24,d24,iop)',/,
     -6x,'call diat34(r34,e34,d34,iop)',/,
     -6x,'call tri123(r12,r13,r23,e123,d123,iop)',/,
     -6x,'call tri123(r12,r14,r24,e124,d124,iop)',/,
     -6x,'call tri134(r13,r14,r34,e134,d134,iop)',/,
     -6x,'call tri234(r23,r24,r34,e234,d234,iop)',/,
     -6x,'call abcc(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)',/,
     -6x,'ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234')
5     format(
     -6x,'subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -6x,'implicit  real * 8 (a-h,o-z)',/,
     -6x,'dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)',/,
     -6x,'call diat12(r12,e12,d12,iop)',/,
     -6x,'call diat13(r13,e13,d13,iop)',/,
     -6x,'call diat14(r14,e14,d14,iop)',/,
     -6x,'call diat23(r23,e23,d23,iop)',/,
     -6x,'call diat24(r24,e24,d24,iop)',/,
     -6x,'call diat34(r34,e34,d34,iop)',/,
     -6x,'call tri123(r12,r13,r23,e123,d123,iop)',/,
     -6x,'call tri124(r12,r14,r24,e124,d124,iop)',/,
     -6x,'call tri134(r13,r14,r34,e134,d134,iop)',/,
     -6x,'call tri234(r23,r24,r34,e234,d234,iop)',/,
     -6x,'call abcd(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)',/,
     -6x,'ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234')
6     format(
     -6x,'der(1)=d12+d123(1)+d124(1)+d1234(1) ',/,
     -6x,'der(2)=d13+d123(2)+d134(1)+d1234(2) ',/,
     -6x,'der(3)=d14+d124(2)+d134(2)+d1234(3) ',/,
     -6x,'der(4)=d23+d123(3)+d234(1)+d1234(4) ',/,
     -6x,'der(5)=d24+d124(3)+d234(2)+d1234(5) ',/,
     -6x,'der(6)=d34+d134(3)+d234(3)+d1234(6) ',/,
     -6x,'return',/,
     -6x,'end')
13    format(
     -6x,'der(1)=d12+d123(3)+d124(3)+d1234(1) ',/,
     -6x,'der(2)=d13+d123(1)+d134(1)+d1234(2) ',/,
     -6x,'der(3)=d14+d124(1)+d134(2)+d1234(3) ',/,
     -6x,'der(4)=d23+d123(2)+d234(1)+d1234(4) ',/,
     -6x,'der(5)=d24+d124(2)+d234(2)+d1234(5) ',/,
     -6x,'der(6)=d34+d134(3)+d234(3)+d1234(6) ',/,
     -6x,'return',/,
     -6x,'end')
      return
      end
*********************************************************************
      subroutine digual(a,na,b,nb,mg)
      implicit real*8 (a-h,o-z)
      dimension a(na),b(nb)
      do 1 i=1,mg
         a(i)=b(i)
1     continue
      return
      end
*************************************************************************
      subroutine triat(npar,mt,r12,r13,r23,vec3,cf,ener)
      implicit real*8 (a-h,o-z)
      include 'dimensions.inc'
      dimension vec2(2),vec3(3),cf(nexp3)
      vec1=vec3(1)
      vec2(1)=vec3(1)
      vec2(2)=vec3(2)
      if (npar.eq.1) call triaaa(mt,r12,r13,r23,vec1,cf,ener)
      if (npar.eq.2) call triabb(mt,r12,r13,r23,vec2,cf,ener)
      if (npar.eq.3) call triabc(mt,r12,r13,r23,vec3,cf,ener)
      return
      end
*************************************************************************
      subroutine ajus3(nombre,nupar,mt12,mt13,mt23,vex12,vex13,vex23,
     ,                 c12,c13,c23,mt123,vex123,c123)
*************************************************************************
      implicit real*8 (a-h,o-z)
      external cal1,cal2,cal3
      character*6 nombre
      
      include 'dimensions.inc'
      
      parameter (bfact=627.51d0)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      dimension vex12(2),vex13(2),vex23(2)
      dimension c12(nexp2),c13(nexp2),c23(nexp2)
      dimension vex123(3),c123(nexp3)
      dimension vex2(2),vex3(3)
      dimension aux1(1),aux2(2),aux3(3)
*************************************************************************
      lim = 500
 1015 format(/,1x,'triatomics fit ',a6)
 1010 format(1x,'iteration',6x,'exponent ',8x,'rms(u.a.)',
     ,       3x,'rms(kcal)',3x,'e max (kcal)'/)
 1020 format(1x,'rms=dsqrt(sum(del**2)/n)= ',f12.6,' u.a. ',1x,f12.4,
     - ' kcal/mol',1x,' emax = ', f12.6,' kcal/mol')
  417 format(1x,'             x               v-inp        v-fit   ',
     *          '     diff(u.a.)  diff(Kcal)',/)
      molecula=3
      npar=nupar
      write(6,1015) nombre
      if (npar.eq.1) then
         read(5,*) n,mt,vex1,e0123
         write(6,*) 'input parameters ',n,mt,vex1,e0123
      endif
      if (npar.eq.2) then
         read(5,*) n,mt,vex2(1),vex2(2),e0123
         write(6,*) 'input parameters ',n,mt,vex2(1),vex2(2),e0123
      endif
      if (npar.eq.3) then
         read(5,*) n,mt,vex3(1),vex3(2),vex3(3),e0123
         write(6,*) 'input parameters ',n,mt,vex3(1),vex3(2),vex3(3),e0123
      endif
      if (n.gt.nmax) stop 'number of points(n) greater than nmax'
      if (mt.gt.n) stop 'polynomial order greater than number of
     , points'
      do 7 i=1,n
         read(5,*) r12(i),r13(i),r23(i),xener(i)
         call diat(mt12,vex12,r12(i),c12,0.d0,e12)
         call diat(mt13,vex13,r13(i),c13,0.d0,e13)
         call diat(mt23,vex23,r23(i),c23,0.d0,e23)
         xener(i)=xener(i)-(e12+e13+e23)-e0123
7     continue
      mt123=mt
      identi=123
      iter = 0
      if (ifor.eq.1.and.numiter.le.1) write(6,1010)
      if (npar.eq.1) then
c         call gifaaa
         nfas=nfas31(mt)
         print*,'number of SAF= ',nfas
         if (numiter.le.1) call minim(cal1,vex1,rmsf,1,lim,aux1,ier)
         vex1=dabs(vex1)
         if (ifor.eq.1) write(6,417)
         ifor1=ifor
         if (ifor.eq.1.or.numiter.gt.1) bux=cal1(vex1)
         ifor1=0
         write(6,*) vex1
         call saaa(nombre,n,mt,vex1,rms,emax)
         vex123(1)=vex1
         vex123(2)=0.d0
         vex123(3)=0.d0
        else
         if (npar.eq.2) then
c            call gifabb
            nfas=nfas32(mt)
            print*,'number of SAF= ',nfas
            if (numiter.le.1) call minim(cal2,vex2,rmsf,2,lim,aux2,ier)
            vex2(1)=dabs(vex2(1))
            vex2(2)=dabs(vex2(2))
            if (ifor.eq.1) write(6,417)
            ifor1=ifor
            if (ifor.eq.1.or.numiter.gt.1) bux=cal2(vex2)
            ifor1=0
            write(6,*) vex2(1),vex2(2)
            call sabb(nombre,n,mt,vex2,rms,emax)
            vex123(1)=vex2(1)
            vex123(2)=vex2(2)
            vex123(3)=0.d0
           else
c            call gifabc
            nfas=nfas33(mt)
            print*,'number of SAF= ',nfas
            if (numiter.le.1) call minim(cal3,vex3,rmsf,3,lim,aux3,ier)
            vex3(1)=dabs(vex3(1))
            vex3(2)=dabs(vex3(2))
            vex3(3)=dabs(vex3(3))
            if (ifor.eq.1) write(6,417)
            ifor1=ifor
            if (ifor.eq.1.or.numiter.gt.1) bux=cal3(vex3)
            ifor1=0
            write(6,*) vex3(1),vex3(2),vex3(3)
            call sabc(nombre,n,mt,vex3,rms,emax)
            vex123(1)=vex3(1)
            vex123(2)=vex3(2)
            vex123(3)=vex3(3)
         endif
      endif
      write(6,1020) rms,rms*bfact,emax*bfact
      call digual(c123,nexp3,sm,nexp,nexp3)
      return
      end
*************************************************************************
      subroutine ajus2(nombre,mt12,vex12,c12)
*************************************************************************
c     needs "nombre" and returns "mt12,vex12,c12"
*************************************************************************
      implicit real*8 (a-h,o-z)
      external cal2
      character*6 nombre
      include 'dimensions.inc'
      parameter (bfact=627.51d0)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax)
      common/info/molecula,npar,mt,n,nfas,iter,ifor,ifor1,numiter
      common/rms/rms,emax,sdv
      dimension vex12(2),c12(nexp2)
      dimension aux2(2)
*************************************************************************
      lim = 500
 1015 format(/,1x,'diatomics fit ',a6)
 1010 format(1x,'iteration',6x,'exponent ',8x,'rms(u.a.)',
     ,       3x,'rms(kcal)',3x,'e max (kcal)'/)
 1020 format(1x,'rms=dsqrt(sum(del**2)/n)= ',f12.6,' u.a. ',1x,f12.4,
     - ' kcal/mol',1x,' emax = ', f12.6,' kcal/mol')
  417 format(1x,'             x               v-inp        v-fit   ',
     *          '     diff(u.a.)  diff(Kcal)',/)
      molecula=2
      write(6,1015) nombre
      read(5,*) n,mt,vex12(1),vex12(2),e012
      print*,'input parameters ', n,mt,vex12(1),vex12(2),e012
      if (n.gt.nmax) stop 'number of points(n) greater than nmax'
      if (mt.gt.n) stop 'polynomial order greater than number of
     , points'
      do 1 i=1,n
         read(5,*) r12(i),xener(i)
         xener(i)=xener(i)-e012
1     continue
      mt12=mt
      identi=12
      iter = 0
      if (ifor.eq.1.and.numiter.le.0) write(6,1010)
      nfas=mt
      print*,'number of SAF= ',nfas
      if (numiter.eq.0) call minim(cal2,vex12,rmsf,2,lim,aux2,ier)
      vex12(1)=dabs(vex12(1))
      vex12(2)=dabs(vex12(2))
      if (ifor.eq.1) write(6,417)
      ifor1=ifor
      if (ifor.eq.1.or.numiter.gt.0) bux=cal2(vex12)
      ifor1=0
      write(6,1020) rms,rms*bfact,emax*bfact
      write(6,*) vex12(1),vex12(2)
      call digual(c12,nexp2,sm,nexp,mt12)
      call scr2(nombre,0.d0,n,mt12,vex12,rms,emax)
      return
      end
*********************************************************************
*************************************************************************
      subroutine diat(mg,expo,r,c,e0,e)
*************************************************************************
      implicit real*8 (a-h,o-z)
      include 'dimensions.inc'
      dimension c(nexp2),expo(2)
      bux=r*dexp(-expo(1)*r)
      aux=1.d0
      e=e0+c(1)*dexp(-expo(2)*r)/r
      do 1 i=2,mg
         aux=aux*bux
         e=e+c(i)*aux
1     continue
      return
      end
*************************************************************************
      subroutine gfab(n,mt,vec2)
************************************************************************
c     SAF generation for a diatomic molecule
c     returns the value of the "j" SAF for each "i" point in x(j,i)
************************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax)
      dimension vec2(2)
      do 1 i=1,n
         xx(1,i) = dexp(-vec2(2)*r12(i))/r12(i)
         bux= dexp(-vec2(1)*r12(i))*r12(i)
         aux=1.d0
         do 1 j=2,mt
            aux=aux*bux
            xx(j,i)=aux
   1  continue
      return
      end
************************************************************************
************************************************************************
      subroutine gfaaa(n,mt,vex1)
************************************************************************
c     returns the value of the "j" SAF for each "i" point in x(j,i)
************************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      do 1 i=1,n
         aux12 = r12(i)*dexp(-vex1*r12(i))
         aux13 = r13(i)*dexp(-vex1*r13(i))
         aux23 = r23(i)*dexp(-vex1*r23(i))
         f12(0)=1.d0
         f13(0)=1.d0
         f23(0)=1.d0
         do 2 k=1,mt-1
            f12(k)=f12(k-1)*aux12
            f13(k)=f13(k-1)*aux13
            f23(k)=f23(k-1)*aux23
   2     continue
         do 3 j=1,nfas31(mt)
            if(ipa31(j,4).eq.1) then
              xx(j,i)=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
            else
            if(ipa31(j,4).eq.3) then
              xux=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
              yux=f12(ipa31(j,2))*f13(ipa31(j,1))*f23(ipa31(j,3))
              zux=f12(ipa31(j,2))*f13(ipa31(j,3))*f23(ipa31(j,1))
              xx(j,i)=xux+yux+zux
            else
            if(ipa31(j,4).eq.6) then
              xux=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
              yux=f12(ipa31(j,1))*f13(ipa31(j,3))*f23(ipa31(j,2))
              zux=f12(ipa31(j,2))*f13(ipa31(j,1))*f23(ipa31(j,3))
              xuy=f12(ipa31(j,2))*f13(ipa31(j,3))*f23(ipa31(j,1))
              yuy=f12(ipa31(j,3))*f13(ipa31(j,1))*f23(ipa31(j,2))
              zuy=f12(ipa31(j,3))*f13(ipa31(j,2))*f23(ipa31(j,1))
              xx(j,i)=xux+yux+zux+xuy+yuy+zuy
            endif
            endif
            endif
   3     continue
   1  continue
      return
      end
****************************************************************
      subroutine triaaa(mt,r12,r13,r23,vex1,cf,ener)
c     computes the three-body term of the triatomic system class AAA
************************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      dimension cf(nexp3)
      aux12 = r12*dexp(-vex1*r12)
      aux13 = r13*dexp(-vex1*r13)
      aux23 = r23*dexp(-vex1*r23)
      f12(0)=1.d0
      f13(0)=1.d0
      f23(0)=1.d0
      do 1 k=1,mt-1
         f12(k)=f12(k-1)*aux12
         f13(k)=f13(k-1)*aux13
         f23(k)=f23(k-1)*aux23
  1   continue
      ener=0.d0
      do 2 j=1,nfas31(mt)
         if(ipa31(j,4).eq.1) then
           aux=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
         else
         if(ipa31(j,4).eq.3) then
           xux=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
           yux=f12(ipa31(j,2))*f13(ipa31(j,1))*f23(ipa31(j,3))
           zux=f12(ipa31(j,2))*f13(ipa31(j,3))*f23(ipa31(j,1))
           aux=xux+yux+zux
         else
         if(ipa31(j,4).eq.6) then
           xux=f12(ipa31(j,1))*f13(ipa31(j,2))*f23(ipa31(j,3))
           yux=f12(ipa31(j,1))*f13(ipa31(j,3))*f23(ipa31(j,2))
           zux=f12(ipa31(j,2))*f13(ipa31(j,1))*f23(ipa31(j,3))
           xuy=f12(ipa31(j,2))*f13(ipa31(j,3))*f23(ipa31(j,1))
           yuy=f12(ipa31(j,3))*f13(ipa31(j,1))*f23(ipa31(j,2))
           zuy=f12(ipa31(j,3))*f13(ipa31(j,2))*f23(ipa31(j,1))
           aux=xux+yux+zux+xuy+yuy+zuy
         endif
         endif
         endif
         ener=ener+cf(j)*aux
   2  continue
      return
      end
************************************************************************
      subroutine gfabb(n,mt,vex2)
c     returns the value of the "j" SAF for each "i" point in x(j,i)
************************************************************************
c     generacion de las fas para cada punto
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      dimension vex2(2)
      do 1 i=1,n
         aux12 = r12(i)*dexp(-vex2(1)*r12(i))
         aux13 = r13(i)*dexp(-vex2(1)*r13(i))
         aux23 = r23(i)*dexp(-vex2(2)*r23(i))
         f12(0)=1.d0
         f13(0)=1.d0
         f23(0)=1.d0
         do 2 k=1,mt-1
            f12(k)=f12(k-1)*aux12
            f13(k)=f13(k-1)*aux13
            f23(k)=f23(k-1)*aux23
   2     continue
         do 3 j=1,nfas32(mt)
            if(ipa32(j,4).eq.1) then
              xx(j,i)=f12(ipa32(j,1))*f13(ipa32(j,3))*f23(ipa32(j,2))
            elseif(ipa32(j,4).eq.2) then
              xux=f12(ipa32(j,1))*f13(ipa32(j,3))
              zux=f12(ipa32(j,3))*f13(ipa32(j,1))
              xx(j,i)=(xux+zux)*f23(ipa32(j,2))
            endif
   3     continue
   1  continue
      return
      end
************************************************************************
      subroutine triabb(mt,r12,r13,r23,vex2,cf,ener)
c     computes the three-body term of the triatomic system class ABB
************************************************************************
c     SAF generation for each point
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      dimension cf(nexp3),vex2(2)
      aux12 = r12*dexp(-vex2(1)*r12)
      aux13 = r13*dexp(-vex2(1)*r13)
      aux23 = r23*dexp(-vex2(2)*r23)
      f12(0)=1.d0
      f13(0)=1.d0
      f23(0)=1.d0
      do 1 k=1,mt-1
         f12(k)=f12(k-1)*aux12
         f13(k)=f13(k-1)*aux13
         f23(k)=f23(k-1)*aux23
   1  continue
      ener=0.d0
      do 2 j=1,nfas32(mt)
         if(ipa32(j,4).eq.1) then
           aux=f12(ipa32(j,1))*f13(ipa32(j,3))*f23(ipa32(j,2))
         elseif(ipa32(j,4).eq.2) then
           xux=f12(ipa32(j,1))*f13(ipa32(j,3))
           zux=f12(ipa32(j,3))*f13(ipa32(j,1))
           aux=(xux+zux)*f23(ipa32(j,2))
         endif
         ener=ener+cf(j)*aux
   2  continue
      return
      end
************************************************************************
      subroutine gfabc(n,mt,vex3)
c     returns the value of the "j" SAF for each "i" point in x(j,i)
************************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      dimension vex3(3)
      do 1 i=1,n
         aux12 = r12(i)*dexp(-vex3(1)*r12(i))
         aux13 = r13(i)*dexp(-vex3(2)*r13(i))
         aux23 = r23(i)*dexp(-vex3(3)*r23(i))
         f12(0)=1.d0
         f13(0)=1.d0
         f23(0)=1.d0
         do 2 k=1,mt-1
            f12(k)=f12(k-1)*aux12
            f13(k)=f13(k-1)*aux13
            f23(k)=f23(k-1)*aux23
   2     continue
         do 3 j=1,nfas33(mt)
            xx(j,i)=f12(ipa33(j,1))*f13(ipa33(j,2))*f23(ipa33(j,3))
   3     continue
   1  continue
      return
      end
************************************************************************
      subroutine triabc(mt,r12,r13,r23,vex3,cf,ener)
************************************************************************
c     computes the three-body term of the triatomic system class ABC
************************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      common/basu3/f12(0:igrad3-1),f13(0:igrad3-1),f23(0:igrad3-1)
      dimension cf(nexp3),vex3(3)
      aux12 = r12*dexp(-vex3(1)*r12)
      aux13 = r13*dexp(-vex3(2)*r13)
      aux23 = r23*dexp(-vex3(3)*r23)
      f12(0)=1.d0
      f13(0)=1.d0
      f23(0)=1.d0
      do 1 k=1,mt-1
         f12(k)=f12(k-1)*aux12
         f13(k)=f13(k-1)*aux13
         f23(k)=f23(k-1)*aux23
   1  continue
      ener=0.d0
      do 2 j=1,nfas33(mt)
         aux=f12(ipa33(j,1))*f13(ipa33(j,2))*f23(ipa33(j,3))
         ener=ener+cf(j)*aux
   2  continue
      return
      end
**********************************************************************
      subroutine gfabcd(n,mt,vex6)
**********************************************************************
c     returns the value of the "j" SAF for each "i" point in xx(j,i)
**********************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot45/ipa45(nexp45,6),nfas45(igrad4)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      dimension vex6(6)
      do 1 i=1,n
         f12(0) = 1.d0
         f13(0) = 1.d0
         f14(0) = 1.d0
         f23(0) = 1.d0
         f24(0) = 1.d0
         f34(0) = 1.d0
         aux12= r12(i)*dexp(-vex6(1)*r12(i))
         aux13= r13(i)*dexp(-vex6(2)*r13(i))
         aux14= r14(i)*dexp(-vex6(3)*r14(i))
         aux23= r23(i)*dexp(-vex6(4)*r23(i))
         aux24= r24(i)*dexp(-vex6(5)*r24(i))
         aux34= r34(i)*dexp(-vex6(6)*r34(i))
         do 2 k=1,mt-2
            f12(k) = aux12*f12(k-1)
            f13(k) = aux13*f13(k-1)
            f14(k) = aux14*f14(k-1)
            f23(k) = aux23*f23(k-1)
            f24(k) = aux24*f24(k-1)
            f34(k) = aux34*f34(k-1)
2        continue
         do 3 j=1,nfas45(mt)
            i1=ipa45(j,1)
            i2=ipa45(j,2)
            i3=ipa45(j,3)
            i4=ipa45(j,4)
            i5=ipa45(j,5)
            i6=ipa45(j,6)
            xx(j,i)=f12(i1)*f13(i2)*f14(i3)*f23(i4)*f24(i5)*f34(i6)
 3       continue
 1    continue
      return
      end
***********************************************************************
      blockdata AA
      parameter (igrad3=10)
      parameter (igrad4=12)
      parameter (nexp3=255)
      parameter (nexp41=865,nexp42=756,nexp43=657)
      parameter (nexp44=639,nexp45=582)
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      data (ipa31(   1,j),j=1,4)/ 0,1,1,3/
      data (ipa31(   2,j),j=1,4)/ 1,1,1,1/
      data (ipa31(   3,j),j=1,4)/ 2,1,0,6/
      data (ipa31(   4,j),j=1,4)/ 2,1,1,3/
      data (ipa31(   5,j),j=1,4)/ 0,2,2,3/
      data (ipa31(   6,j),j=1,4)/ 3,1,0,6/
      data (ipa31(   7,j),j=1,4)/ 1,2,2,3/
      data (ipa31(   8,j),j=1,4)/ 3,1,1,3/
      data (ipa31(   9,j),j=1,4)/ 3,2,0,6/
      data (ipa31(  10,j),j=1,4)/ 4,1,0,6/
      data (ipa31(  11,j),j=1,4)/ 2,2,2,1/
      data (ipa31(  12,j),j=1,4)/ 3,2,1,6/
      data (ipa31(  13,j),j=1,4)/ 0,3,3,3/
      data (ipa31(  14,j),j=1,4)/ 4,1,1,3/
      data (ipa31(  15,j),j=1,4)/ 4,2,0,6/
      data (ipa31(  16,j),j=1,4)/ 5,1,0,6/
      data (ipa31(  17,j),j=1,4)/ 3,2,2,3/
      data (ipa31(  18,j),j=1,4)/ 1,3,3,3/
      data (ipa31(  19,j),j=1,4)/ 4,2,1,6/
      data (ipa31(  20,j),j=1,4)/ 4,3,0,6/
      data (ipa31(  21,j),j=1,4)/ 5,1,1,3/
      data (ipa31(  22,j),j=1,4)/ 5,2,0,6/
      data (ipa31(  23,j),j=1,4)/ 6,1,0,6/
      data (ipa31(  24,j),j=1,4)/ 2,3,3,3/
      data (ipa31(  25,j),j=1,4)/ 4,2,2,3/
      data (ipa31(  26,j),j=1,4)/ 4,3,1,6/
      data (ipa31(  27,j),j=1,4)/ 0,4,4,3/
      data (ipa31(  28,j),j=1,4)/ 5,2,1,6/
      data (ipa31(  29,j),j=1,4)/ 5,3,0,6/
      data (ipa31(  30,j),j=1,4)/ 6,1,1,3/
      data (ipa31(  31,j),j=1,4)/ 6,2,0,6/
      data (ipa31(  32,j),j=1,4)/ 7,1,0,6/
      data (ipa31(  33,j),j=1,4)/ 3,3,3,1/
      data (ipa31(  34,j),j=1,4)/ 4,3,2,6/
      data (ipa31(  35,j),j=1,4)/ 1,4,4,3/
      data (ipa31(  36,j),j=1,4)/ 5,2,2,3/
      data (ipa31(  37,j),j=1,4)/ 5,3,1,6/
      data (ipa31(  38,j),j=1,4)/ 5,4,0,6/
      data (ipa31(  39,j),j=1,4)/ 6,2,1,6/
      data (ipa31(  40,j),j=1,4)/ 6,3,0,6/
      data (ipa31(  41,j),j=1,4)/ 7,1,1,3/
      data (ipa31(  42,j),j=1,4)/ 7,2,0,6/
      data (ipa31(  43,j),j=1,4)/ 8,1,0,6/
      data (ipa31(  44,j),j=1,4)/ 4,3,3,3/
      data (ipa31(  45,j),j=1,4)/ 2,4,4,3/
      data (ipa31(  46,j),j=1,4)/ 5,3,2,6/
      data (ipa31(  47,j),j=1,4)/ 5,4,1,6/
      data (ipa31(  48,j),j=1,4)/ 0,5,5,3/
      data (ipa31(  49,j),j=1,4)/ 6,2,2,3/
      data (ipa31(  50,j),j=1,4)/ 6,3,1,6/
      data (ipa31(  51,j),j=1,4)/ 6,4,0,6/
      data (ipa31(  52,j),j=1,4)/ 7,2,1,6/
      data (ipa31(  53,j),j=1,4)/ 7,3,0,6/
      data (ipa31(  54,j),j=1,4)/ 8,1,1,3/
      data (ipa31(  55,j),j=1,4)/ 8,2,0,6/
      data (ipa31(  56,j),j=1,4)/ 9,1,0,6/
      data (nfas31(i),i=1,10)/0,1,3,6,10,16,23,32,43,56/
      data (ipa32(   1,j),j=1,4)/ 1,1,0,2/
      data (ipa32(   2,j),j=1,4)/ 1,0,1,1/
      data (ipa32(   3,j),j=1,4)/ 1,1,1,1/
      data (ipa32(   4,j),j=1,4)/ 2,1,0,2/
      data (ipa32(   5,j),j=1,4)/ 2,0,1,2/
      data (ipa32(   6,j),j=1,4)/ 0,2,1,2/
      data (ipa32(   7,j),j=1,4)/ 2,1,1,2/
      data (ipa32(   8,j),j=1,4)/ 1,2,1,1/
      data (ipa32(   9,j),j=1,4)/ 2,2,0,2/
      data (ipa32(  10,j),j=1,4)/ 2,0,2,1/
      data (ipa32(  11,j),j=1,4)/ 3,1,0,2/
      data (ipa32(  12,j),j=1,4)/ 3,0,1,2/
      data (ipa32(  13,j),j=1,4)/ 0,3,1,2/
      data (ipa32(  14,j),j=1,4)/ 2,2,1,2/
      data (ipa32(  15,j),j=1,4)/ 2,1,2,1/
      data (ipa32(  16,j),j=1,4)/ 3,1,1,2/
      data (ipa32(  17,j),j=1,4)/ 1,3,1,1/
      data (ipa32(  18,j),j=1,4)/ 3,2,0,2/
      data (ipa32(  19,j),j=1,4)/ 3,0,2,2/
      data (ipa32(  20,j),j=1,4)/ 0,3,2,2/
      data (ipa32(  21,j),j=1,4)/ 4,1,0,2/
      data (ipa32(  22,j),j=1,4)/ 4,0,1,2/
      data (ipa32(  23,j),j=1,4)/ 0,4,1,2/
      data (ipa32(  24,j),j=1,4)/ 2,2,2,1/
      data (ipa32(  25,j),j=1,4)/ 3,2,1,2/
      data (ipa32(  26,j),j=1,4)/ 3,1,2,2/
      data (ipa32(  27,j),j=1,4)/ 1,3,2,2/
      data (ipa32(  28,j),j=1,4)/ 3,3,0,2/
      data (ipa32(  29,j),j=1,4)/ 3,0,3,1/
      data (ipa32(  30,j),j=1,4)/ 4,1,1,2/
      data (ipa32(  31,j),j=1,4)/ 1,4,1,1/
      data (ipa32(  32,j),j=1,4)/ 4,2,0,2/
      data (ipa32(  33,j),j=1,4)/ 4,0,2,2/
      data (ipa32(  34,j),j=1,4)/ 0,4,2,2/
      data (ipa32(  35,j),j=1,4)/ 5,1,0,2/
      data (ipa32(  36,j),j=1,4)/ 5,0,1,2/
      data (ipa32(  37,j),j=1,4)/ 0,5,1,2/
      data (ipa32(  38,j),j=1,4)/ 3,2,2,2/
      data (ipa32(  39,j),j=1,4)/ 2,3,2,1/
      data (ipa32(  40,j),j=1,4)/ 3,3,1,2/
      data (ipa32(  41,j),j=1,4)/ 3,1,3,1/
      data (ipa32(  42,j),j=1,4)/ 4,2,1,2/
      data (ipa32(  43,j),j=1,4)/ 4,1,2,2/
      data (ipa32(  44,j),j=1,4)/ 1,4,2,2/
      data (ipa32(  45,j),j=1,4)/ 4,3,0,2/
      data (ipa32(  46,j),j=1,4)/ 4,0,3,2/
      data (ipa32(  47,j),j=1,4)/ 0,4,3,2/
      data (ipa32(  48,j),j=1,4)/ 5,1,1,2/
      data (ipa32(  49,j),j=1,4)/ 1,5,1,1/
      data (ipa32(  50,j),j=1,4)/ 5,2,0,2/
      data (ipa32(  51,j),j=1,4)/ 5,0,2,2/
      data (ipa32(  52,j),j=1,4)/ 0,5,2,2/
      data (ipa32(  53,j),j=1,4)/ 6,1,0,2/
      data (ipa32(  54,j),j=1,4)/ 6,0,1,2/
      data (ipa32(  55,j),j=1,4)/ 0,6,1,2/
      data (ipa32(  56,j),j=1,4)/ 3,3,2,2/
      data (ipa32(  57,j),j=1,4)/ 3,2,3,1/
      data (ipa32(  58,j),j=1,4)/ 4,2,2,2/
      data (ipa32(  59,j),j=1,4)/ 2,4,2,1/
      data (ipa32(  60,j),j=1,4)/ 4,3,1,2/
      data (ipa32(  61,j),j=1,4)/ 4,1,3,2/
      data (ipa32(  62,j),j=1,4)/ 1,4,3,2/
      data (ipa32(  63,j),j=1,4)/ 4,4,0,2/
      data (ipa32(  64,j),j=1,4)/ 4,0,4,1/
      data (ipa32(  65,j),j=1,4)/ 5,2,1,2/
      data (ipa32(  66,j),j=1,4)/ 5,1,2,2/
      data (ipa32(  67,j),j=1,4)/ 1,5,2,2/
      data (ipa32(  68,j),j=1,4)/ 5,3,0,2/
      data (ipa32(  69,j),j=1,4)/ 5,0,3,2/
      data (ipa32(  70,j),j=1,4)/ 0,5,3,2/
      data (ipa32(  71,j),j=1,4)/ 6,1,1,2/
      data (ipa32(  72,j),j=1,4)/ 1,6,1,1/
      data (ipa32(  73,j),j=1,4)/ 6,2,0,2/
      data (ipa32(  74,j),j=1,4)/ 6,0,2,2/
      data (ipa32(  75,j),j=1,4)/ 0,6,2,2/
      data (ipa32(  76,j),j=1,4)/ 7,1,0,2/
      data (ipa32(  77,j),j=1,4)/ 7,0,1,2/
      data (ipa32(  78,j),j=1,4)/ 0,7,1,2/
      data (ipa32(  79,j),j=1,4)/ 3,3,3,1/
      data (ipa32(  80,j),j=1,4)/ 4,3,2,2/
      data (ipa32(  81,j),j=1,4)/ 4,2,3,2/
      data (ipa32(  82,j),j=1,4)/ 2,4,3,2/
      data (ipa32(  83,j),j=1,4)/ 4,4,1,2/
      data (ipa32(  84,j),j=1,4)/ 4,1,4,1/
      data (ipa32(  85,j),j=1,4)/ 5,2,2,2/
      data (ipa32(  86,j),j=1,4)/ 2,5,2,1/
      data (ipa32(  87,j),j=1,4)/ 5,3,1,2/
      data (ipa32(  88,j),j=1,4)/ 5,1,3,2/
      data (ipa32(  89,j),j=1,4)/ 1,5,3,2/
      data (ipa32(  90,j),j=1,4)/ 5,4,0,2/
      data (ipa32(  91,j),j=1,4)/ 5,0,4,2/
      data (ipa32(  92,j),j=1,4)/ 0,5,4,2/
      data (ipa32(  93,j),j=1,4)/ 6,2,1,2/
      data (ipa32(  94,j),j=1,4)/ 6,1,2,2/
      data (ipa32(  95,j),j=1,4)/ 1,6,2,2/
      data (ipa32(  96,j),j=1,4)/ 6,3,0,2/
      data (ipa32(  97,j),j=1,4)/ 6,0,3,2/
      data (ipa32(  98,j),j=1,4)/ 0,6,3,2/
      data (ipa32(  99,j),j=1,4)/ 7,1,1,2/
      data (ipa32( 100,j),j=1,4)/ 1,7,1,1/
      data (ipa32( 101,j),j=1,4)/ 7,2,0,2/
      data (ipa32( 102,j),j=1,4)/ 7,0,2,2/
      data (ipa32( 103,j),j=1,4)/ 0,7,2,2/
      data (ipa32( 104,j),j=1,4)/ 8,1,0,2/
      data (ipa32( 105,j),j=1,4)/ 8,0,1,2/
      data (ipa32( 106,j),j=1,4)/ 0,8,1,2/
      data (ipa32( 107,j),j=1,4)/ 4,3,3,2/
      data (ipa32( 108,j),j=1,4)/ 3,4,3,1/
      data (ipa32( 109,j),j=1,4)/ 4,4,2,2/
      data (ipa32( 110,j),j=1,4)/ 4,2,4,1/
      data (ipa32( 111,j),j=1,4)/ 5,3,2,2/
      data (ipa32( 112,j),j=1,4)/ 5,2,3,2/
      data (ipa32( 113,j),j=1,4)/ 2,5,3,2/
      data (ipa32( 114,j),j=1,4)/ 5,4,1,2/
      data (ipa32( 115,j),j=1,4)/ 5,1,4,2/
      data (ipa32( 116,j),j=1,4)/ 1,5,4,2/
      data (ipa32( 117,j),j=1,4)/ 5,5,0,2/
      data (ipa32( 118,j),j=1,4)/ 5,0,5,1/
      data (ipa32( 119,j),j=1,4)/ 6,2,2,2/
      data (ipa32( 120,j),j=1,4)/ 2,6,2,1/
      data (ipa32( 121,j),j=1,4)/ 6,3,1,2/
      data (ipa32( 122,j),j=1,4)/ 6,1,3,2/
      data (ipa32( 123,j),j=1,4)/ 1,6,3,2/
      data (ipa32( 124,j),j=1,4)/ 6,4,0,2/
      data (ipa32( 125,j),j=1,4)/ 6,0,4,2/
      data (ipa32( 126,j),j=1,4)/ 0,6,4,2/
      data (ipa32( 127,j),j=1,4)/ 7,2,1,2/
      data (ipa32( 128,j),j=1,4)/ 7,1,2,2/
      data (ipa32( 129,j),j=1,4)/ 1,7,2,2/
      data (ipa32( 130,j),j=1,4)/ 7,3,0,2/
      data (ipa32( 131,j),j=1,4)/ 7,0,3,2/
      data (ipa32( 132,j),j=1,4)/ 0,7,3,2/
      data (ipa32( 133,j),j=1,4)/ 8,1,1,2/
      data (ipa32( 134,j),j=1,4)/ 1,8,1,1/
      data (ipa32( 135,j),j=1,4)/ 8,2,0,2/
      data (ipa32( 136,j),j=1,4)/ 8,0,2,2/
      data (ipa32( 137,j),j=1,4)/ 0,8,2,2/
      data (ipa32( 138,j),j=1,4)/ 9,1,0,2/
      data (ipa32( 139,j),j=1,4)/ 9,0,1,2/
      data (ipa32( 140,j),j=1,4)/ 0,9,1,2/
      data (nfas32(i),i=1,10)/0,2,6,13,23,37,55,78,106,140/
      data (ipa33(   1,j),j=1,3)/ 0,1,1/
      data (ipa33(   2,j),j=1,3)/ 1,0,1/
      data (ipa33(   3,j),j=1,3)/ 1,1,0/
      data (ipa33(   4,j),j=1,3)/ 0,1,2/
      data (ipa33(   5,j),j=1,3)/ 0,2,1/
      data (ipa33(   6,j),j=1,3)/ 1,0,2/
      data (ipa33(   7,j),j=1,3)/ 1,1,1/
      data (ipa33(   8,j),j=1,3)/ 1,2,0/
      data (ipa33(   9,j),j=1,3)/ 2,0,1/
      data (ipa33(  10,j),j=1,3)/ 2,1,0/
      data (ipa33(  11,j),j=1,3)/ 0,1,3/
      data (ipa33(  12,j),j=1,3)/ 0,2,2/
      data (ipa33(  13,j),j=1,3)/ 0,3,1/
      data (ipa33(  14,j),j=1,3)/ 1,0,3/
      data (ipa33(  15,j),j=1,3)/ 1,1,2/
      data (ipa33(  16,j),j=1,3)/ 1,2,1/
      data (ipa33(  17,j),j=1,3)/ 1,3,0/
      data (ipa33(  18,j),j=1,3)/ 2,0,2/
      data (ipa33(  19,j),j=1,3)/ 2,1,1/
      data (ipa33(  20,j),j=1,3)/ 2,2,0/
      data (ipa33(  21,j),j=1,3)/ 3,0,1/
      data (ipa33(  22,j),j=1,3)/ 3,1,0/
      data (ipa33(  23,j),j=1,3)/ 0,1,4/
      data (ipa33(  24,j),j=1,3)/ 0,2,3/
      data (ipa33(  25,j),j=1,3)/ 0,3,2/
      data (ipa33(  26,j),j=1,3)/ 0,4,1/
      data (ipa33(  27,j),j=1,3)/ 1,0,4/
      data (ipa33(  28,j),j=1,3)/ 1,1,3/
      data (ipa33(  29,j),j=1,3)/ 1,2,2/
      data (ipa33(  30,j),j=1,3)/ 1,3,1/
      data (ipa33(  31,j),j=1,3)/ 1,4,0/
      data (ipa33(  32,j),j=1,3)/ 2,0,3/
      data (ipa33(  33,j),j=1,3)/ 2,1,2/
      data (ipa33(  34,j),j=1,3)/ 2,2,1/
      data (ipa33(  35,j),j=1,3)/ 2,3,0/
      data (ipa33(  36,j),j=1,3)/ 3,0,2/
      data (ipa33(  37,j),j=1,3)/ 3,1,1/
      data (ipa33(  38,j),j=1,3)/ 3,2,0/
      data (ipa33(  39,j),j=1,3)/ 4,0,1/
      data (ipa33(  40,j),j=1,3)/ 4,1,0/
      data (ipa33(  41,j),j=1,3)/ 0,1,5/
      data (ipa33(  42,j),j=1,3)/ 0,2,4/
      data (ipa33(  43,j),j=1,3)/ 0,3,3/
      data (ipa33(  44,j),j=1,3)/ 0,4,2/
      data (ipa33(  45,j),j=1,3)/ 0,5,1/
      data (ipa33(  46,j),j=1,3)/ 1,0,5/
      data (ipa33(  47,j),j=1,3)/ 1,1,4/
      data (ipa33(  48,j),j=1,3)/ 1,2,3/
      data (ipa33(  49,j),j=1,3)/ 1,3,2/
      data (ipa33(  50,j),j=1,3)/ 1,4,1/
      data (ipa33(  51,j),j=1,3)/ 1,5,0/
      data (ipa33(  52,j),j=1,3)/ 2,0,4/
      data (ipa33(  53,j),j=1,3)/ 2,1,3/
      data (ipa33(  54,j),j=1,3)/ 2,2,2/
      data (ipa33(  55,j),j=1,3)/ 2,3,1/
      data (ipa33(  56,j),j=1,3)/ 2,4,0/
      data (ipa33(  57,j),j=1,3)/ 3,0,3/
      data (ipa33(  58,j),j=1,3)/ 3,1,2/
      data (ipa33(  59,j),j=1,3)/ 3,2,1/
      data (ipa33(  60,j),j=1,3)/ 3,3,0/
      data (ipa33(  61,j),j=1,3)/ 4,0,2/
      data (ipa33(  62,j),j=1,3)/ 4,1,1/
      data (ipa33(  63,j),j=1,3)/ 4,2,0/
      data (ipa33(  64,j),j=1,3)/ 5,0,1/
      data (ipa33(  65,j),j=1,3)/ 5,1,0/
      data (ipa33(  66,j),j=1,3)/ 0,1,6/
      data (ipa33(  67,j),j=1,3)/ 0,2,5/
      data (ipa33(  68,j),j=1,3)/ 0,3,4/
      data (ipa33(  69,j),j=1,3)/ 0,4,3/
      data (ipa33(  70,j),j=1,3)/ 0,5,2/
      data (ipa33(  71,j),j=1,3)/ 0,6,1/
      data (ipa33(  72,j),j=1,3)/ 1,0,6/
      data (ipa33(  73,j),j=1,3)/ 1,1,5/
      data (ipa33(  74,j),j=1,3)/ 1,2,4/
      data (ipa33(  75,j),j=1,3)/ 1,3,3/
      data (ipa33(  76,j),j=1,3)/ 1,4,2/
      data (ipa33(  77,j),j=1,3)/ 1,5,1/
      data (ipa33(  78,j),j=1,3)/ 1,6,0/
      data (ipa33(  79,j),j=1,3)/ 2,0,5/
      data (ipa33(  80,j),j=1,3)/ 2,1,4/
      data (ipa33(  81,j),j=1,3)/ 2,2,3/
      data (ipa33(  82,j),j=1,3)/ 2,3,2/
      data (ipa33(  83,j),j=1,3)/ 2,4,1/
      data (ipa33(  84,j),j=1,3)/ 2,5,0/
      data (ipa33(  85,j),j=1,3)/ 3,0,4/
      data (ipa33(  86,j),j=1,3)/ 3,1,3/
      data (ipa33(  87,j),j=1,3)/ 3,2,2/
      data (ipa33(  88,j),j=1,3)/ 3,3,1/
      data (ipa33(  89,j),j=1,3)/ 3,4,0/
      data (ipa33(  90,j),j=1,3)/ 4,0,3/
      data (ipa33(  91,j),j=1,3)/ 4,1,2/
      data (ipa33(  92,j),j=1,3)/ 4,2,1/
      data (ipa33(  93,j),j=1,3)/ 4,3,0/
      data (ipa33(  94,j),j=1,3)/ 5,0,2/
      data (ipa33(  95,j),j=1,3)/ 5,1,1/
      data (ipa33(  96,j),j=1,3)/ 5,2,0/
      data (ipa33(  97,j),j=1,3)/ 6,0,1/
      data (ipa33(  98,j),j=1,3)/ 6,1,0/
      data (ipa33(  99,j),j=1,3)/ 0,1,7/
      data (ipa33( 100,j),j=1,3)/ 0,2,6/
      data (ipa33( 101,j),j=1,3)/ 0,3,5/
      data (ipa33( 102,j),j=1,3)/ 0,4,4/
      data (ipa33( 103,j),j=1,3)/ 0,5,3/
      data (ipa33( 104,j),j=1,3)/ 0,6,2/
      data (ipa33( 105,j),j=1,3)/ 0,7,1/
      data (ipa33( 106,j),j=1,3)/ 1,0,7/
      data (ipa33( 107,j),j=1,3)/ 1,1,6/
      data (ipa33( 108,j),j=1,3)/ 1,2,5/
      data (ipa33( 109,j),j=1,3)/ 1,3,4/
      data (ipa33( 110,j),j=1,3)/ 1,4,3/
      data (ipa33( 111,j),j=1,3)/ 1,5,2/
      data (ipa33( 112,j),j=1,3)/ 1,6,1/
      data (ipa33( 113,j),j=1,3)/ 1,7,0/
      data (ipa33( 114,j),j=1,3)/ 2,0,6/
      data (ipa33( 115,j),j=1,3)/ 2,1,5/
      data (ipa33( 116,j),j=1,3)/ 2,2,4/
      data (ipa33( 117,j),j=1,3)/ 2,3,3/
      data (ipa33( 118,j),j=1,3)/ 2,4,2/
      data (ipa33( 119,j),j=1,3)/ 2,5,1/
      data (ipa33( 120,j),j=1,3)/ 2,6,0/
      data (ipa33( 121,j),j=1,3)/ 3,0,5/
      data (ipa33( 122,j),j=1,3)/ 3,1,4/
      data (ipa33( 123,j),j=1,3)/ 3,2,3/
      data (ipa33( 124,j),j=1,3)/ 3,3,2/
      data (ipa33( 125,j),j=1,3)/ 3,4,1/
      data (ipa33( 126,j),j=1,3)/ 3,5,0/
      data (ipa33( 127,j),j=1,3)/ 4,0,4/
      data (ipa33( 128,j),j=1,3)/ 4,1,3/
      data (ipa33( 129,j),j=1,3)/ 4,2,2/
      data (ipa33( 130,j),j=1,3)/ 4,3,1/
      data (ipa33( 131,j),j=1,3)/ 4,4,0/
      data (ipa33( 132,j),j=1,3)/ 5,0,3/
      data (ipa33( 133,j),j=1,3)/ 5,1,2/
      data (ipa33( 134,j),j=1,3)/ 5,2,1/
      data (ipa33( 135,j),j=1,3)/ 5,3,0/
      data (ipa33( 136,j),j=1,3)/ 6,0,2/
      data (ipa33( 137,j),j=1,3)/ 6,1,1/
      data (ipa33( 138,j),j=1,3)/ 6,2,0/
      data (ipa33( 139,j),j=1,3)/ 7,0,1/
      data (ipa33( 140,j),j=1,3)/ 7,1,0/
      data (ipa33( 141,j),j=1,3)/ 0,1,8/
      data (ipa33( 142,j),j=1,3)/ 0,2,7/
      data (ipa33( 143,j),j=1,3)/ 0,3,6/
      data (ipa33( 144,j),j=1,3)/ 0,4,5/
      data (ipa33( 145,j),j=1,3)/ 0,5,4/
      data (ipa33( 146,j),j=1,3)/ 0,6,3/
      data (ipa33( 147,j),j=1,3)/ 0,7,2/
      data (ipa33( 148,j),j=1,3)/ 0,8,1/
      data (ipa33( 149,j),j=1,3)/ 1,0,8/
      data (ipa33( 150,j),j=1,3)/ 1,1,7/
      data (ipa33( 151,j),j=1,3)/ 1,2,6/
      data (ipa33( 152,j),j=1,3)/ 1,3,5/
      data (ipa33( 153,j),j=1,3)/ 1,4,4/
      data (ipa33( 154,j),j=1,3)/ 1,5,3/
      data (ipa33( 155,j),j=1,3)/ 1,6,2/
      data (ipa33( 156,j),j=1,3)/ 1,7,1/
      data (ipa33( 157,j),j=1,3)/ 1,8,0/
      data (ipa33( 158,j),j=1,3)/ 2,0,7/
      data (ipa33( 159,j),j=1,3)/ 2,1,6/
      data (ipa33( 160,j),j=1,3)/ 2,2,5/
      data (ipa33( 161,j),j=1,3)/ 2,3,4/
      data (ipa33( 162,j),j=1,3)/ 2,4,3/
      data (ipa33( 163,j),j=1,3)/ 2,5,2/
      data (ipa33( 164,j),j=1,3)/ 2,6,1/
      data (ipa33( 165,j),j=1,3)/ 2,7,0/
      data (ipa33( 166,j),j=1,3)/ 3,0,6/
      data (ipa33( 167,j),j=1,3)/ 3,1,5/
      data (ipa33( 168,j),j=1,3)/ 3,2,4/
      data (ipa33( 169,j),j=1,3)/ 3,3,3/
      data (ipa33( 170,j),j=1,3)/ 3,4,2/
      data (ipa33( 171,j),j=1,3)/ 3,5,1/
      data (ipa33( 172,j),j=1,3)/ 3,6,0/
      data (ipa33( 173,j),j=1,3)/ 4,0,5/
      data (ipa33( 174,j),j=1,3)/ 4,1,4/
      data (ipa33( 175,j),j=1,3)/ 4,2,3/
      data (ipa33( 176,j),j=1,3)/ 4,3,2/
      data (ipa33( 177,j),j=1,3)/ 4,4,1/
      data (ipa33( 178,j),j=1,3)/ 4,5,0/
      data (ipa33( 179,j),j=1,3)/ 5,0,4/
      data (ipa33( 180,j),j=1,3)/ 5,1,3/
      data (ipa33( 181,j),j=1,3)/ 5,2,2/
      data (ipa33( 182,j),j=1,3)/ 5,3,1/
      data (ipa33( 183,j),j=1,3)/ 5,4,0/
      data (ipa33( 184,j),j=1,3)/ 6,0,3/
      data (ipa33( 185,j),j=1,3)/ 6,1,2/
      data (ipa33( 186,j),j=1,3)/ 6,2,1/
      data (ipa33( 187,j),j=1,3)/ 6,3,0/
      data (ipa33( 188,j),j=1,3)/ 7,0,2/
      data (ipa33( 189,j),j=1,3)/ 7,1,1/
      data (ipa33( 190,j),j=1,3)/ 7,2,0/
      data (ipa33( 191,j),j=1,3)/ 8,0,1/
      data (ipa33( 192,j),j=1,3)/ 8,1,0/
      data (ipa33( 193,j),j=1,3)/ 0,1,9/
      data (ipa33( 194,j),j=1,3)/ 0,2,8/
      data (ipa33( 195,j),j=1,3)/ 0,3,7/
      data (ipa33( 196,j),j=1,3)/ 0,4,6/
      data (ipa33( 197,j),j=1,3)/ 0,5,5/
      data (ipa33( 198,j),j=1,3)/ 0,6,4/
      data (ipa33( 199,j),j=1,3)/ 0,7,3/
      data (ipa33( 200,j),j=1,3)/ 0,8,2/
      data (ipa33( 201,j),j=1,3)/ 0,9,1/
      data (ipa33( 202,j),j=1,3)/ 1,0,9/
      data (ipa33( 203,j),j=1,3)/ 1,1,8/
      data (ipa33( 204,j),j=1,3)/ 1,2,7/
      data (ipa33( 205,j),j=1,3)/ 1,3,6/
      data (ipa33( 206,j),j=1,3)/ 1,4,5/
      data (ipa33( 207,j),j=1,3)/ 1,5,4/
      data (ipa33( 208,j),j=1,3)/ 1,6,3/
      data (ipa33( 209,j),j=1,3)/ 1,7,2/
      data (ipa33( 210,j),j=1,3)/ 1,8,1/
      data (ipa33( 211,j),j=1,3)/ 1,9,0/
      data (ipa33( 212,j),j=1,3)/ 2,0,8/
      data (ipa33( 213,j),j=1,3)/ 2,1,7/
      data (ipa33( 214,j),j=1,3)/ 2,2,6/
      data (ipa33( 215,j),j=1,3)/ 2,3,5/
      data (ipa33( 216,j),j=1,3)/ 2,4,4/
      data (ipa33( 217,j),j=1,3)/ 2,5,3/
      data (ipa33( 218,j),j=1,3)/ 2,6,2/
      data (ipa33( 219,j),j=1,3)/ 2,7,1/
      data (ipa33( 220,j),j=1,3)/ 2,8,0/
      data (ipa33( 221,j),j=1,3)/ 3,0,7/
      data (ipa33( 222,j),j=1,3)/ 3,1,6/
      data (ipa33( 223,j),j=1,3)/ 3,2,5/
      data (ipa33( 224,j),j=1,3)/ 3,3,4/
      data (ipa33( 225,j),j=1,3)/ 3,4,3/
      data (ipa33( 226,j),j=1,3)/ 3,5,2/
      data (ipa33( 227,j),j=1,3)/ 3,6,1/
      data (ipa33( 228,j),j=1,3)/ 3,7,0/
      data (ipa33( 229,j),j=1,3)/ 4,0,6/
      data (ipa33( 230,j),j=1,3)/ 4,1,5/
      data (ipa33( 231,j),j=1,3)/ 4,2,4/
      data (ipa33( 232,j),j=1,3)/ 4,3,3/
      data (ipa33( 233,j),j=1,3)/ 4,4,2/
      data (ipa33( 234,j),j=1,3)/ 4,5,1/
      data (ipa33( 235,j),j=1,3)/ 4,6,0/
      data (ipa33( 236,j),j=1,3)/ 5,0,5/
      data (ipa33( 237,j),j=1,3)/ 5,1,4/
      data (ipa33( 238,j),j=1,3)/ 5,2,3/
      data (ipa33( 239,j),j=1,3)/ 5,3,2/
      data (ipa33( 240,j),j=1,3)/ 5,4,1/
      data (ipa33( 241,j),j=1,3)/ 5,5,0/
      data (ipa33( 242,j),j=1,3)/ 6,0,4/
      data (ipa33( 243,j),j=1,3)/ 6,1,3/
      data (ipa33( 244,j),j=1,3)/ 6,2,2/
      data (ipa33( 245,j),j=1,3)/ 6,3,1/
      data (ipa33( 246,j),j=1,3)/ 6,4,0/
      data (ipa33( 247,j),j=1,3)/ 7,0,3/
      data (ipa33( 248,j),j=1,3)/ 7,1,2/
      data (ipa33( 249,j),j=1,3)/ 7,2,1/
      data (ipa33( 250,j),j=1,3)/ 7,3,0/
      data (ipa33( 251,j),j=1,3)/ 8,0,2/
      data (ipa33( 252,j),j=1,3)/ 8,1,1/
      data (ipa33( 253,j),j=1,3)/ 8,2,0/
      data (ipa33( 254,j),j=1,3)/ 9,0,1/
      data (ipa33( 255,j),j=1,3)/ 9,1,0/
      data (nfas33(i),i=1,10)/0,3,10,22,40,65,98,140,192,255/
      end
**********************************************************************
      subroutine gfaaaa(n,mt,vex1)
**********************************************************************
*     devuelve el valor de la fas j para cada punto i en xx(j,i)
**********************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/ipot41/ipa41(nexp41,7),nfas41(igrad4)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      do 1 i=1,n
         f12(0) = 1.d0
         f13(0) = 1.d0
         f14(0) = 1.d0
         f23(0) = 1.d0
         f24(0) = 1.d0
         f34(0) = 1.d0
         aux12= r12(i)*dexp(-vex1*r12(i))
         aux13= r13(i)*dexp(-vex1*r13(i))
         aux14= r14(i)*dexp(-vex1*r14(i))
         aux23= r23(i)*dexp(-vex1*r23(i))
         aux24= r24(i)*dexp(-vex1*r24(i))
         aux34= r34(i)*dexp(-vex1*r34(i))
*     como el grado es mt, la potencia de orden mayor es mt-3
         do 2 k=1,mt-2
            f12(k) = aux12*f12(k-1)
            f13(k) = aux13*f13(k-1)
            f14(k) = aux14*f14(k-1)
            f23(k) = aux23*f23(k-1)
            f24(k) = aux24*f24(k-1)
            f34(k) = aux34*f34(k-1)
2        continue
         do 3 j=1,nfas41(mt)
            i1=ipa41(j,1)
            i2=ipa41(j,2)
            i3=ipa41(j,3)
            i4=ipa41(j,4)
            i5=ipa41(j,5)
            i6=ipa41(j,6)
            i7=ipa41(j,7)
            xx(j,i)=faaaa(i1,i2,i3,i4,i5,i6,i7)
 3       continue
 1    continue
      return
      end
**********************************************************************
      function faaaa(i1,i2,i3,i4,i5,i6,i7)
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
      if (i7.eq.1) then
         faaaa=g(i1,i1,i1,i1,i1,i1)
         return
      else
      if (i7.eq.3) then
         faaaa=g(i1,i1,i3,i3,i1,i1)+g(i1,i3,i1,i1,i3,i1)+
     ,         g(i3,i1,i1,i1,i1,i3)
         return
      else
      if (i7.eq.4) then
         faaaa=g(i1,i1,i1,i4,i4,i4)+g(i1,i4,i4,i1,i1,i4)+
     ,         g(i4,i1,i4,i1,i4,i1)+g(i4,i4,i1,i4,i1,i1)
         return
      else
      if (i7.eq.6) then
         faaaa=g(i1,i2,i3,i4,i5,i6)+g(i1,i3,i2,i5,i4,i6)+
     ,         g(i3,i1,i2,i5,i6,i4)+g(i4,i5,i6,i1,i2,i3)+
     ,         g(i5,i4,i6,i1,i3,i2)+g(i5,i6,i4,i3,i1,i2)
         return
      else
      if (i7.eq.11) then
         faaaa=g(i1,i2,i3,i4,i2,i1)+g(i1,i2,i4,i3,i2,i1)+
     ,         g(i1,i3,i2,i2,i4,i1)+g(i1,i4,i2,i2,i3,i1)+
     ,         g(i2,i1,i3,i4,i1,i2)+g(i2,i1,i4,i3,i1,i2)+
     ,         g(i2,i3,i1,i1,i4,i2)+g(i2,i4,i1,i1,i3,i2)+
     ,         g(i3,i1,i2,i2,i1,i4)+g(i3,i2,i1,i1,i2,i4)+
     ,         g(i4,i1,i2,i2,i1,i3)+g(i4,i2,i1,i1,i2,i3)
         return
      else
      if (i7.eq.13) then
         faaaa=g(i1,i1,i3,i4,i6,i6)+g(i1,i6,i4,i3,i1,i6)+
     ,         g(i1,i3,i1,i6,i4,i6)+g(i1,i4,i6,i1,i3,i6)+
     ,         g(i6,i1,i4,i3,i6,i1)+g(i6,i6,i3,i4,i1,i1)+
     ,         g(i6,i3,i6,i1,i4,i1)+g(i6,i4,i1,i6,i3,i1)+
     ,         g(i3,i1,i1,i6,i6,i4)+g(i3,i6,i6,i1,i1,i4)+
     ,         g(i4,i1,i6,i1,i6,i3)+g(i4,i6,i1,i6,i1,i3)
         return
      else
         faaaa=
     *  g(i6,i5,i3,i4,i2,i1)+g(i6,i4,i2,i5,i3,i1)+g(i6,i3,i5,i2,i4,i1)+
     *  g(i6,i2,i4,i3,i5,i1)+g(i5,i6,i3,i4,i1,i2)+g(i5,i4,i1,i6,i3,i2)+
     *  g(i5,i3,i6,i1,i4,i2)+g(i5,i1,i4,i3,i6,i2)+g(i4,i6,i2,i5,i1,i3)+
     *  g(i4,i5,i1,i6,i2,i3)+g(i4,i2,i6,i1,i5,i3)+g(i4,i1,i5,i2,i6,i3)+
     *  g(i3,i6,i5,i2,i1,i4)+g(i3,i5,i6,i1,i2,i4)+g(i3,i2,i1,i6,i5,i4)+
     *  g(i3,i1,i2,i5,i6,i4)+g(i2,i6,i4,i3,i1,i5)+g(i2,i4,i6,i1,i3,i5)+
     *  g(i2,i3,i1,i6,i4,i5)+g(i2,i1,i3,i4,i6,i5)+g(i1,i5,i4,i3,i2,i6)+
     *  g(i1,i4,i5,i2,i3,i6)+g(i1,i3,i2,i5,i4,i6)+g(i1,i2,i3,i4,i5,i6)
         return
      endif
      endif
      endif
      endif
      endif
      endif
      return
      end
**********************************************************************
      subroutine dfaaaa(i1,i2,i3,i4,i5,i6,i7,d12,d13,d14,d23,d24,d34)
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
      if (i7.eq.1) then
         d12=dble(i1)*g(i1-1,i1,i1,i1,i1,i1)
         d13=dble(i1)*g(i1,i1-1,i1,i1,i1,i1)
         d14=dble(i1)*g(i1,i1,i1-1,i1,i1,i1)
         d23=dble(i1)*g(i1,i1,i1,i1-1,i1,i1)
         d24=dble(i1)*g(i1,i1,i1,i1,i1-1,i1)
         d34=dble(i1)*g(i1,i1,i1,i1,i1,i1-1)
         return
      else
      if (i7.eq.3) then
c         faaaa=g(i1,i1,i3,i3,i1,i1)+g(i1,i3,i1,i1,i3,i1)+
c     ,         g(i3,i1,i1,i1,i1,i3)
         d12=dble(i1)*g(i1-1,i1,i3,i3,i1,i1)+
     ,       dble(i1)*g(i1-1,i3,i1,i1,i3,i1)+
     ,       dble(i3)*g(i3-1,i1,i1,i1,i1,i3)
         d13=dble(i1)*g(i1,i1-1,i3,i3,i1,i1)+
     ,       dble(i3)*g(i1,i3-1,i1,i1,i3,i1)+
     ,       dble(i1)*g(i3,i1-1,i1,i1,i1,i3)
         d14=dble(i3)*g(i1,i1,i3-1,i3,i1,i1)+
     ,       dble(i1)*g(i1,i3,i1-1,i1,i3,i1)+
     ,       dble(i1)*g(i3,i1,i1-1,i1,i1,i3)
         d23=dble(i3)*g(i1,i1,i3,i3-1,i1,i1)+
     ,       dble(i1)*g(i1,i3,i1,i1-1,i3,i1)+
     ,       dble(i1)*g(i3,i1,i1,i1-1,i1,i3) 
         d24=dble(i1)*g(i1,i1,i3,i3,i1-1,i1)+
     ,       dble(i3)*g(i1,i3,i1,i1,i3-1,i1)+
     ,       dble(i1)*g(i3,i1,i1,i1,i1-1,i3) 
         d34=dble(i1)*g(i1,i1,i3,i3,i1,i1-1)+
     ,       dble(i1)*g(i1,i3,i1,i1,i3,i1-1)+
     ,       dble(i3)*g(i3,i1,i1,i1,i1,i3-1) 
         return
      else
      if (i7.eq.4) then
c         faaaa=g(i1,i1,i1,i4,i4,i4)+g(i1,i4,i4,i1,i1,i4)+
c     ,         g(i4,i1,i4,i1,i4,i1)+g(i4,i4,i1,i4,i1,i1)
         d12=dble(i1)*g(i1-1,i1,i1,i4,i4,i4)+
     ,       dble(i1)*g(i1-1,i4,i4,i1,i1,i4)+
     ,       dble(i4)*g(i4-1,i1,i4,i1,i4,i1)+
     ,       dble(i4)*g(i4-1,i4,i1,i4,i1,i1)
         d13=dble(i1)*g(i1,i1-1,i1,i4,i4,i4)+
     ,       dble(i4)*g(i1,i4-1,i4,i1,i1,i4)+
     ,       dble(i1)*g(i4,i1-1,i4,i1,i4,i1)+
     ,       dble(i4)*g(i4,i4-1,i1,i4,i1,i1)
         d14=dble(i1)*g(i1,i1,i1-1,i4,i4,i4)+
     ,       dble(i4)*g(i1,i4,i4-1,i1,i1,i4)+
     ,       dble(i4)*g(i4,i1,i4-1,i1,i4,i1)+
     ,       dble(i1)*g(i4,i4,i1-1,i4,i1,i1)
         d23=dble(i4)*g(i1,i1,i1,i4-1,i4,i4)+
     ,       dble(i1)*g(i1,i4,i4,i1-1,i1,i4)+
     ,       dble(i1)*g(i4,i1,i4,i1-1,i4,i1)+
     ,       dble(i4)*g(i4,i4,i1,i4-1,i1,i1)
         d24=dble(i4)*g(i1,i1,i1,i4,i4-1,i4)+
     ,       dble(i1)*g(i1,i4,i4,i1,i1-1,i4)+
     ,       dble(i4)*g(i4,i1,i4,i1,i4-1,i1)+
     ,       dble(i1)*g(i4,i4,i1,i4,i1-1,i1)
         d34=dble(i4)*g(i1,i1,i1,i4,i4,i4-1)+
     ,       dble(i4)*g(i1,i4,i4,i1,i1,i4-1)+
     ,       dble(i1)*g(i4,i1,i4,i1,i4,i1-1)+
     ,       dble(i1)*g(i4,i4,i1,i4,i1,i1-1)
         return
      else
      if (i7.eq.6) then
c         faaaa=g(i1,i2,i3,i4,i5,i6)+g(i1,i3,i2,i5,i4,i6)+
c     ,         g(i3,i1,i2,i5,i6,i4)+g(i4,i5,i6,i1,i2,i3)+
c     ,         g(i5,i4,i6,i1,i3,i2)+g(i5,i6,i4,i3,i1,i2)
         d12=dble(i1)*g(i1-1,i2,i3,i4,i5,i6)+
     ,       dble(i1)*g(i1-1,i3,i2,i5,i4,i6)+
     ,       dble(i3)*g(i3-1,i1,i2,i5,i6,i4)+
     ,       dble(i4)*g(i4-1,i5,i6,i1,i2,i3)+
     ,       dble(i5)*g(i5-1,i4,i6,i1,i3,i2)+
     ,       dble(i5)*g(i5-1,i6,i4,i3,i1,i2)
         d13=dble(i2)*g(i1,i2-1,i3,i4,i5,i6)+
     ,       dble(i3)*g(i1,i3-1,i2,i5,i4,i6)+
     ,       dble(i1)*g(i3,i1-1,i2,i5,i6,i4)+
     ,       dble(i5)*g(i4,i5-1,i6,i1,i2,i3)+
     ,       dble(i4)*g(i5,i4-1,i6,i1,i3,i2)+
     ,       dble(i6)*g(i5,i6-1,i4,i3,i1,i2)
         d14=dble(i3)*g(i1,i2,i3-1,i4,i5,i6)+
     ,       dble(i2)*g(i1,i3,i2-1,i5,i4,i6)+
     ,       dble(i2)*g(i3,i1,i2-1,i5,i6,i4)+
     ,       dble(i6)*g(i4,i5,i6-1,i1,i2,i3)+
     ,       dble(i6)*g(i5,i4,i6-1,i1,i3,i2)+
     ,       dble(i4)*g(i5,i6,i4-1,i3,i1,i2)
         d23=dble(i4)*g(i1,i2,i3,i4-1,i5,i6)+
     ,       dble(i5)*g(i1,i3,i2,i5-1,i4,i6)+
     ,       dble(i5)*g(i3,i1,i2,i5-1,i6,i4)+
     ,       dble(i1)*g(i4,i5,i6,i1-1,i2,i3)+
     ,       dble(i1)*g(i5,i4,i6,i1-1,i3,i2)+
     ,       dble(i3)*g(i5,i6,i4,i3-1,i1,i2)
         d24=dble(i5)*g(i1,i2,i3,i4,i5-1,i6)+
     ,       dble(i4)*g(i1,i3,i2,i5,i4-1,i6)+
     ,       dble(i6)*g(i3,i1,i2,i5,i6-1,i4)+
     ,       dble(i2)*g(i4,i5,i6,i1,i2-1,i3)+
     ,       dble(i3)*g(i5,i4,i6,i1,i3-1,i2)+
     ,       dble(i1)*g(i5,i6,i4,i3,i1-1,i2)
         d34=dble(i6)*g(i1,i2,i3,i4,i5,i6-1)+
     ,       dble(i6)*g(i1,i3,i2,i5,i4,i6-1)+
     ,       dble(i4)*g(i3,i1,i2,i5,i6,i4-1)+
     ,       dble(i3)*g(i4,i5,i6,i1,i2,i3-1)+
     ,       dble(i2)*g(i5,i4,i6,i1,i3,i2-1)+
     ,       dble(i2)*g(i5,i6,i4,i3,i1,i2-1)
         return
      else
      if (i7.eq.11) then
c         faaaa=g(i1,i2,i3,i4,i2,i1)+g(i1,i2,i4,i3,i2,i1)+
c     ,         g(i1,i3,i2,i2,i4,i1)+g(i1,i4,i2,i2,i3,i1)+
c     ,         g(i2,i1,i3,i4,i1,i2)+g(i2,i1,i4,i3,i1,i2)+
c     ,         g(i2,i3,i1,i1,i4,i2)+g(i2,i4,i1,i1,i3,i2)+
c     ,         g(i3,i1,i2,i2,i1,i4)+g(i3,i2,i1,i1,i2,i4)+
c     ,         g(i4,i1,i2,i2,i1,i3)+g(i4,i2,i1,i1,i2,i3)
         d12=dble(i1)*g(i1-1,i2,i3,i4,i2,i1)+
     ,       dble(i1)*g(i1-1,i2,i4,i3,i2,i1)+
     ,       dble(i1)*g(i1-1,i3,i2,i2,i4,i1)+
     ,       dble(i1)*g(i1-1,i4,i2,i2,i3,i1)+
     ,       dble(i2)*g(i2-1,i1,i3,i4,i1,i2)+
     ,       dble(i2)*g(i2-1,i1,i4,i3,i1,i2)+
     ,       dble(i2)*g(i2-1,i3,i1,i1,i4,i2)+
     ,       dble(i2)*g(i2-1,i4,i1,i1,i3,i2)+
     ,       dble(i3)*g(i3-1,i1,i2,i2,i1,i4)+
     ,       dble(i3)*g(i3-1,i2,i1,i1,i2,i4)+
     ,       dble(i4)*g(i4-1,i1,i2,i2,i1,i3)+
     ,       dble(i4)*g(i4-1,i2,i1,i1,i2,i3)
         d13=dble(i2)*g(i1,i2-1,i3,i4,i2,i1)+
     ,       dble(i2)*g(i1,i2-1,i4,i3,i2,i1)+
     ,       dble(i3)*g(i1,i3-1,i2,i2,i4,i1)+
     ,       dble(i4)*g(i1,i4-1,i2,i2,i3,i1)+
     ,       dble(i1)*g(i2,i1-1,i3,i4,i1,i2)+
     ,       dble(i1)*g(i2,i1-1,i4,i3,i1,i2)+
     ,       dble(i3)*g(i2,i3-1,i1,i1,i4,i2)+
     ,       dble(i4)*g(i2,i4-1,i1,i1,i3,i2)+
     ,       dble(i1)*g(i3,i1-1,i2,i2,i1,i4)+
     ,       dble(i2)*g(i3,i2-1,i1,i1,i2,i4)+
     ,       dble(i1)*g(i4,i1-1,i2,i2,i1,i3)+
     ,       dble(i2)*g(i4,i2-1,i1,i1,i2,i3)
         d14=dble(i3)*g(i1,i2,i3-1,i4,i2,i1)+
     ,       dble(i4)*g(i1,i2,i4-1,i3,i2,i1)+
     ,       dble(i2)*g(i1,i3,i2-1,i2,i4,i1)+
     ,       dble(i2)*g(i1,i4,i2-1,i2,i3,i1)+
     ,       dble(i3)*g(i2,i1,i3-1,i4,i1,i2)+
     ,       dble(i4)*g(i2,i1,i4-1,i3,i1,i2)+
     ,       dble(i1)*g(i2,i3,i1-1,i1,i4,i2)+
     ,       dble(i1)*g(i2,i4,i1-1,i1,i3,i2)+
     ,       dble(i2)*g(i3,i1,i2-1,i2,i1,i4)+
     ,       dble(i1)*g(i3,i2,i1-1,i1,i2,i4)+
     ,       dble(i2)*g(i4,i1,i2-1,i2,i1,i3)+
     ,       dble(i1)*g(i4,i2,i1-1,i1,i2,i3)
         d23=dble(i4)*g(i1,i2,i3,i4-1,i2,i1)+
     ,       dble(i3)*g(i1,i2,i4,i3-1,i2,i1)+
     ,       dble(i2)*g(i1,i3,i2,i2-1,i4,i1)+
     ,       dble(i2)*g(i1,i4,i2,i2-1,i3,i1)+
     ,       dble(i4)*g(i2,i1,i3,i4-1,i1,i2)+
     ,       dble(i3)*g(i2,i1,i4,i3-1,i1,i2)+
     ,       dble(i1)*g(i2,i3,i1,i1-1,i4,i2)+
     ,       dble(i1)*g(i2,i4,i1,i1-1,i3,i2)+
     ,       dble(i2)*g(i3,i1,i2,i2-1,i1,i4)+
     ,       dble(i1)*g(i3,i2,i1,i1-1,i2,i4)+
     ,       dble(i2)*g(i4,i1,i2,i2-1,i1,i3)+
     ,       dble(i1)*g(i4,i2,i1,i1-1,i2,i3)
         d24=dble(i2)*g(i1,i2,i3,i4,i2-1,i1)+
     ,       dble(i2)*g(i1,i2,i4,i3,i2-1,i1)+
     ,       dble(i4)*g(i1,i3,i2,i2,i4-1,i1)+
     ,       dble(i3)*g(i1,i4,i2,i2,i3-1,i1)+
     ,       dble(i1)*g(i2,i1,i3,i4,i1-1,i2)+
     ,       dble(i1)*g(i2,i1,i4,i3,i1-1,i2)+
     ,       dble(i4)*g(i2,i3,i1,i1,i4-1,i2)+
     ,       dble(i3)*g(i2,i4,i1,i1,i3-1,i2)+
     ,       dble(i1)*g(i3,i1,i2,i2,i1-1,i4)+
     ,       dble(i2)*g(i3,i2,i1,i1,i2-1,i4)+
     ,       dble(i1)*g(i4,i1,i2,i2,i1-1,i3)+
     ,       dble(i2)*g(i4,i2,i1,i1,i2-1,i3)
         d34=dble(i1)*g(i1,i2,i3,i4,i2,i1-1)+
     ,       dble(i1)*g(i1,i2,i4,i3,i2,i1-1)+
     ,       dble(i1)*g(i1,i3,i2,i2,i4,i1-1)+
     ,       dble(i1)*g(i1,i4,i2,i2,i3,i1-1)+
     ,       dble(i2)*g(i2,i1,i3,i4,i1,i2-1)+
     ,       dble(i2)*g(i2,i1,i4,i3,i1,i2-1)+
     ,       dble(i2)*g(i2,i3,i1,i1,i4,i2-1)+
     ,       dble(i2)*g(i2,i4,i1,i1,i3,i2-1)+
     ,       dble(i4)*g(i3,i1,i2,i2,i1,i4-1)+
     ,       dble(i4)*g(i3,i2,i1,i1,i2,i4-1)+
     ,       dble(i3)*g(i4,i1,i2,i2,i1,i3-1)+
     ,       dble(i3)*g(i4,i2,i1,i1,i2,i3-1)
         return
      else
      if (i7.eq.13) then
c         faaaa=g(i1,i1,i3,i4,i6,i6)+g(i1,i6,i4,i3,i1,i6)+
c     ,         g(i1,i3,i1,i6,i4,i6)+g(i1,i4,i6,i1,i3,i6)+
c     ,         g(i6,i1,i4,i3,i6,i1)+g(i6,i6,i3,i4,i1,i1)+
c     ,         g(i6,i3,i6,i1,i4,i1)+g(i6,i4,i1,i6,i3,i1)+
c     ,         g(i3,i1,i1,i6,i6,i4)+g(i3,i6,i6,i1,i1,i4)+
c     ,         g(i4,i1,i6,i1,i6,i3)+g(i4,i6,i1,i6,i1,i3)
        d12=dble(i1)*g(i1-1,i1,i3,i4,i6,i6)+
     ,      dble(i1)*g(i1-1,i6,i4,i3,i1,i6)+
     ,      dble(i1)*g(i1-1,i3,i1,i6,i4,i6)+
     ,      dble(i1)*g(i1-1,i4,i6,i1,i3,i6)+
     ,      dble(i6)*g(i6-1,i1,i4,i3,i6,i1)+
     ,      dble(i6)*g(i6-1,i6,i3,i4,i1,i1)+
     ,      dble(i6)*g(i6-1,i3,i6,i1,i4,i1)+
     ,      dble(i6)*g(i6-1,i4,i1,i6,i3,i1)+
     ,      dble(i3)*g(i3-1,i1,i1,i6,i6,i4)+
     ,      dble(i3)*g(i3-1,i6,i6,i1,i1,i4)+
     ,      dble(i4)*g(i4-1,i1,i6,i1,i6,i3)+
     ,      dble(i4)*g(i4-1,i6,i1,i6,i1,i3)
        d13=dble(i1)*g(i1,i1-1,i3,i4,i6,i6)+
     ,      dble(i6)*g(i1,i6-1,i4,i3,i1,i6)+
     ,      dble(i3)*g(i1,i3-1,i1,i6,i4,i6)+
     ,      dble(i4)*g(i1,i4-1,i6,i1,i3,i6)+
     ,      dble(i1)*g(i6,i1-1,i4,i3,i6,i1)+
     ,      dble(i6)*g(i6,i6-1,i3,i4,i1,i1)+
     ,      dble(i3)*g(i6,i3-1,i6,i1,i4,i1)+
     ,      dble(i4)*g(i6,i4-1,i1,i6,i3,i1)+
     ,      dble(i1)*g(i3,i1-1,i1,i6,i6,i4)+
     ,      dble(i6)*g(i3,i6-1,i6,i1,i1,i4)+
     ,      dble(i1)*g(i4,i1-1,i6,i1,i6,i3)+
     ,      dble(i6)*g(i4,i6-1,i1,i6,i1,i3)
        d14=dble(i3)*g(i1,i1,i3-1,i4,i6,i6)+
     ,      dble(i4)*g(i1,i6,i4-1,i3,i1,i6)+
     ,      dble(i1)*g(i1,i3,i1-1,i6,i4,i6)+
     ,      dble(i6)*g(i1,i4,i6-1,i1,i3,i6)+
     ,      dble(i4)*g(i6,i1,i4-1,i3,i6,i1)+
     ,      dble(i3)*g(i6,i6,i3-1,i4,i1,i1)+
     ,      dble(i6)*g(i6,i3,i6-1,i1,i4,i1)+
     ,      dble(i1)*g(i6,i4,i1-1,i6,i3,i1)+
     ,      dble(i1)*g(i3,i1,i1-1,i6,i6,i4)+
     ,      dble(i6)*g(i3,i6,i6-1,i1,i1,i4)+
     ,      dble(i6)*g(i4,i1,i6-1,i1,i6,i3)+
     ,      dble(i1)*g(i4,i6,i1-1,i6,i1,i3)
        d23=dble(i4)*g(i1,i1,i3,i4-1,i6,i6)+
     ,      dble(i3)*g(i1,i6,i4,i3-1,i1,i6)+
     ,      dble(i6)*g(i1,i3,i1,i6-1,i4,i6)+
     ,      dble(i1)*g(i1,i4,i6,i1-1,i3,i6)+
     ,      dble(i3)*g(i6,i1,i4,i3-1,i6,i1)+
     ,      dble(i4)*g(i6,i6,i3,i4-1,i1,i1)+
     ,      dble(i1)*g(i6,i3,i6,i1-1,i4,i1)+
     ,      dble(i6)*g(i6,i4,i1,i6-1,i3,i1)+
     ,      dble(i6)*g(i3,i1,i1,i6-1,i6,i4)+
     ,      dble(i1)*g(i3,i6,i6,i1-1,i1,i4)+
     ,      dble(i1)*g(i4,i1,i6,i1-1,i6,i3)+
     ,      dble(i6)*g(i4,i6,i1,i6-1,i1,i3)
        d24=dble(i6)*g(i1,i1,i3,i4,i6-1,i6)+
     ,      dble(i1)*g(i1,i6,i4,i3,i1-1,i6)+
     ,      dble(i4)*g(i1,i3,i1,i6,i4-1,i6)+
     ,      dble(i3)*g(i1,i4,i6,i1,i3-1,i6)+
     ,      dble(i6)*g(i6,i1,i4,i3,i6-1,i1)+
     ,      dble(i1)*g(i6,i6,i3,i4,i1-1,i1)+
     ,      dble(i4)*g(i6,i3,i6,i1,i4-1,i1)+
     ,      dble(i3)*g(i6,i4,i1,i6,i3-1,i1)+
     ,      dble(i6)*g(i3,i1,i1,i6,i6-1,i4)+
     ,      dble(i1)*g(i3,i6,i6,i1,i1-1,i4)+
     ,      dble(i6)*g(i4,i1,i6,i1,i6-1,i3)+
     ,      dble(i1)*g(i4,i6,i1,i6,i1-1,i3)
        d34=dble(i6)*g(i1,i1,i3,i4,i6,i6-1)+
     ,      dble(i6)*g(i1,i6,i4,i3,i1,i6-1)+
     ,      dble(i6)*g(i1,i3,i1,i6,i4,i6-1)+
     ,      dble(i6)*g(i1,i4,i6,i1,i3,i6-1)+
     ,      dble(i1)*g(i6,i1,i4,i3,i6,i1-1)+
     ,      dble(i1)*g(i6,i6,i3,i4,i1,i1-1)+
     ,      dble(i1)*g(i6,i3,i6,i1,i4,i1-1)+
     ,      dble(i1)*g(i6,i4,i1,i6,i3,i1-1)+
     ,      dble(i4)*g(i3,i1,i1,i6,i6,i4-1)+
     ,      dble(i4)*g(i3,i6,i6,i1,i1,i4-1)+
     ,      dble(i3)*g(i4,i1,i6,i1,i6,i3-1)+
     ,      dble(i3)*g(i4,i6,i1,i6,i1,i3-1)
         return
      else
c         faaaa=
c     *  g(i6,i5,i3,i4,i2,i1)+g(i6,i4,i2,i5,i3,i1)+g(i6,i3,i5,i2,i4,i1)+
c     *  g(i6,i2,i4,i3,i5,i1)+g(i5,i6,i3,i4,i1,i2)+g(i5,i4,i1,i6,i3,i2)+
c     *  g(i5,i3,i6,i1,i4,i2)+g(i5,i1,i4,i3,i6,i2)+g(i4,i6,i2,i5,i1,i3)+
c     *  g(i4,i5,i1,i6,i2,i3)+g(i4,i2,i6,i1,i5,i3)+g(i4,i1,i5,i2,i6,i3)+
c     *  g(i3,i6,i5,i2,i1,i4)+g(i3,i5,i6,i1,i2,i4)+g(i3,i2,i1,i6,i5,i4)+
c     *  g(i3,i1,i2,i5,i6,i4)+g(i2,i6,i4,i3,i1,i5)+g(i2,i4,i6,i1,i3,i5)+
c     *  g(i2,i3,i1,i6,i4,i5)+g(i2,i1,i3,i4,i6,i5)+g(i1,i5,i4,i3,i2,i6)+
c     *  g(i1,i4,i5,i2,i3,i6)+g(i1,i3,i2,i5,i4,i6)+g(i1,i2,i3,i4,i5,i6)
       d12=dble(i6)*g(i6-1,i5,i3,i4,i2,i1)+
     ,     dble(i6)*g(i6-1,i4,i2,i5,i3,i1)+
     ,     dble(i6)*g(i6-1,i3,i5,i2,i4,i1)+
     ,     dble(i6)*g(i6-1,i2,i4,i3,i5,i1)+
     ,     dble(i5)*g(i5-1,i6,i3,i4,i1,i2)+
     ,     dble(i5)*g(i5-1,i4,i1,i6,i3,i2)+
     ,     dble(i5)*g(i5-1,i3,i6,i1,i4,i2)+
     ,     dble(i5)*g(i5-1,i1,i4,i3,i6,i2)+
     ,     dble(i4)*g(i4-1,i6,i2,i5,i1,i3)+
     ,     dble(i4)*g(i4-1,i5,i1,i6,i2,i3)+
     ,     dble(i4)*g(i4-1,i2,i6,i1,i5,i3)+
     ,     dble(i4)*g(i4-1,i1,i5,i2,i6,i3)+
     ,     dble(i3)*g(i3-1,i6,i5,i2,i1,i4)+
     ,     dble(i3)*g(i3-1,i5,i6,i1,i2,i4)+
     ,     dble(i3)*g(i3-1,i2,i1,i6,i5,i4)+
     ,     dble(i3)*g(i3-1,i1,i2,i5,i6,i4)+
     ,     dble(i2)*g(i2-1,i6,i4,i3,i1,i5)+
     ,     dble(i2)*g(i2-1,i4,i6,i1,i3,i5)+
     ,     dble(i2)*g(i2-1,i3,i1,i6,i4,i5)+
     ,     dble(i2)*g(i2-1,i1,i3,i4,i6,i5)+
     ,     dble(i1)*g(i1-1,i5,i4,i3,i2,i6)+
     ,     dble(i1)*g(i1-1,i4,i5,i2,i3,i6)+
     ,     dble(i1)*g(i1-1,i3,i2,i5,i4,i6)+
     ,     dble(i1)*g(i1-1,i2,i3,i4,i5,i6)
       d13=dble(i5)*g(i6,i5-1,i3,i4,i2,i1)+
     ,     dble(i4)*g(i6,i4-1,i2,i5,i3,i1)+
     ,     dble(i3)*g(i6,i3-1,i5,i2,i4,i1)+
     ,     dble(i2)*g(i6,i2-1,i4,i3,i5,i1)+
     ,     dble(i6)*g(i5,i6-1,i3,i4,i1,i2)+
     ,     dble(i4)*g(i5,i4-1,i1,i6,i3,i2)+
     ,     dble(i3)*g(i5,i3-1,i6,i1,i4,i2)+
     ,     dble(i1)*g(i5,i1-1,i4,i3,i6,i2)+
     ,     dble(i6)*g(i4,i6-1,i2,i5,i1,i3)+
     ,     dble(i5)*g(i4,i5-1,i1,i6,i2,i3)+
     ,     dble(i2)*g(i4,i2-1,i6,i1,i5,i3)+
     ,     dble(i1)*g(i4,i1-1,i5,i2,i6,i3)+
     ,     dble(i6)*g(i3,i6-1,i5,i2,i1,i4)+
     ,     dble(i5)*g(i3,i5-1,i6,i1,i2,i4)+
     ,     dble(i2)*g(i3,i2-1,i1,i6,i5,i4)+
     ,     dble(i1)*g(i3,i1-1,i2,i5,i6,i4)+
     ,     dble(i6)*g(i2,i6-1,i4,i3,i1,i5)+
     ,     dble(i4)*g(i2,i4-1,i6,i1,i3,i5)+
     ,     dble(i3)*g(i2,i3-1,i1,i6,i4,i5)+
     ,     dble(i1)*g(i2,i1-1,i3,i4,i6,i5)+
     ,     dble(i5)*g(i1,i5-1,i4,i3,i2,i6)+
     ,     dble(i4)*g(i1,i4-1,i5,i2,i3,i6)+
     ,     dble(i3)*g(i1,i3-1,i2,i5,i4,i6)+
     ,     dble(i2)*g(i1,i2-1,i3,i4,i5,i6)
       d14=dble(i3)*g(i6,i5,i3-1,i4,i2,i1)+
     ,     dble(i2)*g(i6,i4,i2-1,i5,i3,i1)+
     ,     dble(i5)*g(i6,i3,i5-1,i2,i4,i1)+
     ,     dble(i4)*g(i6,i2,i4-1,i3,i5,i1)+
     ,     dble(i3)*g(i5,i6,i3-1,i4,i1,i2)+
     ,     dble(i1)*g(i5,i4,i1-1,i6,i3,i2)+
     ,     dble(i6)*g(i5,i3,i6-1,i1,i4,i2)+
     ,     dble(i4)*g(i5,i1,i4-1,i3,i6,i2)+
     ,     dble(i2)*g(i4,i6,i2-1,i5,i1,i3)+
     ,     dble(i1)*g(i4,i5,i1-1,i6,i2,i3)+
     ,     dble(i6)*g(i4,i2,i6-1,i1,i5,i3)+
     ,     dble(i5)*g(i4,i1,i5-1,i2,i6,i3)+
     ,     dble(i5)*g(i3,i6,i5-1,i2,i1,i4)+
     ,     dble(i6)*g(i3,i5,i6-1,i1,i2,i4)+
     ,     dble(i1)*g(i3,i2,i1-1,i6,i5,i4)+
     ,     dble(i2)*g(i3,i1,i2-1,i5,i6,i4)+
     ,     dble(i4)*g(i2,i6,i4-1,i3,i1,i5)+
     ,     dble(i6)*g(i2,i4,i6-1,i1,i3,i5)+
     ,     dble(i1)*g(i2,i3,i1-1,i6,i4,i5)+
     ,     dble(i3)*g(i2,i1,i3-1,i4,i6,i5)+
     ,     dble(i4)*g(i1,i5,i4-1,i3,i2,i6)+
     ,     dble(i5)*g(i1,i4,i5-1,i2,i3,i6)+
     ,     dble(i2)*g(i1,i3,i2-1,i5,i4,i6)+
     ,     dble(i3)*g(i1,i2,i3-1,i4,i5,i6)
       d23=dble(i4)*g(i6,i5,i3,i4-1,i2,i1)+
     ,     dble(i5)*g(i6,i4,i2,i5-1,i3,i1)+
     ,     dble(i2)*g(i6,i3,i5,i2-1,i4,i1)+
     ,     dble(i3)*g(i6,i2,i4,i3-1,i5,i1)+
     ,     dble(i4)*g(i5,i6,i3,i4-1,i1,i2)+
     ,     dble(i6)*g(i5,i4,i1,i6-1,i3,i2)+
     ,     dble(i1)*g(i5,i3,i6,i1-1,i4,i2)+
     ,     dble(i3)*g(i5,i1,i4,i3-1,i6,i2)+
     ,     dble(i5)*g(i4,i6,i2,i5-1,i1,i3)+
     ,     dble(i6)*g(i4,i5,i1,i6-1,i2,i3)+
     ,     dble(i1)*g(i4,i2,i6,i1-1,i5,i3)+
     ,     dble(i2)*g(i4,i1,i5,i2-1,i6,i3)+
     ,     dble(i2)*g(i3,i6,i5,i2-1,i1,i4)+
     ,     dble(i1)*g(i3,i5,i6,i1-1,i2,i4)+
     ,     dble(i6)*g(i3,i2,i1,i6-1,i5,i4)+
     ,     dble(i5)*g(i3,i1,i2,i5-1,i6,i4)+
     ,     dble(i3)*g(i2,i6,i4,i3-1,i1,i5)+
     ,     dble(i1)*g(i2,i4,i6,i1-1,i3,i5)+
     ,     dble(i6)*g(i2,i3,i1,i6-1,i4,i5)+
     ,     dble(i4)*g(i2,i1,i3,i4-1,i6,i5)+
     ,     dble(i3)*g(i1,i5,i4,i3-1,i2,i6)+
     ,     dble(i2)*g(i1,i4,i5,i2-1,i3,i6)+
     ,     dble(i5)*g(i1,i3,i2,i5-1,i4,i6)+
     ,     dble(i4)*g(i1,i2,i3,i4-1,i5,i6)
       d24=dble(i2)*g(i6,i5,i3,i4,i2-1,i1)+
     ,     dble(i3)*g(i6,i4,i2,i5,i3-1,i1)+
     ,     dble(i4)*g(i6,i3,i5,i2,i4-1,i1)+
     ,     dble(i5)*g(i6,i2,i4,i3,i5-1,i1)+
     ,     dble(i1)*g(i5,i6,i3,i4,i1-1,i2)+
     ,     dble(i3)*g(i5,i4,i1,i6,i3-1,i2)+
     ,     dble(i4)*g(i5,i3,i6,i1,i4-1,i2)+
     ,     dble(i6)*g(i5,i1,i4,i3,i6-1,i2)+
     ,     dble(i1)*g(i4,i6,i2,i5,i1-1,i3)+
     ,     dble(i2)*g(i4,i5,i1,i6,i2-1,i3)+
     ,     dble(i5)*g(i4,i2,i6,i1,i5-1,i3)+
     ,     dble(i6)*g(i4,i1,i5,i2,i6-1,i3)+
     ,     dble(i1)*g(i3,i6,i5,i2,i1-1,i4)+
     ,     dble(i2)*g(i3,i5,i6,i1,i2-1,i4)+
     ,     dble(i5)*g(i3,i2,i1,i6,i5-1,i4)+
     ,     dble(i6)*g(i3,i1,i2,i5,i6-1,i4)+
     ,     dble(i1)*g(i2,i6,i4,i3,i1-1,i5)+
     ,     dble(i3)*g(i2,i4,i6,i1,i3-1,i5)+
     ,     dble(i4)*g(i2,i3,i1,i6,i4-1,i5)+
     ,     dble(i6)*g(i2,i1,i3,i4,i6-1,i5)+
     ,     dble(i2)*g(i1,i5,i4,i3,i2-1,i6)+
     ,     dble(i3)*g(i1,i4,i5,i2,i3-1,i6)+
     ,     dble(i4)*g(i1,i3,i2,i5,i4-1,i6)+
     ,     dble(i5)*g(i1,i2,i3,i4,i5-1,i6)
       d34=dble(i1)*g(i6,i5,i3,i4,i2,i1-1)+
     ,     dble(i1)*g(i6,i4,i2,i5,i3,i1-1)+
     ,     dble(i1)*g(i6,i3,i5,i2,i4,i1-1)+
     ,     dble(i1)*g(i6,i2,i4,i3,i5,i1-1)+
     ,     dble(i2)*g(i5,i6,i3,i4,i1,i2-1)+
     ,     dble(i2)*g(i5,i4,i1,i6,i3,i2-1)+
     ,     dble(i2)*g(i5,i3,i6,i1,i4,i2-1)+
     ,     dble(i2)*g(i5,i1,i4,i3,i6,i2-1)+
     ,     dble(i3)*g(i4,i6,i2,i5,i1,i3-1)+
     ,     dble(i3)*g(i4,i5,i1,i6,i2,i3-1)+
     ,     dble(i3)*g(i4,i2,i6,i1,i5,i3-1)+
     ,     dble(i3)*g(i4,i1,i5,i2,i6,i3-1)+
     ,     dble(i4)*g(i3,i6,i5,i2,i1,i4-1)+
     ,     dble(i4)*g(i3,i5,i6,i1,i2,i4-1)+
     ,     dble(i4)*g(i3,i2,i1,i6,i5,i4-1)+
     ,     dble(i4)*g(i3,i1,i2,i5,i6,i4-1)+
     ,     dble(i5)*g(i2,i6,i4,i3,i1,i5-1)+
     ,     dble(i5)*g(i2,i4,i6,i1,i3,i5-1)+
     ,     dble(i5)*g(i2,i3,i1,i6,i4,i5-1)+
     ,     dble(i5)*g(i2,i1,i3,i4,i6,i5-1)+
     ,     dble(i6)*g(i1,i5,i4,i3,i2,i6-1)+
     ,     dble(i6)*g(i1,i4,i5,i2,i3,i6-1)+
     ,     dble(i6)*g(i1,i3,i2,i5,i4,i6-1)+
     ,     dble(i6)*g(i1,i2,i3,i4,i5,i6-1)
         return
      endif
      endif
      endif
      endif
      endif
      endif
      return
      end
****************************************************************
      subroutine saaaa(nombre,n,mt,vex1,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie AAAA
c     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*4 nombre
      include 'dimensions.inc'
      common/ipot41/ipa41(nexp41,7),nfas41(igrad4)
      common sm(nexp)
      is=nfas41(mt)
      write(7,1) nombre,n,rms*627.51d0,emax*627.51d0,
     -is,mt-2,vex1
      do 2 k=1,is
         write(7,3) k,sm(k)
2     continue
      do k=1,is
         write(7,5) k,ipa41(k,1),k,ipa41(k,2),k,ipa41(k,3),
     -         k,ipa41(k,4),k,ipa41(k,5),k,ipa41(k,6),k,ipa41(k,7)
      enddo
!      write(7,4) vex1
!4     format(6x,'vex1=',d23.15)
      write(7,10)
      write(7,11)
      write(7,12)
      write(7,13) mt-2
      write(7,14) mt-2
1     format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 6D PES      ',/,
     -'*     for the AAAA system class fitted to ',i4,' points      ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_2=',i2,')                      ',/,
     -'      parameter(vex1=',d23.15,')                             ',/,
     -'      dimension der(6)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),',
     ,'i7(is),cf(is)',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/         ',/,
     -'      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/            ',/,
     -'      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            ')
3     format(
     -'      data cf(',i3,')/',d23.15,'/                            ')
5     format(
     -'      data i1(',i3,')/',i2,'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,
     -'/,i4(',i3,')/',i2,'/   ',/,
     -'      data i5(',i3,')/',i2,'/,i6(',i3,')/',i2,'/,',
     -'i7(',i3,')/',i2,'/')
10    format(
     -'      ener=0.d0          ',/,
     -'      der12=0.d0         ',/,
     -'      der13=0.d0         ',/,
     -'      der14=0.d0         ',/,
     -'      der23=0.d0         ',/,
     -'      der24=0.d0         ',/,
     -'      der34=0.d0         ',/,
     -'      pux12= vex1*r12 ',/,
     -'      pux13= vex1*r13 ',/,
     -'      pux14= vex1*r14 ',/,
     -'      pux23= vex1*r23 ',/,
     -'      pux24= vex1*r24 ',/,
     -'      pux34= vex1*r34 ',/,
     -'      qux12= dexp(-pux12) ',/,
     -'      qux13= dexp(-pux13) ',/,
     -'      qux14= dexp(-pux14) ',/,
     -'      qux23= dexp(-pux23) ',/,
     -'      qux24= dexp(-pux24) ',/,
     -'      qux34= dexp(-pux34) ',/,
     -'      aux12= r12*qux12 ',/,
     -'      aux13= r13*qux13 ',/,
     -'      aux14= r14*qux14 ',/,
     -'      aux23= r23*qux23 ',/,
     -'      aux24= r24*qux24 ',/,
     -'      aux34= r34*qux34 ')
11    format(
     -'      do 1 i=1,mt_2               ',/,
     -'         f12(i) = aux12*f12(i-1)  ',/,
     -'         f13(i) = aux13*f13(i-1)  ',/,
     -'         f14(i) = aux14*f14(i-1)  ',/,
     -'         f23(i) = aux23*f23(i-1)  ',/,
     -'         f24(i) = aux24*f24(i-1)  ',/,
     -'         f34(i) = aux34*f34(i-1)  ',/,
     -'1     continue                 ')
12    format(
     -'      do 2 j=1,is     ',/,
     -'         au=faaaa(i1(j),i2(j),i3(j),i4(j),i5(j),i6(j),i7(j))',/,
     -'         ener=ener+cf(j)*au                                 ',/,
     -'         if (iop.eq.1) then                                 ',/,
     -'            call dfaaaa(                                    ',/,
     -'     ,            i1(j),i2(j),i3(j),i4(j),i5(j),i6(j),i7(j),',/,
     -'     ,            d12,d13,d14,d23,d24,d34)      ',/,
     -'            der12=der12+cf(j)*d12                          ',/,
     -'            der13=der13+cf(j)*d13                          ',/,
     -'            der14=der14+cf(j)*d14                          ',/,
     -'            der23=der23+cf(j)*d23                          ',/,
     -'            der24=der24+cf(j)*d24                          ',/,
     -'            der34=der34+cf(j)*d34                          ',/,
     -'         endif                                             ',/,
     -'2     continue                                              ',/,
     -'      if (iop.eq.1) then                 ',/,
     -'         der(1)=der12*(1.d0-pux12)*qux12    ',/,
     -'         der(2)=der13*(1.d0-pux13)*qux13    ',/,
     -'         der(3)=der14*(1.d0-pux14)*qux14    ',/,
     -'         der(4)=der23*(1.d0-pux23)*qux23    ',/,
     -'         der(5)=der24*(1.d0-pux24)*qux24    ',/,
     -'         der(6)=der34*(1.d0-pux34)*qux34    ',/,
     -'      endif                                 ',/,
     -'      return                                                ',/,
     -'      end                                                   ')
**********************************************************************
13    format(
     -'*************************************************************',/,
     -'      function faaaa(i1,i2,i3,i4,i5,i6,i7)                   ',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(mt_2=',i2,')                      ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
     -',/,
     -'      if (i7.eq.1) then                                      ',/,
     -'         faaaa=g(i1,i1,i1,i1,i1,i1)                          ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.3) then                                  ',/,
     -'         faaaa=g(i1,i1,i3,i3,i1,i1)+g(i1,i3,i1,i1,i3,i1)+    ',/,
     -'     ,         g(i3,i1,i1,i1,i1,i3)                          ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.4) then                                  ',/,
     -'         faaaa=g(i1,i1,i1,i4,i4,i4)+g(i1,i4,i4,i1,i1,i4)+    ',/,
     -'     ,         g(i4,i1,i4,i1,i4,i1)+g(i4,i4,i1,i4,i1,i1)     ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.6) then                                  ',/,
     -'         faaaa=g(i1,i2,i3,i4,i5,i6)+g(i1,i3,i2,i5,i4,i6)+    ',/,
     -'     ,         g(i3,i1,i2,i5,i6,i4)+g(i4,i5,i6,i1,i2,i3)+    ',/,
     -'     ,         g(i5,i4,i6,i1,i3,i2)+g(i5,i6,i4,i3,i1,i2)     ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.11) then                                 ',/,
     -'         faaaa=g(i1,i2,i3,i4,i2,i1)+g(i1,i2,i4,i3,i2,i1)+    ',/,
     -'     ,         g(i1,i3,i2,i2,i4,i1)+g(i1,i4,i2,i2,i3,i1)+    ',/,
     -'     ,         g(i2,i1,i3,i4,i1,i2)+g(i2,i1,i4,i3,i1,i2)+    ',/,
     -'     ,         g(i2,i3,i1,i1,i4,i2)+g(i2,i4,i1,i1,i3,i2)+    ',/,
     -'     ,         g(i3,i1,i2,i2,i1,i4)+g(i3,i2,i1,i1,i2,i4)+    ',/,
     -'     ,         g(i4,i1,i2,i2,i1,i3)+g(i4,i2,i1,i1,i2,i3)     ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.13) then                                 ',/,
     -'         faaaa=g(i1,i1,i3,i4,i6,i6)+g(i1,i6,i4,i3,i1,i6)+    ',/,
     -'     ,         g(i1,i3,i1,i6,i4,i6)+g(i1,i4,i6,i1,i3,i6)+    ',/,
     -'     ,         g(i6,i1,i4,i3,i6,i1)+g(i6,i6,i3,i4,i1,i1)+    ',/,
     -'     ,         g(i6,i3,i6,i1,i4,i1)+g(i6,i4,i1,i6,i3,i1)+    ',/,
     -'     ,         g(i3,i1,i1,i6,i6,i4)+g(i3,i6,i6,i1,i1,i4)+    ',/,
     -'     ,         g(i4,i1,i6,i1,i6,i3)+g(i4,i6,i1,i6,i1,i3)     ',/,
     -'         return                                              ',/,
     -'      else                                                   ',/,
     -'         faaaa=                                              ',/,
     -'     *  g(i6,i5,i3,i4,i2,i1)+g(i6,i4,i2,i5,i3,i1)+'
     -,'g(i6,i3,i5,i2,i4,i1)+                                       ',/,
     -'     *  g(i6,i2,i4,i3,i5,i1)+g(i5,i6,i3,i4,i1,i2)+'
     -,'g(i5,i4,i1,i6,i3,i2)+                                       ',/,
     -'     *  g(i5,i3,i6,i1,i4,i2)+g(i5,i1,i4,i3,i6,i2)+'
     -,'g(i4,i6,i2,i5,i1,i3)+                                       ',/,
     -'     *  g(i4,i5,i1,i6,i2,i3)+g(i4,i2,i6,i1,i5,i3)+'
     -,'g(i4,i1,i5,i2,i6,i3)+                                       ',/,
     -'     *  g(i3,i6,i5,i2,i1,i4)+g(i3,i5,i6,i1,i2,i4)+'
     -,'g(i3,i2,i1,i6,i5,i4)+                                       ',/,
     -'     *  g(i3,i1,i2,i5,i6,i4)+g(i2,i6,i4,i3,i1,i5)+'
     -,'g(i2,i4,i6,i1,i3,i5)+                                       ',/,
     -'     *  g(i2,i3,i1,i6,i4,i5)+g(i2,i1,i3,i4,i6,i5)+'
     -,'g(i1,i5,i4,i3,i2,i6)+                                       ',/,
     -'     *  g(i1,i4,i5,i2,i3,i6)+g(i1,i3,i2,i5,i4,i6)+'
     -,'g(i1,i2,i3,i4,i5,i6)                                        ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      return                                                 ',/,
     -'      end                                                    ')
14    format(
     -'*************************************************************',/,
     -'      subroutine dfaaaa(i1,i2,i3,i4,i5,i6,i7,                ',/,
     -'     ,            d12,d13,d14,d23,d24,d34)                   ',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(mt_2=',i2,')                      ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
     -',/,
     -'      if (i7.eq.1) then                                      ',/,
     -'!        faaaa=g(i1,i1,i1,i1,i1,i1)                          ',/,
     -'         d12=i1*g(i1-1,i1,i1,i1,i1,i1)                       ',/,
     -'         d13=i1*g(i1,i1-1,i1,i1,i1,i1)                       ',/,
     -'         d14=i1*g(i1,i1,i1-1,i1,i1,i1)                       ',/,
     -'         d23=i1*g(i1,i1,i1,i1-1,i1,i1)                       ',/,
     -'         d24=i1*g(i1,i1,i1,i1,i1-1,i1)                       ',/,
     -'         d34=i1*g(i1,i1,i1,i1,i1,i1-1)                       ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.3) then                                  ',/,
     -'!        faaaa=g(i1,i1,i3,i3,i1,i1)+g(i1,i3,i1,i1,i3,i1)+    ',/,
     -'!    ,         g(i3,i1,i1,i1,i1,i3)                          ',/,
     -'         d12=i1*g(i1-1,i1,i3,i3,i1,i1)+                      ',/,
     -'     ,       i1*g(i1-1,i3,i1,i1,i3,i1)+                      ',/,
     -'     ,       i3*g(i3-1,i1,i1,i1,i1,i3)                       ',/,
     -'         d13=i1*g(i1,i1-1,i3,i3,i1,i1)+                      ',/,
     -'     ,       i3*g(i1,i3-1,i1,i1,i3,i1)+                      ',/,
     -'     ,       i1*g(i3,i1-1,i1,i1,i1,i3)                       ',/,
     -'         d14=i3*g(i1,i1,i3-1,i3,i1,i1)+                      ',/,
     -'     ,       i1*g(i1,i3,i1-1,i1,i3,i1)+                      ',/,
     -'     ,       i1*g(i3,i1,i1-1,i1,i1,i3)                       ',/,
     -'         d23=i3*g(i1,i1,i3,i3-1,i1,i1)+                      ',/,
     -'     ,       i1*g(i1,i3,i1,i1-1,i3,i1)+                      ',/,
     -'     ,       i1*g(i3,i1,i1,i1-1,i1,i3)                       ',/,
     -'         d24=i1*g(i1,i1,i3,i3,i1-1,i1)+                      ',/,
     -'     ,       i3*g(i1,i3,i1,i1,i3-1,i1)+                      ',/,
     -'     ,       i1*g(i3,i1,i1,i1,i1-1,i3)                       ',/,
     -'         d34=i1*g(i1,i1,i3,i3,i1,i1-1)+                      ',/,
     -'     ,       i1*g(i1,i3,i1,i1,i3,i1-1)+                      ',/,
     -'     ,       i3*g(i3,i1,i1,i1,i1,i3-1)                       ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.4) then                                  ',/,
     -'!        faaaa=g(i1,i1,i1,i4,i4,i4)+g(i1,i4,i4,i1,i1,i4)+    ',/,
     -'!    ,         g(i4,i1,i4,i1,i4,i1)+g(i4,i4,i1,i4,i1,i1)     ',/,
     -'         d12=i1*g(i1-1,i1,i1,i4,i4,i4)+                      ',/,
     -'     ,       i1*g(i1-1,i4,i4,i1,i1,i4)+                      ',/,
     -'     ,       i4*g(i4-1,i1,i4,i1,i4,i1)+                      ',/,
     -'     ,       i4*g(i4-1,i4,i1,i4,i1,i1)                       ',/,
     -'         d13=i1*g(i1,i1-1,i1,i4,i4,i4)+                      ',/,
     -'     ,       i4*g(i1,i4-1,i4,i1,i1,i4)+                      ',/,
     -'     ,       i1*g(i4,i1-1,i4,i1,i4,i1)+                      ',/,
     -'     ,       i4*g(i4,i4-1,i1,i4,i1,i1)                       ',/,
     -'         d14=i1*g(i1,i1,i1-1,i4,i4,i4)+                      ',/,
     -'     ,       i4*g(i1,i4,i4-1,i1,i1,i4)+                      ',/,
     -'     ,       i4*g(i4,i1,i4-1,i1,i4,i1)+                      ',/,
     -'     ,       i1*g(i4,i4,i1-1,i4,i1,i1)                       ',/,
     -'         d23=i4*g(i1,i1,i1,i4-1,i4,i4)+                      ',/,
     -'     ,       i1*g(i1,i4,i4,i1-1,i1,i4)+                      ',/,
     -'     ,       i1*g(i4,i1,i4,i1-1,i4,i1)+                      ',/,
     -'     ,       i4*g(i4,i4,i1,i4-1,i1,i1)                       ',/,
     -'         d24=i4*g(i1,i1,i1,i4,i4-1,i4)+                      ',/,
     -'     ,       i1*g(i1,i4,i4,i1,i1-1,i4)+                      ',/,
     -'     ,       i4*g(i4,i1,i4,i1,i4-1,i1)+                      ',/,
     -'     ,       i1*g(i4,i4,i1,i4,i1-1,i1)                       ',/,
     -'         d34=i4*g(i1,i1,i1,i4,i4,i4-1)+                      ',/,
     -'     ,       i4*g(i1,i4,i4,i1,i1,i4-1)+                      ',/,
     -'     ,       i1*g(i4,i1,i4,i1,i4,i1-1)+                      ',/,
     -'     ,       i1*g(i4,i4,i1,i4,i1,i1-1)                       ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.6) then                                  ',/,
     -'!        faaaa=g(i1,i2,i3,i4,i5,i6)+g(i1,i3,i2,i5,i4,i6)+    ',/,
     -'!    ,         g(i3,i1,i2,i5,i6,i4)+g(i4,i5,i6,i1,i2,i3)+    ',/,
     -'!    ,         g(i5,i4,i6,i1,i3,i2)+g(i5,i6,i4,i3,i1,i2)     ',/,
     -'         d12=i1*g(i1-1,i2,i3,i4,i5,i6)+                      ',/,
     -'     ,       i1*g(i1-1,i3,i2,i5,i4,i6)+                      ',/,
     -'     ,       i3*g(i3-1,i1,i2,i5,i6,i4)+                      ',/,
     -'     ,       i4*g(i4-1,i5,i6,i1,i2,i3)+                      ',/,
     -'     ,       i5*g(i5-1,i4,i6,i1,i3,i2)+                      ',/,
     -'     ,       i5*g(i5-1,i6,i4,i3,i1,i2)                       ',/,
     -'         d13=i2*g(i1,i2-1,i3,i4,i5,i6)+                      ',/,
     -'     ,       i3*g(i1,i3-1,i2,i5,i4,i6)+                      ',/,
     -'     ,       i1*g(i3,i1-1,i2,i5,i6,i4)+                      ',/,
     -'     ,       i5*g(i4,i5-1,i6,i1,i2,i3)+                      ',/,
     -'     ,       i4*g(i5,i4-1,i6,i1,i3,i2)+                      ',/,
     -'     ,       i6*g(i5,i6-1,i4,i3,i1,i2)                       ',/,
     -'         d14=i3*g(i1,i2,i3-1,i4,i5,i6)+                      ',/,
     -'     ,       i2*g(i1,i3,i2-1,i5,i4,i6)+                      ',/,
     -'     ,       i2*g(i3,i1,i2-1,i5,i6,i4)+                      ',/,
     -'     ,       i6*g(i4,i5,i6-1,i1,i2,i3)+                      ',/,
     -'     ,       i6*g(i5,i4,i6-1,i1,i3,i2)+                      ',/,
     -'     ,       i4*g(i5,i6,i4-1,i3,i1,i2)                       ',/,
     -'         d23=i4*g(i1,i2,i3,i4-1,i5,i6)+                      ',/,
     -'     ,       i5*g(i1,i3,i2,i5-1,i4,i6)+                      ',/,
     -'     ,       i5*g(i3,i1,i2,i5-1,i6,i4)+                      ',/,
     -'     ,       i1*g(i4,i5,i6,i1-1,i2,i3)+                      ',/,
     -'     ,       i1*g(i5,i4,i6,i1-1,i3,i2)+                      ',/,
     -'     ,       i3*g(i5,i6,i4,i3-1,i1,i2)                       ',/,
     -'         d24=i5*g(i1,i2,i3,i4,i5-1,i6)+                      ',/,
     -'     ,       i4*g(i1,i3,i2,i5,i4-1,i6)+                      ',/,
     -'     ,       i6*g(i3,i1,i2,i5,i6-1,i4)+                      ',/,
     -'     ,       i2*g(i4,i5,i6,i1,i2-1,i3)+                      ',/,
     -'     ,       i3*g(i5,i4,i6,i1,i3-1,i2)+                      ',/,
     -'     ,       i1*g(i5,i6,i4,i3,i1-1,i2)                       ',/,
     -'         d34=i6*g(i1,i2,i3,i4,i5,i6-1)+                      ',/,
     -'     ,       i6*g(i1,i3,i2,i5,i4,i6-1)+                      ',/,
     -'     ,       i4*g(i3,i1,i2,i5,i6,i4-1)+                      ',/,
     -'     ,       i3*g(i4,i5,i6,i1,i2,i3-1)+                      ',/,
     -'     ,       i2*g(i5,i4,i6,i1,i3,i2-1)+                      ',/,
     -'     ,       i2*g(i5,i6,i4,i3,i1,i2-1)                       ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.11) then                                 ',/,
     -'!        faaaa=g(i1,i2,i3,i4,i2,i1)+g(i1,i2,i4,i3,i2,i1)+    ',/,
     -'!    ,         g(i1,i3,i2,i2,i4,i1)+g(i1,i4,i2,i2,i3,i1)+    ',/,
     -'!    ,         g(i2,i1,i3,i4,i1,i2)+g(i2,i1,i4,i3,i1,i2)+    ',/,
     -'!    ,         g(i2,i3,i1,i1,i4,i2)+g(i2,i4,i1,i1,i3,i2)+    ',/,
     -'!    ,         g(i3,i1,i2,i2,i1,i4)+g(i3,i2,i1,i1,i2,i4)+    ',/,
     -'!    ,         g(i4,i1,i2,i2,i1,i3)+g(i4,i2,i1,i1,i2,i3)     ',/,
     -'         d12=i1*g(i1-1,i2,i3,i4,i2,i1)+                      ',/,
     -'     ,       i1*g(i1-1,i2,i4,i3,i2,i1)+                      ',/,
     -'     ,       i1*g(i1-1,i3,i2,i2,i4,i1)+                      ',/,
     -'     ,       i1*g(i1-1,i4,i2,i2,i3,i1)+                      ',/,
     -'     ,       i2*g(i2-1,i1,i3,i4,i1,i2)+                      ',/,
     -'     ,       i2*g(i2-1,i1,i4,i3,i1,i2)+                      ',/,
     -'     ,       i2*g(i2-1,i3,i1,i1,i4,i2)+                      ',/,
     -'     ,       i2*g(i2-1,i4,i1,i1,i3,i2)+                      ',/,
     -'     ,       i3*g(i3-1,i1,i2,i2,i1,i4)+                      ',/,
     -'     ,       i3*g(i3-1,i2,i1,i1,i2,i4)+                      ',/,
     -'     ,       i4*g(i4-1,i1,i2,i2,i1,i3)+                      ',/,
     -'     ,       i4*g(i4-1,i2,i1,i1,i2,i3)                       ',/,
     -'         d13=i2*g(i1,i2-1,i3,i4,i2,i1)+                      ',/,
     -'     ,       i2*g(i1,i2-1,i4,i3,i2,i1)+                      ',/,
     -'     ,       i3*g(i1,i3-1,i2,i2,i4,i1)+                      ',/,
     -'     ,       i4*g(i1,i4-1,i2,i2,i3,i1)+                      ',/,
     -'     ,       i1*g(i2,i1-1,i3,i4,i1,i2)+                      ',/,
     -'     ,       i1*g(i2,i1-1,i4,i3,i1,i2)+                      ',/,
     -'     ,       i3*g(i2,i3-1,i1,i1,i4,i2)+                      ',/,
     -'     ,       i4*g(i2,i4-1,i1,i1,i3,i2)+                      ',/,
     -'     ,       i1*g(i3,i1-1,i2,i2,i1,i4)+                      ',/,
     -'     ,       i2*g(i3,i2-1,i1,i1,i2,i4)+                      ',/,
     -'     ,       i1*g(i4,i1-1,i2,i2,i1,i3)+                      ',/,
     -'     ,       i2*g(i4,i2-1,i1,i1,i2,i3)                       ',/,
     -'         d14=i3*g(i1,i2,i3-1,i4,i2,i1)+                      ',/,
     -'     ,       i4*g(i1,i2,i4-1,i3,i2,i1)+                      ',/,
     -'     ,       i2*g(i1,i3,i2-1,i2,i4,i1)+                      ',/,
     -'     ,       i2*g(i1,i4,i2-1,i2,i3,i1)+                      ',/,
     -'     ,       i3*g(i2,i1,i3-1,i4,i1,i2)+                      ',/,
     -'     ,       i4*g(i2,i1,i4-1,i3,i1,i2)+                      ',/,
     -'     ,       i1*g(i2,i3,i1-1,i1,i4,i2)+                      ',/,
     -'     ,       i1*g(i2,i4,i1-1,i1,i3,i2)+                      ',/,
     -'     ,       i2*g(i3,i1,i2-1,i2,i1,i4)+                      ',/,
     -'     ,       i1*g(i3,i2,i1-1,i1,i2,i4)+                      ',/,
     -'     ,       i2*g(i4,i1,i2-1,i2,i1,i3)+                      ',/,
     -'     ,       i1*g(i4,i2,i1-1,i1,i2,i3)                       ',/,
     -'         d23=i4*g(i1,i2,i3,i4-1,i2,i1)+                      ',/,
     -'     ,       i3*g(i1,i2,i4,i3-1,i2,i1)+                      ',/,
     -'     ,       i2*g(i1,i3,i2,i2-1,i4,i1)+                      ',/,
     -'     ,       i2*g(i1,i4,i2,i2-1,i3,i1)+                      ',/,
     -'     ,       i4*g(i2,i1,i3,i4-1,i1,i2)+                      ',/,
     -'     ,       i3*g(i2,i1,i4,i3-1,i1,i2)+                      ',/,
     -'     ,       i1*g(i2,i3,i1,i1-1,i4,i2)+                      ',/,
     -'     ,       i1*g(i2,i4,i1,i1-1,i3,i2)+                      ',/,
     -'     ,       i2*g(i3,i1,i2,i2-1,i1,i4)+                      ',/,
     -'     ,       i1*g(i3,i2,i1,i1-1,i2,i4)+                      ',/,
     -'     ,       i2*g(i4,i1,i2,i2-1,i1,i3)+                      ',/,
     -'     ,       i1*g(i4,i2,i1,i1-1,i2,i3)                       ',/,
     -'         d24=i2*g(i1,i2,i3,i4,i2-1,i1)+                      ',/,
     -'     ,       i2*g(i1,i2,i4,i3,i2-1,i1)+                      ',/,
     -'     ,       i4*g(i1,i3,i2,i2,i4-1,i1)+                      ',/,
     -'     ,       i3*g(i1,i4,i2,i2,i3-1,i1)+                      ',/,
     -'     ,       i1*g(i2,i1,i3,i4,i1-1,i2)+                      ',/,
     -'     ,       i1*g(i2,i1,i4,i3,i1-1,i2)+                      ',/,
     -'     ,       i4*g(i2,i3,i1,i1,i4-1,i2)+                      ',/,
     -'     ,       i3*g(i2,i4,i1,i1,i3-1,i2)+                      ',/,
     -'     ,       i1*g(i3,i1,i2,i2,i1-1,i4)+                      ',/,
     -'     ,       i2*g(i3,i2,i1,i1,i2-1,i4)+                      ',/,
     -'     ,       i1*g(i4,i1,i2,i2,i1-1,i3)+                      ',/,
     -'     ,       i2*g(i4,i2,i1,i1,i2-1,i3)                       ',/,
     -'         d34=i1*g(i1,i2,i3,i4,i2,i1-1)+                      ',/,
     -'     ,       i1*g(i1,i2,i4,i3,i2,i1-1)+                      ',/,
     -'     ,       i1*g(i1,i3,i2,i2,i4,i1-1)+                      ',/,
     -'     ,       i1*g(i1,i4,i2,i2,i3,i1-1)+                      ',/,
     -'     ,       i2*g(i2,i1,i3,i4,i1,i2-1)+                      ',/,
     -'     ,       i2*g(i2,i1,i4,i3,i1,i2-1)+                      ',/,
     -'     ,       i2*g(i2,i3,i1,i1,i4,i2-1)+                      ',/,
     -'     ,       i2*g(i2,i4,i1,i1,i3,i2-1)+                      ',/,
     -'     ,       i4*g(i3,i1,i2,i2,i1,i4-1)+                      ',/,
     -'     ,       i4*g(i3,i2,i1,i1,i2,i4-1)+                      ',/,
     -'     ,       i3*g(i4,i1,i2,i2,i1,i3-1)+                      ',/,
     -'     ,       i3*g(i4,i2,i1,i1,i2,i3-1)                       ',/,
     -'         return                                              ',/,
     -'      elseif (i7.eq.13) then                                 ',/,
     -'!        faaaa=g(i1,i1,i3,i4,i6,i6)+g(i1,i6,i4,i3,i1,i6)+    ',/,
     -'!    ,         g(i1,i3,i1,i6,i4,i6)+g(i1,i4,i6,i1,i3,i6)+    ',/,
     -'!    ,         g(i6,i1,i4,i3,i6,i1)+g(i6,i6,i3,i4,i1,i1)+    ',/,
     -'!    ,         g(i6,i3,i6,i1,i4,i1)+g(i6,i4,i1,i6,i3,i1)+    ',/,
     -'!    ,         g(i3,i1,i1,i6,i6,i4)+g(i3,i6,i6,i1,i1,i4)+    ',/,
     -'!    ,         g(i4,i1,i6,i1,i6,i3)+g(i4,i6,i1,i6,i1,i3)     ',/,
     -'         d12=i1*g(i1-1,i1,i3,i4,i6,i6)+                      ',/,
     -'     ,       i1*g(i1-1,i6,i4,i3,i1,i6)+                      ',/,
     -'     ,       i1*g(i1-1,i3,i1,i6,i4,i6)+                      ',/,
     -'     ,       i1*g(i1-1,i4,i6,i1,i3,i6)+                      ',/,
     -'     ,       i6*g(i6-1,i1,i4,i3,i6,i1)+                      ',/,
     -'     ,       i6*g(i6-1,i6,i3,i4,i1,i1)+                      ',/,
     -'     ,       i6*g(i6-1,i3,i6,i1,i4,i1)+                      ',/,
     -'     ,       i6*g(i6-1,i4,i1,i6,i3,i1)+                      ',/,
     -'     ,       i3*g(i3-1,i1,i1,i6,i6,i4)+                      ',/,
     -'     ,       i3*g(i3-1,i6,i6,i1,i1,i4)+                      ',/,
     -'     ,       i4*g(i4-1,i1,i6,i1,i6,i3)+                      ',/,
     -'     ,       i4*g(i4-1,i6,i1,i6,i1,i3)                       ',/,
     -'         d13=i1*g(i1,i1-1,i3,i4,i6,i6)+                      ',/,
     -'     ,       i6*g(i1,i6-1,i4,i3,i1,i6)+                      ',/,
     -'     ,       i3*g(i1,i3-1,i1,i6,i4,i6)+                      ',/,
     -'     ,       i4*g(i1,i4-1,i6,i1,i3,i6)+                      ',/,
     -'     ,       i1*g(i6,i1-1,i4,i3,i6,i1)+                      ',/,
     -'     ,       i6*g(i6,i6-1,i3,i4,i1,i1)+                      ',/,
     -'     ,       i3*g(i6,i3-1,i6,i1,i4,i1)+                      ',/,
     -'     ,       i4*g(i6,i4-1,i1,i6,i3,i1)+                      ',/,
     -'     ,       i1*g(i3,i1-1,i1,i6,i6,i4)+                      ',/,
     -'     ,       i6*g(i3,i6-1,i6,i1,i1,i4)+                      ',/,
     -'     ,       i1*g(i4,i1-1,i6,i1,i6,i3)+                      ',/,
     -'     ,       i6*g(i4,i6-1,i1,i6,i1,i3)                       ',/,
     -'         d14=i3*g(i1,i1,i3-1,i4,i6,i6)+                      ',/,
     -'     ,       i4*g(i1,i6,i4-1,i3,i1,i6)+                      ',/,
     -'     ,       i1*g(i1,i3,i1-1,i6,i4,i6)+                      ',/,
     -'     ,       i6*g(i1,i4,i6-1,i1,i3,i6)+                      ',/,
     -'     ,       i4*g(i6,i1,i4-1,i3,i6,i1)+                      ',/,
     -'     ,       i3*g(i6,i6,i3-1,i4,i1,i1)+                      ',/,
     -'     ,       i6*g(i6,i3,i6-1,i1,i4,i1)+                      ',/,
     -'     ,       i1*g(i6,i4,i1-1,i6,i3,i1)+                      ',/,
     -'     ,       i1*g(i3,i1,i1-1,i6,i6,i4)+                      ',/,
     -'     ,       i6*g(i3,i6,i6-1,i1,i1,i4)+                      ',/,
     -'     ,       i6*g(i4,i1,i6-1,i1,i6,i3)+                      ',/,
     -'     ,       i1*g(i4,i6,i1-1,i6,i1,i3)                       ',/,
     -'         d23=i4*g(i1,i1,i3,i4-1,i6,i6)+                      ',/,
     -'     ,       i3*g(i1,i6,i4,i3-1,i1,i6)+                      ',/,
     -'     ,       i6*g(i1,i3,i1,i6-1,i4,i6)+                      ',/,
     -'     ,       i1*g(i1,i4,i6,i1-1,i3,i6)+                      ',/,
     -'     ,       i3*g(i6,i1,i4,i3-1,i6,i1)+                      ',/,
     -'     ,       i4*g(i6,i6,i3,i4-1,i1,i1)+                      ',/,
     -'     ,       i1*g(i6,i3,i6,i1-1,i4,i1)+                      ',/,
     -'     ,       i6*g(i6,i4,i1,i6-1,i3,i1)+                      ',/,
     -'     ,       i6*g(i3,i1,i1,i6-1,i6,i4)+                      ',/,
     -'     ,       i1*g(i3,i6,i6,i1-1,i1,i4)+                      ',/,
     -'     ,       i1*g(i4,i1,i6,i1-1,i6,i3)+                      ',/,
     -'     ,       i6*g(i4,i6,i1,i6-1,i1,i3)                       ',/,
     -'         d24=i6*g(i1,i1,i3,i4,i6-1,i6)+                      ',/,
     -'     ,       i1*g(i1,i6,i4,i3,i1-1,i6)+                      ',/,
     -'     ,       i4*g(i1,i3,i1,i6,i4-1,i6)+                      ',/,
     -'     ,       i3*g(i1,i4,i6,i1,i3-1,i6)+                      ',/,
     -'     ,       i6*g(i6,i1,i4,i3,i6-1,i1)+                      ',/,
     -'     ,       i1*g(i6,i6,i3,i4,i1-1,i1)+                      ',/,
     -'     ,       i4*g(i6,i3,i6,i1,i4-1,i1)+                      ',/,
     -'     ,       i3*g(i6,i4,i1,i6,i3-1,i1)+                      ',/,
     -'     ,       i6*g(i3,i1,i1,i6,i6-1,i4)+                      ',/,
     -'     ,       i1*g(i3,i6,i6,i1,i1-1,i4)+                      ',/,
     -'     ,       i6*g(i4,i1,i6,i1,i6-1,i3)+                      ',/,
     -'     ,       i1*g(i4,i6,i1,i6,i1-1,i3)                       ',/,
     -'         d34=i6*g(i1,i1,i3,i4,i6,i6-1)+                      ',/,
     -'     ,       i6*g(i1,i6,i4,i3,i1,i6-1)+                      ',/,
     -'     ,       i6*g(i1,i3,i1,i6,i4,i6-1)+                      ',/,
     -'     ,       i6*g(i1,i4,i6,i1,i3,i6-1)+                      ',/,
     -'     ,       i1*g(i6,i1,i4,i3,i6,i1-1)+                      ',/,
     -'     ,       i1*g(i6,i6,i3,i4,i1,i1-1)+                      ',/,
     -'     ,       i1*g(i6,i3,i6,i1,i4,i1-1)+                      ',/,
     -'     ,       i1*g(i6,i4,i1,i6,i3,i1-1)+                      ',/,
     -'     ,       i4*g(i3,i1,i1,i6,i6,i4-1)+                      ',/,
     -'     ,       i4*g(i3,i6,i6,i1,i1,i4-1)+                      ',/,
     -'     ,       i3*g(i4,i1,i6,i1,i6,i3-1)+                      ',/,
     -'     ,       i3*g(i4,i6,i1,i6,i1,i3-1)                       ',/,
     -'         return                                              ',/,
     -'      else                                                   ',/,
     -'!         faaaa=                                             ',/,
     -'!    *  g(i6,i5,i3,i4,i2,i1)+g(i6,i4,i2,i5,i3,i1)+           ',/,
     -'!    *  g(i6,i3,i5,i2,i4,i1)+                                ',/,
     -'!    *  g(i6,i2,i4,i3,i5,i1)+g(i5,i6,i3,i4,i1,i2)+           ',/,
     -'!    *  g(i5,i4,i1,i6,i3,i2)+                                ',/,
     -'!    *  g(i5,i3,i6,i1,i4,i2)+g(i5,i1,i4,i3,i6,i2)+           ',/,
     -'!    *  g(i4,i6,i2,i5,i1,i3)+                                ',/,      
     -'!    *  g(i4,i5,i1,i6,i2,i3)+g(i4,i2,i6,i1,i5,i3)+           ',/,
     -'!    *  g(i4,i1,i5,i2,i6,i3)+                                ',/,      
     -'!    *  g(i3,i6,i5,i2,i1,i4)+g(i3,i5,i6,i1,i2,i4)+           ',/,
     -'!    *  g(i3,i2,i1,i6,i5,i4)+                                ',/,      
     -'!    *  g(i3,i1,i2,i5,i6,i4)+g(i2,i6,i4,i3,i1,i5)+           ',/,
     -'!    *  g(i2,i4,i6,i1,i3,i5)+                                ',/,      
     -'!    *  g(i2,i3,i1,i6,i4,i5)+g(i2,i1,i3,i4,i6,i5)+           ',/,
     -'!    *  g(i1,i5,i4,i3,i2,i6)+                                ',/,      
     -'!    *  g(i1,i4,i5,i2,i3,i6)+g(i1,i3,i2,i5,i4,i6)+           ',/,
     -'!    *  g(i1,i2,i3,i4,i5,i6)                                 ',/,      
     -'          d12=                                               ',/,
     -'     *   i6*g(i6-1,i5,i3,i4,i2,i1)+i6*g(i6-1,i4,i2,i5,i3,i1)+',/,
     -'     *   i6*g(i6-1,i3,i5,i2,i4,i1)+i6*g(i6-1,i2,i4,i3,i5,i1)+',/,
     -'     *   i5*g(i5-1,i6,i3,i4,i1,i2)+i5*g(i5-1,i4,i1,i6,i3,i2)+',/,
     -'     *   i5*g(i5-1,i3,i6,i1,i4,i2)+i5*g(i5-1,i1,i4,i3,i6,i2)+',/,
     -'     *   i4*g(i4-1,i6,i2,i5,i1,i3)+i4*g(i4-1,i5,i1,i6,i2,i3)+',/,
     -'     *   i4*g(i4-1,i2,i6,i1,i5,i3)+i4*g(i4-1,i1,i5,i2,i6,i3)+',/,
     -'     *   i3*g(i3-1,i6,i5,i2,i1,i4)+i3*g(i3-1,i5,i6,i1,i2,i4)+',/,
     -'     *   i3*g(i3-1,i2,i1,i6,i5,i4)+i3*g(i3-1,i1,i2,i5,i6,i4)+',/,
     -'     *   i2*g(i2-1,i6,i4,i3,i1,i5)+i2*g(i2-1,i4,i6,i1,i3,i5)+',/,
     -'     *   i2*g(i2-1,i3,i1,i6,i4,i5)+i2*g(i2-1,i1,i3,i4,i6,i5)+',/,
     -'     *   i1*g(i1-1,i5,i4,i3,i2,i6)+i1*g(i1-1,i4,i5,i2,i3,i6)+',/,
     -'     *   i1*g(i1-1,i3,i2,i5,i4,i6)+i1*g(i1-1,i2,i3,i4,i5,i6) ',/,
     -'          d13=                                               ',/,
     -'     *   i5*g(i6,i5-1,i3,i4,i2,i1)+i4*g(i6,i4-1,i2,i5,i3,i1)+',/,
     -'     *   i3*g(i6,i3-1,i5,i2,i4,i1)+i2*g(i6,i2-1,i4,i3,i5,i1)+',/,
     -'     *   i6*g(i5,i6-1,i3,i4,i1,i2)+i4*g(i5,i4-1,i1,i6,i3,i2)+',/,
     -'     *   i3*g(i5,i3-1,i6,i1,i4,i2)+i1*g(i5,i1-1,i4,i3,i6,i2)+',/,
     -'     *   i6*g(i4,i6-1,i2,i5,i1,i3)+i5*g(i4,i5-1,i1,i6,i2,i3)+',/,
     -'     *   i2*g(i4,i2-1,i6,i1,i5,i3)+i1*g(i4,i1-1,i5,i2,i6,i3)+',/,
     -'     *   i6*g(i3,i6-1,i5,i2,i1,i4)+i5*g(i3,i5-1,i6,i1,i2,i4)+',/,
     -'     *   i2*g(i3,i2-1,i1,i6,i5,i4)+i1*g(i3,i1-1,i2,i5,i6,i4)+',/,
     -'     *   i6*g(i2,i6-1,i4,i3,i1,i5)+i4*g(i2,i4-1,i6,i1,i3,i5)+',/,
     -'     *   i3*g(i2,i3-1,i1,i6,i4,i5)+i1*g(i2,i1-1,i3,i4,i6,i5)+',/,
     -'     *   i5*g(i1,i5-1,i4,i3,i2,i6)+i4*g(i1,i4-1,i5,i2,i3,i6)+',/,
     -'     *   i3*g(i1,i3-1,i2,i5,i4,i6)+i2*g(i1,i2-1,i3,i4,i5,i6) ',/,
     -'          d14=                                               ',/,
     -'     *   i3*g(i6,i5,i3-1,i4,i2,i1)+i2*g(i6,i4,i2-1,i5,i3,i1)+',/,
     -'     *   i5*g(i6,i3,i5-1,i2,i4,i1)+i4*g(i6,i2,i4-1,i3,i5,i1)+',/,
     -'     *   i3*g(i5,i6,i3-1,i4,i1,i2)+i1*g(i5,i4,i1-1,i6,i3,i2)+',/,
     -'     *   i6*g(i5,i3,i6-1,i1,i4,i2)+i4*g(i5,i1,i4-1,i3,i6,i2)+',/,
     -'     *   i2*g(i4,i6,i2-1,i5,i1,i3)+i1*g(i4,i5,i1-1,i6,i2,i3)+',/,
     -'     *   i6*g(i4,i2,i6-1,i1,i5,i3)+i5*g(i4,i1,i5-1,i2,i6,i3)+',/,
     -'     *   i5*g(i3,i6,i5-1,i2,i1,i4)+i6*g(i3,i5,i6-1,i1,i2,i4)+',/,
     -'     *   i1*g(i3,i2,i1-1,i6,i5,i4)+i2*g(i3,i1,i2-1,i5,i6,i4)+',/,
     -'     *   i4*g(i2,i6,i4-1,i3,i1,i5)+i6*g(i2,i4,i6-1,i1,i3,i5)+',/,
     -'     *   i1*g(i2,i3,i1-1,i6,i4,i5)+i3*g(i2,i1,i3-1,i4,i6,i5)+',/,
     -'     *   i4*g(i1,i5,i4-1,i3,i2,i6)+i5*g(i1,i4,i5-1,i2,i3,i6)+',/,
     -'     *   i2*g(i1,i3,i2-1,i5,i4,i6)+i3*g(i1,i2,i3-1,i4,i5,i6) ',/,
     -'          d23=                                               ',/,
     -'     *   i4*g(i6,i5,i3,i4-1,i2,i1)+i5*g(i6,i4,i2,i5-1,i3,i1)+',/,
     -'     *   i2*g(i6,i3,i5,i2-1,i4,i1)+i3*g(i6,i2,i4,i3-1,i5,i1)+',/,
     -'     *   i4*g(i5,i6,i3,i4-1,i1,i2)+i6*g(i5,i4,i1,i6-1,i3,i2)+',/,
     -'     *   i1*g(i5,i3,i6,i1-1,i4,i2)+i3*g(i5,i1,i4,i3-1,i6,i2)+',/,
     -'     *   i5*g(i4,i6,i2,i5-1,i1,i3)+i6*g(i4,i5,i1,i6-1,i2,i3)+',/,
     -'     *   i1*g(i4,i2,i6,i1-1,i5,i3)+i2*g(i4,i1,i5,i2-1,i6,i3)+',/,
     -'     *   i2*g(i3,i6,i5,i2-1,i1,i4)+i1*g(i3,i5,i6,i1-1,i2,i4)+',/,
     -'     *   i6*g(i3,i2,i1,i6-1,i5,i4)+i5*g(i3,i1,i2,i5-1,i6,i4)+',/,
     -'     *   i3*g(i2,i6,i4,i3-1,i1,i5)+i1*g(i2,i4,i6,i1-1,i3,i5)+',/,
     -'     *   i6*g(i2,i3,i1,i6-1,i4,i5)+i4*g(i2,i1,i3,i4-1,i6,i5)+',/,
     -'     *   i3*g(i1,i5,i4,i3-1,i2,i6)+i2*g(i1,i4,i5,i2-1,i3,i6)+',/,
     -'     *   i5*g(i1,i3,i2,i5-1,i4,i6)+i4*g(i1,i2,i3,i4-1,i5,i6) ',/,
     -'          d24=                                               ',/,
     -'     *   i2*g(i6,i5,i3,i4,i2-1,i1)+i3*g(i6,i4,i2,i5,i3-1,i1)+',/,
     -'     *   i4*g(i6,i3,i5,i2,i4-1,i1)+i5*g(i6,i2,i4,i3,i5-1,i1)+',/,
     -'     *   i1*g(i5,i6,i3,i4,i1-1,i2)+i3*g(i5,i4,i1,i6,i3-1,i2)+',/,
     -'     *   i4*g(i5,i3,i6,i1,i4-1,i2)+i6*g(i5,i1,i4,i3,i6-1,i2)+',/,
     -'     *   i1*g(i4,i6,i2,i5,i1-1,i3)+i2*g(i4,i5,i1,i6,i2-1,i3)+',/,
     -'     *   i5*g(i4,i2,i6,i1,i5-1,i3)+i6*g(i4,i1,i5,i2,i6-1,i3)+',/,
     -'     *   i1*g(i3,i6,i5,i2,i1-1,i4)+i2*g(i3,i5,i6,i1,i2-1,i4)+',/,
     -'     *   i5*g(i3,i2,i1,i6,i5-1,i4)+i6*g(i3,i1,i2,i5,i6-1,i4)+',/,
     -'     *   i1*g(i2,i6,i4,i3,i1-1,i5)+i3*g(i2,i4,i6,i1,i3-1,i5)+',/,
     -'     *   i4*g(i2,i3,i1,i6,i4-1,i5)+i6*g(i2,i1,i3,i4,i6-1,i5)+',/,
     -'     *   i2*g(i1,i5,i4,i3,i2-1,i6)+i3*g(i1,i4,i5,i2,i3-1,i6)+',/,
     -'     *   i4*g(i1,i3,i2,i5,i4-1,i6)+i5*g(i1,i2,i3,i4,i5-1,i6) ',/,
     -'          d34=                                               ',/,
     -'     *   i1*g(i6,i5,i3,i4,i2,i1-1)+i1*g(i6,i4,i2,i5,i3,i1-1)+',/,
     -'     *   i1*g(i6,i3,i5,i2,i4,i1-1)+i1*g(i6,i2,i4,i3,i5,i1-1)+',/,
     -'     *   i2*g(i5,i6,i3,i4,i1,i2-1)+i2*g(i5,i4,i1,i6,i3,i2-1)+',/,
     -'     *   i2*g(i5,i3,i6,i1,i4,i2-1)+i2*g(i5,i1,i4,i3,i6,i2-1)+',/,
     -'     *   i3*g(i4,i6,i2,i5,i1,i3-1)+i3*g(i4,i5,i1,i6,i2,i3-1)+',/,
     -'     *   i3*g(i4,i2,i6,i1,i5,i3-1)+i3*g(i4,i1,i5,i2,i6,i3-1)+',/,
     -'     *   i4*g(i3,i6,i5,i2,i1,i4-1)+i4*g(i3,i5,i6,i1,i2,i4-1)+',/,
     -'     *   i4*g(i3,i2,i1,i6,i5,i4-1)+i4*g(i3,i1,i2,i5,i6,i4-1)+',/,
     -'     *   i5*g(i2,i6,i4,i3,i1,i5-1)+i5*g(i2,i4,i6,i1,i3,i5-1)+',/,
     -'     *   i5*g(i2,i3,i1,i6,i4,i5-1)+i5*g(i2,i1,i3,i4,i6,i5-1)+',/,
     -'     *   i6*g(i1,i5,i4,i3,i2,i6-1)+i6*g(i1,i4,i5,i2,i3,i6-1)+',/,
     -'     *   i6*g(i1,i3,i2,i5,i4,i6-1)+i6*g(i1,i2,i3,i4,i5,i6-1) ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      return                                                 ',/,
     -'      end                                                    ')
      return
      end
**********************************************************************
      subroutine gfabbb(n,mt,vex2)
**********************************************************************
*     devuelve el valor de la fas j para cada punto i en xx(j,i)
**********************************************************************
      implicit real*8(a-h,o-z)
      parameter (nmax=6200,nexp=865)
      parameter (igrad4=12,nexp42=756)
      common/ipot42/ipa42(nexp42,7),nfas42(igrad4)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      dimension vex2(2)
      do 1 i=1,n
         f12(0) = 1.d0
         f13(0) = 1.d0
         f14(0) = 1.d0
         f23(0) = 1.d0
         f24(0) = 1.d0
         f34(0) = 1.d0
         aux12= r12(i)*dexp(-vex2(1)*r12(i))
         aux13= r13(i)*dexp(-vex2(1)*r13(i))
         aux14= r14(i)*dexp(-vex2(1)*r14(i))
         aux23= r23(i)*dexp(-vex2(2)*r23(i))
         aux24= r24(i)*dexp(-vex2(2)*r24(i))
         aux34= r34(i)*dexp(-vex2(2)*r34(i))
         do 2 k=1,mt-2
            f12(k) = aux12*f12(k-1)
            f13(k) = aux13*f13(k-1)
            f14(k) = aux14*f14(k-1)
            f23(k) = aux23*f23(k-1)
            f24(k) = aux24*f24(k-1)
            f34(k) = aux34*f34(k-1)
2        continue
         do 3 j=1,nfas42(mt)
            i1=ipa42(j,1)
            i2=ipa42(j,2)
            i3=ipa42(j,3)
            i4=ipa42(j,4)
            i5=ipa42(j,5)
            i6=ipa42(j,6)
            i7=ipa42(j,7)
            xx(j,i)=fabbb(i1,i2,i3,i4,i5,i6,i7)
 3       continue
 1    continue
      end
**********************************************************************
      function fabbb(k1,k2,k3,k4,k5,k6,jt)
      implicit real*8(a-h,o-z)
      parameter (igrad4=12)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
      if (jt.eq.1) then
         fabbb=g(k1,k2,k3,k4,k5,k6)
         return
      endif
      if (jt.eq.2) then
         fabbb=(g(k1,k1,k3,k4,k5,k5)+g(k1,k3,k1,k5,k4,k5)+
     *   g(k3,k1,k1,k5,k5,k4))
         return
      endif
      if (jt.eq.3) then
         fabbb=g(k1,k2,k3,k4,k5,k6)+g(k2,k3,k1,k6,k4,k5)+
     *         g(k3,k1,k2,k5,k6,k4)+g(k3,k2,k1,k6,k5,k4)+
     *         g(k1,k3,k2,k5,k4,k6)+g(k2,k1,k3,k4,k6,k5)
         return
      endif
      end
****************************************************************
      subroutine sabbb(nombre,n,mt,vex2,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie ABBB
c     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*4 nombre
      include 'dimensions.inc'
      common/ipot42/ipa42(nexp42,7),nfas42(igrad4)
      common sm(nexp)
      dimension vex2(2)
      is=nfas42(mt)
      write(7,1) nombre,n,rms*627.51d0,emax*627.51d0,
     -is,mt-2
      do 2 k=1,is
         write(7,3) k,sm(k),k,ipa42(k,1),k,ipa42(k,2),k,ipa42(k,3),
     -         k,ipa42(k,4),k,ipa42(k,5),k,ipa42(k,6),k,ipa42(k,7)
2     continue
      do i=1,2
         write(7,4) i,vex2(i)
4        format(6x,'vex2(',i1,')=',d23.15)
      enddo
      write(7,10)
      write(7,11)
      write(7,12)
      write(7,13) mt-2
1     format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 6D PES      ',/,
     -'*     for the ABBB system class fitted to ',i4,' points      ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real * 8 (a-h,o-z)                            ',/,
     -'      parameter(is=',i5,', mt_2=',i2,')                      ',/,
     -'      dimension vex2(2)                                      ',/,
     -'      dimension der(6)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),',
     ,'i7(is),cf(is)',/,
     -'      common/basu4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),   ',/, 
     -'     ,             f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)    ',/, 
     -'      data der12/0.d0/,der13/0.d0/,der14/0.d0/    ',/,
     -'      data der23/0.d0/,der24/0.d0/,der34/0.d0/    ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/      ',/,
     -'      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/      ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/         ',/,
     -'      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            ')
3     format(
     -'      data cf(',i5,')/',d23.15,'/                           ',/,
     -'      data i1(',i5,')/',i2,'/,i2(',i5,')/',i2,'/,i3(',i5,')/',i2,
     -'/,i4(',i5,')/',i2,'/                                        ',/,
     -'      data i5(',i5,')/',i2,'/,i6(',i5,')/',i2,'/,',
     -'i7(',i5,')/',i2,'/')
10    format(
     -'      ener=0.d0          ',/,
     -'      der12=0.d0     ',/,
     -'      der13=0.d0     ',/,
     -'      der14=0.d0     ',/,
     -'      der23=0.d0     ',/,
     -'      der24=0.d0     ',/,
     -'      der34=0.d0     ',/,
     -'      pux12= vex2(1)*r12 ',/,
     -'      pux13= vex2(1)*r13 ',/,
     -'      pux14= vex2(1)*r14 ',/,
     -'      pux23= vex2(2)*r23 ',/,
     -'      pux24= vex2(2)*r24 ',/,
     -'      pux34= vex2(2)*r34 ',/,
     -'      qux12=  exp(-pux12) ',/,
     -'      qux13=  exp(-pux13) ',/,
     -'      qux14=  exp(-pux14) ',/,
     -'      qux23=  exp(-pux23) ',/,
     -'      qux24=  exp(-pux24) ',/,
     -'      qux34=  exp(-pux34) ',/,
     -'      aux12= r12*qux12 ',/,
     -'      aux13= r13*qux13 ',/,
     -'      aux14= r14*qux14 ',/,
     -'      aux23= r23*qux23 ',/,
     -'      aux24= r24*qux24 ',/,
     -'      aux34= r34*qux34 ')
11    format(
     -'      do 1 i=1,mt_2               ',/,
     -'         f12(i) = aux12*f12(i-1)  ',/,
     -'         f13(i) = aux13*f13(i-1)  ',/,
     -'         f14(i) = aux14*f14(i-1)  ',/,
     -'         f23(i) = aux23*f23(i-1)  ',/,
     -'         f24(i) = aux24*f24(i-1)  ',/,
     -'         f34(i) = aux34*f34(i-1)  ',/,
     -'1     continue                 ')
12    format(
     -'      do 2 j=1,is                                           ',/,
     -'         au=fabbb(i1(j),i2(j),i3(j),i4(j),i5(j),i6(j),i7(j))',/,
     -'         ener=ener+cf(j)*au                                 ',/,
     -'         if (iop.eq.1) then                                 ',/,
     -'            call dfabbb(                                    ',/,
     -'     ,            i1(j),i2(j),i3(j),i4(j),i5(j),i6(j),i7(j),',/,
     -'     ,            d12,d13,d14,d23,d24,d34)      ',/,
     -'            der12=der12+cf(j)*d12                          ',/,
     -'            der13=der13+cf(j)*d13                          ',/,
     -'            der14=der14+cf(j)*d14                          ',/,
     -'            der23=der23+cf(j)*d23                          ',/,
     -'            der24=der24+cf(j)*d24                          ',/,
     -'            der34=der34+cf(j)*d34                          ',/,
     -'         endif                                             ',/,
     -'2     continue                                              ',/,
     -'      if (iop.eq.1) then                 ',/,
     -'         der(1)=der12*(1.d0-pux12)*qux12    ',/,
     -'         der(2)=der13*(1.d0-pux13)*qux13    ',/,
     -'         der(3)=der14*(1.d0-pux14)*qux14    ',/,
     -'         der(4)=der23*(1.d0-pux23)*qux23    ',/,
     -'         der(5)=der24*(1.d0-pux24)*qux24    ',/,
     -'         der(6)=der34*(1.d0-pux34)*qux34    ',/,
     -'      endif                                 ',/,
     -'      return                                                ',/,
     -'      end                                                   ')
**********************************************************************
13    format(
     -'*************************************************************',/,
     -'      function fabbb(k1,k2,k3,k4,k5,k6,jt)                   ',/,
     -'      implicit real * 8 (a-h,o-z)                            ',/,
     -'      parameter(mt_2=',i2,')                      ',/,
     -'      common/basu4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),   ',/,
     -'     ,             f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)    ',/,
     -'      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
     -',/,
     -'      if (jt.eq.1) then                                      ',/,
     -'         fabbb=g(k1,k2,k3,k4,k5,k6)                          ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      if (jt.eq.2) then                                      ',/,
     -'         fabbb=(g(k1,k1,k3,k4,k5,k5)+g(k1,k3,k1,k5,k4,k5)+   ',/,
     -'     *   g(k3,k1,k1,k5,k5,k4))                               ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      if (jt.eq.3) then                                      ',/,
     -'         fabbb=g(k1,k2,k3,k4,k5,k6)+g(k2,k3,k1,k6,k4,k5)+    ',/,
     -'     *         g(k3,k1,k2,k5,k6,k4)+g(k3,k2,k1,k6,k5,k4)+    ',/,
     -'     *         g(k1,k3,k2,k5,k4,k6)+g(k2,k1,k3,k4,k6,k5)     ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      end                                                    ')

      write(7,88) mt-2
88    format(
     -/,5x,' subroutine dfabbb ',
     -/,5x,',  (k1,k2,k3,k4,k5,k6,jt,d12,d13,d14,d23,d24,d34) ',
     -/,5x,' implicit real * 8 (a-h,o-z)                          ',
     -/,5x,' parameter(mt_2=',i2,')                            ',
     -/,5x,' common/basu4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),',
     -/,5x,',             f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2) ',
     -/,5x,' g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)',
     -/,5x,' if (jt.eq.1) then                                   ',
     -/,5x,'    d12=dble(k1)*g(k1-1,k2,k3,k4,k5,k6)              ',
     -/,5x,'    d13=dble(k2)*g(k1,k2-1,k3,k4,k5,k6)              ',
     -/,5x,'    d14=dble(k3)*g(k1,k2,k3-1,k4,k5,k6)              ',
     -/,5x,'    d23=dble(k4)*g(k1,k2,k3,k4-1,k5,k6)              ',
     -/,5x,'    d24=dble(k5)*g(k1,k2,k3,k4,k5-1,k6)              ',
     -/,5x,'    d34=dble(k6)*g(k1,k2,k3,k4,k5,k6-1)              ',
     -/,5x,'    return                                           ',
     -/,5x,' endif                                               ',
     -/,5x,' if (jt.eq.2) then                                   ',
     -/,5x,'    d12= dble(k1)*g(k1-1,k1,k3,k4,k5,k5)+            ',
     -/,5x,',        dble(k1)*g(k1-1,k3,k1,k5,k4,k5)+            ',
     -/,5x,',        dble(k3)*g(k3-1,k1,k1,k5,k5,k4)             ',
     -/,5x,'    d13= dble(k1)*g(k1,k1-1,k3,k4,k5,k5)+            ',
     -/,5x,',        dble(k3)*g(k1,k3-1,k1,k5,k4,k5)+            ',
     -/,5x,',        dble(k1)*g(k3,k1-1,k1,k5,k5,k4)             ',
     -/,5x,'    d14= dble(k3)*g(k1,k1,k3-1,k4,k5,k5)+            ',
     -/,5x,',        dble(k1)*g(k1,k3,k1-1,k5,k4,k5)+            ',
     -/,5x,',        dble(k1)*g(k3,k1,k1-1,k5,k5,k4)             ',
     -/,5x,'    d23= dble(k4)*g(k1,k1,k3,k4-1,k5,k5)+            ',
     -/,5x,',        dble(k5)*g(k1,k3,k1,k5-1,k4,k5)+            ',
     -/,5x,',        dble(k5)*g(k3,k1,k1,k5-1,k5,k4)             ',
     -/,5x,'    d24= dble(k5)*g(k1,k1,k3,k4,k5-1,k5)+            ',
     -/,5x,',        dble(k4)*g(k1,k3,k1,k5,k4-1,k5)+            ',
     -/,5x,',        dble(k5)*g(k3,k1,k1,k5,k5-1,k4)             ',
     -/,5x,'    d34= dble(k5)*g(k1,k1,k3,k4,k5,k5-1)+            ',
     -/,5x,',        dble(k5)*g(k1,k3,k1,k5,k4,k5-1)+            ',
     -/,5x,',        dble(k4)*g(k3,k1,k1,k5,k5,k4-1)             ',
     -/,5x,'    return                                           ',
     -/,5x,' endif                                               ',
     -/,5x,' if (jt.eq.3) then                                   ',
     -/,5x,'    d12=dble(k1)*g(k1-1,k2,k3,k4,k5,k6)+             ',
     -/,5x,',       dble(k2)*g(k2-1,k3,k1,k6,k4,k5)+             ',
     -/,5x,',       dble(k3)*g(k3-1,k1,k2,k5,k6,k4)+             ',
     -/,5x,',       dble(k3)*g(k3-1,k2,k1,k6,k5,k4)+             ',
     -/,5x,',       dble(k1)*g(k1-1,k3,k2,k5,k4,k6)+             ',
     -/,5x,',       dble(k2)*g(k2-1,k1,k3,k4,k6,k5)              ',
     -/,5x,'    d13=dble(k2)*g(k1,k2-1,k3,k4,k5,k6)+             ',
     -/,5x,',       dble(k3)*g(k2,k3-1,k1,k6,k4,k5)+             ',
     -/,5x,',       dble(k1)*g(k3,k1-1,k2,k5,k6,k4)+             ',
     -/,5x,',       dble(k2)*g(k3,k2-1,k1,k6,k5,k4)+             ',
     -/,5x,',       dble(k3)*g(k1,k3-1,k2,k5,k4,k6)+             ',
     -/,5x,',       dble(k1)*g(k2,k1-1,k3,k4,k6,k5)              ',
     -/,5x,'    d14=dble(k3)*g(k1,k2,k3-1,k4,k5,k6)+             ',
     -/,5x,',       dble(k1)*g(k2,k3,k1-1,k6,k4,k5)+             ',
     -/,5x,',       dble(k2)*g(k3,k1,k2-1,k5,k6,k4)+             ',
     -/,5x,',       dble(k1)*g(k3,k2,k1-1,k6,k5,k4)+             ',
     -/,5x,',       dble(k2)*g(k1,k3,k2-1,k5,k4,k6)+             ',
     -/,5x,',       dble(k3)*g(k2,k1,k3-1,k4,k6,k5)              ',
     -/,5x,'    d23=dble(k4)*g(k1,k2,k3,k4-1,k5,k6)+             ',
     -/,5x,',       dble(k6)*g(k2,k3,k1,k6-1,k4,k5)+             ',
     -/,5x,',       dble(k5)*g(k3,k1,k2,k5-1,k6,k4)+             ',
     -/,5x,',       dble(k6)*g(k3,k2,k1,k6-1,k5,k4)+             ',
     -/,5x,',       dble(k5)*g(k1,k3,k2,k5-1,k4,k6)+             ',
     -/,5x,',       dble(k4)*g(k2,k1,k3,k4-1,k6,k5)              ',
     -/,5x,'    d24=dble(k5)*g(k1,k2,k3,k4,k5-1,k6)+             ',
     -/,5x,',       dble(k4)*g(k2,k3,k1,k6,k4-1,k5)+             ',
     -/,5x,',       dble(k6)*g(k3,k1,k2,k5,k6-1,k4)+             ',
     -/,5x,',       dble(k5)*g(k3,k2,k1,k6,k5-1,k4)+             ',
     -/,5x,',       dble(k4)*g(k1,k3,k2,k5,k4-1,k6)+             ',
     -/,5x,',       dble(k6)*g(k2,k1,k3,k4,k6-1,k5)              ',
     -/,5x,'    d34=dble(k6)*g(k1,k2,k3,k4,k5,k6-1)+             ',
     -/,5x,',       dble(k5)*g(k2,k3,k1,k6,k4,k5-1)+             ',
     -/,5x,',       dble(k4)*g(k3,k1,k2,k5,k6,k4-1)+             ',
     -/,5x,',       dble(k4)*g(k3,k2,k1,k6,k5,k4-1)+             ',
     -/,5x,',       dble(k6)*g(k1,k3,k2,k5,k4,k6-1)+             ',
     -/,5x,',       dble(k5)*g(k2,k1,k3,k4,k6,k5-1)              ',
     -/,5x,'    return                                           ',
     -/,5x,' endif                                               ',
     -/,5x,' end')

      end
****************************************************************
**********************************************************************
C     *** MINIM *******************************************************
C     *
C     *   MINIMIZATION OF A FUNCTION OF SEVERAL VARIABLES
C     *   USING NR'S ALGORITHM
C     *
C     ******************************************************************
C
      subroutine minim(f,x,y,n,lim,x0,ier)
      implicit real*8 (a-h,o-z)
      external f
c     dimension x(*),x0(*),der(3)
      dimension x(*),x0(*),der(6)
c
cesar eps=1.d-6
      eps=1.d-6
cesar dereps=1.d-3
      dereps=1.d-2
      dereps2=1.d-3
      pasoeps=1.d-4
      do 5 i=1,n
         x0(i)=x(i)
  5   continue
      ymin=1.d100
                        ymin=f(x)
                        y=ymin
c                       goto 88
      if(n.eq.1) then
           xmin1=x(1)
c          do 1 x1=x0(1)-.6d0,x0(1)+.6d0,.3d0
           do 1 x1=x0(1)-.5d0,x0(1)+.5d0,.5d0
                x(1)=x1
                y=f(x)
                if(y.lt.ymin) then
                  ymin=y
                  xmin1=x(1)
                endif
1          continue
           x(1)=xmin1
      elseif(n.eq.2) then
           xmin1=x(1)
           xmin2=x(2)
           do 2 x1=x0(1)-.5d0,x0(1)+.5d0,.5d0
           do 2 x2=x0(2)-.5d0,x0(2)+.5d0,.5d0
                x(1)=x1
                x(2)=x2
                y=f(x)
                if(y.lt.ymin) then
                  ymin=y
                  xmin1=x(1)
                  xmin2=x(2)
                endif
2          continue
           x(1)=xmin1
           x(2)=xmin2
      elseif(n.eq.3) then
           xmin1=x(1)
           xmin2=x(2)
           xmin3=x(3)
           do 3 x1=x0(1)-.5d0,x0(1)+.5d0,.5d0
           do 3 x2=x0(2)-.5d0,x0(2)+.5d0,.5d0
           do 3 x3=x0(3)-.5d0,x0(3)+.5d0,.5d0
                x(1)=x1
                x(2)=x2
                x(3)=x3
                y=f(x)
                if(y.lt.ymin) then
                  ymin=y
                  xmin1=x(1)
                  xmin2=x(2)
                  xmin3=x(3)
                endif
3          continue
           x(1)=xmin1
           x(2)=xmin2
           x(3)=xmin3
      elseif(n.eq.6) then
           xmin1=x(1)
           xmin2=x(2)
           xmin3=x(3)
           xmin4=x(4)
           xmin5=x(5)
           xmin6=x(6)
           do 6 x1=x0(1)-.5d0,x0(1)+.5d0,.5d0
           do 6 x2=x0(2)-.5d0,x0(2)+.5d0,.5d0
           do 6 x3=x0(3)-.5d0,x0(3)+.5d0,.5d0
           do 6 x4=x0(4)-.5d0,x0(4)+.5d0,.5d0
           do 6 x5=x0(5)-.5d0,x0(5)+.5d0,.5d0
           do 6 x6=x0(6)-.5d0,x0(6)+.5d0,.5d0
                x(1)=x1
                x(2)=x2
                x(3)=x3
                x(4)=x4
                x(5)=x5
                x(6)=x6
cc              aux=(x(1)-x0(1))**2+(x(2)-x0(2))**2+(x(3)-x0(3))**2
cc              aux=aux+(x(4)-x0(4))**2+(x(5)-x0(5))**2+(x(6)-x0(6))**2
cc              if (aux.gt.1.d-5) then
                   y=f(x)
                   if(y.lt.ymin) then
                     ymin=y
                     xmin1=x(1)
                     xmin2=x(2)
                     xmin3=x(3)
                     xmin4=x(4)
                     xmin5=x(5)
                     xmin6=x(6)
                   endif
c               endif
6       continue
        x(1)=xmin1
        x(2)=xmin2
        x(3)=xmin3
        x(4)=xmin4
        x(5)=xmin5
        x(6)=xmin6
      endif
88                         continue
*----->
      do 20 j=1,lim
*-->
      do 10 i=1,n
         xmin=x(i)
         x(i)=x(i)+eps
         y2=f(x)
         if(y2.lt.ymin) then
           ymin=y2
           xmin=x(i)
         endif
         x(i)=x(i)-2.d0*eps
         y3=f(x)
         x(i)=x(i)+eps
         if(y3.lt.ymin) then
           ymin=y3
           xmin=x(i)
         endif
         der1=(y2-y3)/(2.d0*eps)
         derder=(y2-y)/eps
         derizq=(y-y3)/eps
         der2=(y2+y3-2.d0*y)/(eps*eps)
         paso0=-y/der1
         paso1=-der1/der2
         if(dabs(der1).gt.dereps2) then
           paso=paso0
         else
           paso=paso1
         endif
         x(i)=x(i)+paso
         y=f(x)
*- if -->
         if(y.lt.ymin) then
             ymin=y
             xmin=x(i)
         elseif(dabs(paso).gt.pasoeps) then
c        else
            do 100 k=1,2
               x(i)=x(i)-paso
               paso=paso*0.1d0
               x(i)=x(i)+paso
               y=f(x)
              if(y.lt.ymin) then
                ymin=y
                xmin=x(i)
              endif
100         continue
         endif
*- if --<
         der(i)=der1
         x(i)=xmin
         y=ymin
 10   continue
*--<
      dif=0.d0
      deri=0.d0
      do 15 i=1,n
         dif=dif+dabs(x0(i)-x(i))
         deri=deri+der(i)**2
         x0(i)=x(i)
 15   continue
      deri=dsqrt(deri)
      if(dif.lt.pasoeps.and.deri.lt.dereps) then
         ier=0
         return
      endif
 20   continue
      print*,' Warning: Number of iterations exceded'
      ier=1
      return
      end
************************************************************************
c
      subroutine gauss(a,b,n,ks)
      implicit real * 8 (a-h,o-z)
c
      include 'dimensions.inc'
c
      dimension a(nexp*nexp),b(nexp)
      tol=0.d0
      ks=0
      jj=-n
      do 10 j=1,n
      jy=j+1
      jj=jj+n+1
      biga=0.d0
      it=jj-j
      do 30 i=j,n
      ij=it+i
      if(dabs(biga)-dabs(a(ij))) 20,30,30
   20 biga=a(ij)
      imax=i
   30 continue
      if(dabs(biga)-tol) 35,35,40
   35 ks=1
      return
   40 i1=j+n*(j-2)
      it=imax-j
      do 50 k=j,n
      i1=i1+n
      i2=i1+it
      save=a(i1)
      a(i1)=a(i2)
      a(i2)=save
   50 a(i1)=a(i1)/biga
      save=b(imax)
      b(imax)=b(j)
      b(j)=save/biga
      if(j-n) 55,70,55
   55 iqs=n*(j-1)
      do 10 ix=jy,n
      ixj=iqs+ix
      it=j-ix
      do 60 jx=jy,n
      ixjx=n*(jx-1)+ix
      jjx=ixjx+it
   60 a(ixjx)=a(ixjx)-(a(ixj)*a(jjx))
   10 b(ix)=b(ix)-(b(j)*a(ixj))
   70 ny=n-1
      it=n*n
      do 80 j=1,ny
      ia=it-j
      ib=n-j
      ic=n
      do 80 k=1,j
      b(ib)=b(ib)-a(ia)*b(ic)
      ia=ia-n
   80 ic=ic-1
      return
      end
*************************************************************************
      subroutine scr2(nombre,e0,n,mt,vec2,rms,emax)
****************************************************************
*     Escribe el programa fortran de las diat¢micas
*     (en la unidad 7)
****************************************************************
      implicit real*8 (a-h,o-z)
      character*6 nombre
      include 'dimensions.inc'
      common sm(nexp)
      dimension vec2(2)
      write(7,3010) nombre,n,rms*627.51d0,emax*627.51d0,mt,
     ,              e0,vec2(1),vec2(2)
 3010 format(72('*'),/,6x,'subroutine ',a6,'(r,ener,der,iop)',/,
     -72('*'),/,
     -'*     This subroutine computes ',
     -'the energies of a diatomic potential ',/,
     -'*     fitted to',i6,' points',/,
     -'*     rms = ', f15.8,' kcal/mol',/,
     -'*     emax = ', f15.8,' kcal/mol',/,
     -72('*'),/,
     -6x,'implicit real*8 (a-h,o-z)',/,
     -6x,'parameter(mt=',i3,')',/,
     -6x,'parameter(e0=',d14.7,')',/,
     -6x,'parameter(vex1=',d23.15,')',/,
     -6x,'parameter(vex2=',d23.15,')',/,
     -6x,'dimension cf(mt)')

      do 17 k=1,mt
         write(7,3012) k,sm(k)
   17 continue
 3012 format(6x,'data cf(',i3,')/',d23.15,'/')

      write(7,3)
3     format(
     -6x,'der = 0.d0',/,
     -6x,'aux = 1.d0/r',/,
     -6x,'bux = dexp(-vex2*r)*aux',/,
     -6x,'cux = dexp(-vex1*r)')
      write(7,3014)
 3014 format(
     -6x,'ener=e0+cf(1)*bux',/,
     -6x,'dux=1.d0',/,
     -6x,'eux=r*cux')

      write(7,3015)
 3015 format(
     -6x,'do 1 i=2,mt',/,
     -9x,'if (iop.eq.1) der=der+(i-1)*cf(i)*dux',/,
     -9x,'dux=dux*eux',/,
     -9x,'ener=ener+cf(i)*dux',
     -/,4x,'1 continue',
     -/,6x,'if (iop.eq.1) then',
     -/,9x,'der=der*(1.d0-vex1*r)*cux',
     -/,9x,'der=der-cf(1)*(vex2+aux)*bux',
     -/,6x,'endif',
     -/,6x,'return',/,6x,'end')
      return
      end
****************************************************************
      subroutine saaa(nombre,n,mt,vex1,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie AAA
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*6 nombre
      include 'dimensions.inc'
      common/ipot31/ipa31(nexp3,4),nfas31(igrad3)
      common sm(nexp)
      is=nfas31(mt)
      write(7,3010) nombre,n,rms*627.51d0,emax*627.51d0,is,mt-1,vex1
3010  format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r23,ener,der,iop)          ',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 3D PES      ',/,
     -'*     for the AAA system class fitted to ',i4,' points       ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_1=',i2,')                      ',/,
     -'      parameter(vex1=',d23.15,')                             ',/,
     -'      dimension der(3)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),cf(is)           ',/,
     -'      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            ')
      do 17 k=1,is
      write(7,1)
     , k,sm(k),k,ipa31(k,1),k,ipa31(k,2),k,ipa31(k,3),k,ipa31(k,4)
   17 continue
 1    format(6x,'data cf(',i3,')/',d23.15,'/',/,
     -6x,'data i1(',i3,')/',i2,
     -'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,'/,i4(',i3,')/',i2,'/')

      write(7,3)
3     format(
     -6x,'ener=0.d0',/,
     -6x,'der12=0.d0',/,
     -6x,'der13=0.d0',/,
     -6x,'der23=0.d0',/,
     -6x,'pux12=vex1*r12',/,
     -6x,'pux13=vex1*r13',/,
     -6x,'pux23=vex1*r23',/,
     -6x,'qux12=dexp(-pux12)',/,
     -6x,'qux13=dexp(-pux13)',/,
     -6x,'qux23=dexp(-pux23)',/,
     -6x,'bux12=r12*qux12',/,
     -6x,'bux13=r13*qux13',/,
     -6x,'bux23=r23*qux23',/,
     -6x,'do 1 i=1,mt_1',/,
     -9x,'f12(i)=f12(i-1)*bux12',/,
     -9x,'f13(i)=f13(i-1)*bux13',/,
     -9x,'f23(i)=f23(i-1)*bux23',/,
     -'1',5x,'continue')
      write(7,4)   ! is
4     format(
     -/,6x,'do 2 l=1,is',
     -/,6x,'   if (i4(l).eq.1) then                           ',
     -/,6x,'      aux=f12(i1(l))*f13(i2(l))*f23(i3(l))        ',
     -/,6x,'      if (iop.eq.1) then                          ',
     -/,6x,'         dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))',
     -/,6x,'         dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))',
     -/,6x,'         dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)',
     -/,6x,'      endif                                    ',
     -/,6x,'      elseif (i4(l).eq.3) then                        ',
     -/,6x,'         aux1=f12(i1(l))*f13(i2(l))*f23(i3(l))    ',
     -/,6x,'         aux2=f12(i3(l))*f13(i1(l))*f23(i2(l))    ',
     -/,6x,'         aux3=f12(i3(l))*f13(i2(l))*f23(i1(l))    ',
     -/,6x,'         aux=aux1+aux2+aux3                    ',
     -/,6x,'         if (iop.eq.1) then                    ',
     -/,6x,'            dux1=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l)) ',
     -/,6x,'            dux2=i3(l)*f12(i3(l)-1)*f13(i1(l))*f23(i2(l)) ',
     -/,6x,'            dux3=i2(l)*f12(i2(l)-1)*f13(i3(l))*f23(i1(l)) ',
     -/,6x,'            dux12=dux1+dux2+dux3                          ',
     -/,6x,'            dux1=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l)) ',
     -/,6x,'            dux2=i1(l)*f12(i3(l))*f13(i1(l)-1)*f23(i2(l)) ',
     -/,6x,'            dux3=i3(l)*f12(i2(l))*f13(i3(l)-1)*f23(i1(l)) ',
     -/,6x,'            dux13=dux1+dux2+dux3                          ',
     -/,6x,'            dux1=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1) ',
     -/,6x,'            dux2=i2(l)*f12(i3(l))*f13(i1(l))*f23(i2(l)-1) ',
     -/,6x,'            dux3=i1(l)*f12(i2(l))*f13(i3(l))*f23(i1(l)-1) ',
     -/,6x,'            dux23=dux1+dux2+dux3                          ',
     -/,6x,'         endif                                 ',
     -/,6x,'         elseif (i4(l).eq.6) then                     ',
     -/,6x,'            aux1=f12(i1(l))*f13(i2(l))*f23(i3(l)) ',
     -/,6x,'            aux2=f12(i1(l))*f13(i3(l))*f23(i2(l)) ',
     -/,6x,'            aux3=f12(i2(l))*f13(i1(l))*f23(i3(l)) ',
     -/,6x,'            aux4=f12(i2(l))*f13(i3(l))*f23(i1(l)) ',
     -/,6x,'            aux5=f12(i3(l))*f13(i1(l))*f23(i2(l)) ',
     -/,6x,'            aux6=f12(i3(l))*f13(i2(l))*f23(i1(l)) ',
     -/,6x,'            aux=aux1+aux2+aux3+aux4+aux5+aux6  ',
     -/,6x,'            if (iop.eq.1) then                 ',
     -/,6x,'             dux1=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))',
     -/,6x,'             dux2=i1(l)*f12(i1(l)-1)*f13(i3(l))*f23(i2(l))',
     -/,6x,'             dux3=i2(l)*f12(i2(l)-1)*f13(i1(l))*f23(i3(l))',
     -/,6x,'             dux4=i2(l)*f12(i2(l)-1)*f13(i3(l))*f23(i1(l))',
     -/,6x,'             dux5=i3(l)*f12(i3(l)-1)*f13(i1(l))*f23(i2(l))',
     -/,6x,'             dux6=i3(l)*f12(i3(l)-1)*f13(i2(l))*f23(i1(l))',
     -/,6x,'             dux12=dux1+dux2+dux3+dux4+dux5+dux6          ',
     -/,6x,'             dux1=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))',
     -/,6x,'             dux2=i3(l)*f12(i1(l))*f13(i3(l)-1)*f23(i2(l))',
     -/,6x,'             dux3=i1(l)*f12(i2(l))*f13(i1(l)-1)*f23(i3(l))',
     -/,6x,'             dux4=i3(l)*f12(i2(l))*f13(i3(l)-1)*f23(i1(l))',
     -/,6x,'             dux5=i1(l)*f12(i3(l))*f13(i1(l)-1)*f23(i2(l))',
     -/,6x,'             dux6=i2(l)*f12(i3(l))*f13(i2(l)-1)*f23(i1(l))',
     -/,6x,'             dux13=dux1+dux2+dux3+dux4+dux5+dux6          ',
     -/,6x,'             dux1=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)',
     -/,6x,'             dux2=i2(l)*f12(i1(l))*f13(i3(l))*f23(i2(l)-1)',
     -/,6x,'             dux3=i3(l)*f12(i2(l))*f13(i1(l))*f23(i3(l)-1)',
     -/,6x,'             dux4=i1(l)*f12(i2(l))*f13(i3(l))*f23(i1(l)-1)',
     -/,6x,'             dux5=i2(l)*f12(i3(l))*f13(i1(l))*f23(i2(l)-1)',
     -/,6x,'             dux6=i1(l)*f12(i3(l))*f13(i2(l))*f23(i1(l)-1)',
     -/,6x,'             dux23=dux1+dux2+dux3+dux4+dux5+dux6          ',
     -/,6x,'            endif                                       ',
     -/,6x,'   endif                                       ',
     -/,6x,'   ener=ener+cf(l)*aux                         ',
     -/,6x,'   if (iop.eq.1) then                          ',
     -/,6x,'      der12=der12+cf(l)*dux12                  ',
     -/,6x,'      der13=der13+cf(l)*dux13                  ',
     -/,6x,'      der23=der23+cf(l)*dux23                  ',
     -/,6x,'   endif                                       ',
     -/,4x,'2 continue                                     ')
      write(7,5)
5     format(
     -/,6x,' if (iop.eq.1) then                       ',
     -/,6x,'    der(1)=der12*(1.d0-pux12)*qux12       ',
     -/,6x,'    der(2)=der13*(1.d0-pux13)*qux13       ',
     -/,6x,'    der(3)=der23*(1.d0-pux23)*qux23       ',
     -/,6x,' endif                                    ',
     -/,6x,'return                                    ',
     -/,6x,'end')
      return
      end
************************************************************************
****************************************************************
*     Escribe el programa fortran que calcula la superficie ABB
c     necesita la matriz sm
****************************************************************
      subroutine sabb(nombre,n,mt,vex2,rms,emax)
      implicit  real * 8 (a-h,o-z)
      character*6 nombre
      include 'dimensions.inc'
      common/ipot32/ipa32(nexp3,4),nfas32(igrad3)
      common sm(nexp)
      dimension vex2(2)
      is=nfas32(mt)
      write(7,3010) nombre,n,rms*627.51d0,emax*627.51d0,is,mt-1,
     ,vex2(1),vex2(2)
3010  format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r23,ener,der,iop)          ',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 3D PES      ',/,
     -'*     for the ABB system class fitted to ',i4,' points       ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_1=',i2,')                      ',/,
     -'      parameter(vex1=',d23.15,')                             ',/,
     -'      parameter(vex2=',d23.15,')                             ',/,
     -'      dimension der(3)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),cf(is)           ',/,
     -'      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            ')


      do 17 k=1,is
      write(7,1)
     , k,sm(k),k,ipa32(k,1),k,ipa32(k,2),k,ipa32(k,3),k,ipa32(k,4)
   17 continue
 1    format(6x,'data cf(',i3,')/',d23.15,'/',/,
     -6x,'data i1(',i3,')/',i2,
     -'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,'/,i4(',i3,')/',i2,'/')

      write(7,3)
3     format(
     -6x,'ener=0.d0',/,
     -6x,'der12=0.d0',/,
     -6x,'der13=0.d0',/,
     -6x,'der23=0.d0',/,
     -6x,'pux12=vex1*r12',/,
     -6x,'pux13=vex1*r13',/,
     -6x,'pux23=vex2*r23',/,
     -6x,'qux12=dexp(-pux12)',/,
     -6x,'qux13=dexp(-pux13)',/,
     -6x,'qux23=dexp(-pux23)',/,
     -6x,'bux12=r12*qux12',/,
     -6x,'bux13=r13*qux13',/,
     -6x,'bux23=r23*qux23',/,
     -6x,'do 1 i=1,mt_1',/,
     -9x,'f12(i)=f12(i-1)*bux12',/,
     -9x,'f13(i)=f13(i-1)*bux13',/,
     -9x,'f23(i)=f23(i-1)*bux23',/,
     -'1',5x,'continue')

      write(7,4)
4     format(
     -/,6x,'do 2 l=1,is',
     -/,6x,'   if (i4(l).eq.1) then                           ',
     -/,6x,'      aux=f12(i1(l))*f13(i3(l))*f23(i2(l))        ',
     -/,6x,'      if (iop.eq.1) then                          ',
     -/,6x,'         dux12=i1(l)*f12(i1(l)-1)*f13(i3(l))*f23(i2(l))',
     -/,6x,'         dux13=i3(l)*f12(i1(l))*f13(i3(l)-1)*f23(i2(l))',
     -/,6x,'         dux23=i2(l)*f12(i1(l))*f13(i3(l))*f23(i2(l)-1)',
     -/,6x,'      endif                          ',
     -/,6x,'     else                                      ',
     -/,6x,'      aux1=f12(i1(l))*f13(i3(l))                 ',
     -/,6x,'      aux2=f12(i3(l))*f13(i1(l))                 ',
     -/,6x,'      aux=(aux1+aux2)*f23(i2(l))                ',
     -/,6x,'      if (iop.eq.1) then                          ',
     -/,6x,'         dux23=(aux1+aux2)*i2(l)*f23(i2(l)-1)          ',
     -/,6x,'         dux1=i1(l)*f12(i1(l)-1)*f13(i3(l))            ',
     -/,6x,'         dux2=i3(l)*f12(i3(l)-1)*f13(i1(l))            ',
     -/,6x,'         dux12=(dux1+dux2)*f23(i2(l))                  ',
     -/,6x,'         dux1=i3(l)*f12(i1(l))*f13(i3(l)-1)            ',
     -/,6x,'         dux2=i1(l)*f12(i3(l))*f13(i1(l)-1)            ',
     -/,6x,'         dux13=(dux1+dux2)*f23(i2(l))                  ',
     -/,6x,'      endif                                       ',
     -/,6x,'   endif                                       ',
     -/,6x,'   ener=ener+cf(l)*aux                         ',
     -/,6x,'   if (iop.eq.1) then                          ',
     -/,6x,'      der12=der12+cf(l)*dux12                  ',
     -/,6x,'      der13=der13+cf(l)*dux13                  ',
     -/,6x,'      der23=der23+cf(l)*dux23                  ',
     -/,6x,'   endif                                       ',
     -/,4x,'2 continue                                     ',
     -/,6x,'if (iop.eq.1) then                          ',
     -/,6x,'   der(1)=der12*(1.d0-pux12)*qux12        ',
     -/,6x,'   der(2)=der13*(1.d0-pux13)*qux13        ',
     -/,6x,'   der(3)=der23*(1.d0-pux23)*qux23        ',
     -/,6x,'endif                         ',
     -/,6x,'return',/,6x,'end')

      return
      end
************************************************************************
****************************************************************
      subroutine sabc(nombre,n,mt,vex3,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie ABC
*     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*6 nombre
      include 'dimensions.inc'
      common/ipot33/ipa33(nexp3,3),nfas33(igrad3)
      common sm(nexp)
      dimension vex3(3)
      is=nfas33(mt)
      write(7,3010) nombre,n,rms*627.51d0,emax*627.51d0,is,mt-1,
     ,vex3(1),vex3(2),vex3(3)
3010  format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r23,ener,der,iop)          ',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 3D PES      ',/,
     -'*     for the ABC system class fitted to ',i4,' points       ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_1=',i2,')                      ',/,
     -'      parameter(vex1=',d23.15,')                             ',/,
     -'      parameter(vex2=',d23.15,')                             ',/,
     -'      parameter(vex3=',d23.15,')                             ',/,
     -'      dimension der(3)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),cf(is)                  ',/,
     -'      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            ')
      do 17 k=1,is
      write(7,1) k,sm(k),k,ipa33(k,1),k,ipa33(k,2),k,ipa33(k,3)
   17 continue
 1    format(6x,'data cf(',i3,')/',d23.15,'/',/,
     -6x,'data i1(',i3,')/',i2,
     -'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,'/')

      write(7,3)
3     format(
     -6x,'ener=0.d0',/,
     -6x,'der12=0.d0',/,
     -6x,'der13=0.d0',/,
     -6x,'der23=0.d0',/,
     -6x,'pux12=vex1*r12',/,
     -6x,'pux13=vex2*r13',/,
     -6x,'pux23=vex3*r23',/,
     -6x,'qux12=dexp(-pux12)',/,
     -6x,'qux13=dexp(-pux13)',/,
     -6x,'qux23=dexp(-pux23)',/,
     -6x,'bux12=r12*qux12',/,
     -6x,'bux13=r13*qux13',/,
     -6x,'bux23=r23*qux23',/,
     -6x,'do 1 i=1,mt_1',/,
     -9x,'f12(i)=f12(i-1)*bux12',/,
     -9x,'f13(i)=f13(i-1)*bux13',/,
     -9x,'f23(i)=f23(i-1)*bux23',/,
     -'1',5x,'continue')

      write(7,4)
4     format(
     -/,6x,'do 2 l=1,is',
     -/,6x,'aux=f12(i1(l))*f13(i2(l))*f23(i3(l))           ',
     -/,6x,'ener=ener+cf(l)*aux                         ',
     -/,6x,'if (iop.eq.1) then                          ',
     -/,6x,'   dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))   ',
     -/,6x,'   dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))   ',
     -/,6x,'   dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)   ',
     -/,6x,'   der12=der12+cf(l)*dux12                          ',
     -/,6x,'   der13=der13+cf(l)*dux13                          ',
     -/,6x,'   der23=der23+cf(l)*dux23                          ',
     -/,6x,'endif                         ',
     -/,4x,'2 continue                                     ',
     -/,6x,'if (iop.eq.1) then                          ',
     -/,6x,'   der(1)=der12*(1.d0-pux12)*qux12        ',
     -/,6x,'   der(2)=der13*(1.d0-pux13)*qux13        ',
     -/,6x,'   der(3)=der23*(1.d0-pux23)*qux23        ',
     -/,6x,'endif                         ',
     -/,6x,'return',/,6x,'end')

      return
      end
****************************************************************
      subroutine sabcd(nombre,n,mt,vex6,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie ABCD
c     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*4 nombre
      include 'dimensions.inc'
      common/ipot45/ipa45(nexp45,6),nfas45(igrad4)
      common sm(nexp)
      dimension vex6(6)
      is=nfas45(mt)

      write(7,1) nombre,n,rms*627.51d0,emax*627.51d0,
     -is,mt-2,(vex6(i),i=1,6)
      do 2 k=1,is
         write(7,3) k,sm(k),k,ipa45(k,1),k,ipa45(k,2),k,ipa45(k,3),
     -                      k,ipa45(k,4),k,ipa45(k,5),k,ipa45(k,6)
2     continue
      write(7,10)
      write(7,11)
      write(7,12)
1     format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 6D PES      ',/,
     -'*     for the ABCD system class fitted to ',i4,' points      ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_2=',i2,')                      ',/,
     -'      dimension vex6(6)                                      ',/,
     -'      dimension der(6)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),cf(is)
     ,       ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      data vex6(1)/',d23.15,'/                               ',/,
     -'      data vex6(2)/',d23.15,'/                               ',/,
     -'      data vex6(3)/',d23.15,'/                               ',/,
     -'      data vex6(4)/',d23.15,'/                               ',/,
     -'      data vex6(5)/',d23.15,'/                               ',/,
     -'      data vex6(6)/',d23.15,'/                               ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/         ',/,
     -'      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/            ',/,
     -'      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            ')



3     format(
     -'      data cf(',i3,')/',d23.15,'/                           ',/,
     -'      data i1(',i3,')/',i2,'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,
     -'/'   ,/,
     -'      data i4(',i3,')/',i2,'/,i5(',i3,')/',i2,'/,i6(',i3,')/',i2,
     -'/')

10    format(
     -'      ener=0.d0          ',/,
     -'      der12=0.d0         ',/,
     -'      der13=0.d0         ',/,
     -'      der14=0.d0         ',/,
     -'      der23=0.d0         ',/,
     -'      der24=0.d0         ',/,
     -'      der34=0.d0         ',/,
     -'      pux12= vex6(1)*r12 ',/,
     -'      pux13= vex6(2)*r13 ',/,
     -'      pux14= vex6(3)*r14 ',/,
     -'      pux23= vex6(4)*r23 ',/,
     -'      pux24= vex6(5)*r24 ',/,
     -'      pux34= vex6(6)*r34 ',/,
     -'      qux12= dexp(-pux12) ',/,
     -'      qux13= dexp(-pux13) ',/,
     -'      qux14= dexp(-pux14) ',/,
     -'      qux23= dexp(-pux23) ',/,
     -'      qux24= dexp(-pux24) ',/,
     -'      qux34= dexp(-pux34) ',/,
     -'      aux12= r12*qux12 ',/,
     -'      aux13= r13*qux13 ',/,
     -'      aux14= r14*qux14 ',/,
     -'      aux23= r23*qux23 ',/,
     -'      aux24= r24*qux24 ',/,
     -'      aux34= r34*qux34 ')
11    format(
     -'      do 1 i=1,mt_2               ',/,
     -'         f12(i) = aux12*f12(i-1)  ',/,
     -'         f13(i) = aux13*f13(i-1)  ',/,
     -'         f14(i) = aux14*f14(i-1)  ',/,
     -'         f23(i) = aux23*f23(i-1)  ',/,
     -'         f24(i) = aux24*f24(i-1)  ',/,
     -'         f34(i) = aux34*f34(i-1)  ',/,
     -'1     continue                 ')
12    format(
     -'      do 2 j=1,is                               ',/,
     -'         f1213=f12(i1(j))*f13(i2(j))            ',/,
     -'         f1423=f14(i3(j))*f23(i4(j))            ',/,
     -'         f2434=f24(i5(j))*f34(i6(j))            ',/,
     -'         au=f1213*f1423*f2434                ',/,
     -'         ener=ener+cf(j)*au                     ',/,
     -'         if (iop.eq.1) then                     ',/,
     -'            f14232434=f1423*f2434               ',/,
     -'            f12132434=f1213*f2434               ',/,
     -'            f12131423=f1213*f1423               ',/,
     -'            d12=i1(j)*f12(i1(j)-1)*f13(i2(j))*f14232434  ',/,
     -'            d13=i2(j)*f12(i1(j))*f13(i2(j)-1)*f14232434  ',/,
     -'            d14=i3(j)*f14(i3(j)-1)*f23(i4(j))*f12132434  ',/,
     -'            d23=i4(j)*f14(i3(j))*f23(i4(j)-1)*f12132434  ',/,
     -'            d24=i5(j)*f12131423*f24(i5(j)-1)*f34(i6(j))  ',/,
     -'            d34=i6(j)*f12131423*f24(i5(j))*f34(i6(j)-1)  ',/,
     -'            der12=der12+cf(j)*d12                          ',/,
     -'            der13=der13+cf(j)*d13                          ',/,
     -'            der14=der14+cf(j)*d14                          ',/,
     -'            der23=der23+cf(j)*d23                          ',/,
     -'            der24=der24+cf(j)*d24                          ',/,
     -'            der34=der34+cf(j)*d34                          ',/,
     -'         endif                                             ',/,
     -'2     continue                                              ',/,
     -'      if (iop.eq.1) then                 ',/,
     -'         der(1)=der12*(1.d0-pux12)*qux12    ',/,
     -'         der(2)=der13*(1.d0-pux13)*qux13    ',/,
     -'         der(3)=der14*(1.d0-pux14)*qux14    ',/,
     -'         der(4)=der23*(1.d0-pux23)*qux23    ',/,
     -'         der(5)=der24*(1.d0-pux24)*qux24    ',/,
     -'         der(6)=der34*(1.d0-pux34)*qux34    ',/,
     -'      endif                                        ',/,
     -'      return                                                ',/,
     -'      end                                                   ')
      return
      end
***********************************************************************
**********************************************************************
      subroutine gfaabb(n,mt,vex3)
**********************************************************************
*     devuelve el valor de la fas j para cada punto i en xx(j,i)
**********************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot43/ipa43(nexp43,7),nfas43(igrad4)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      dimension vex3(3)
      do 1 i=1,n
         f12(0) = 1.d0
         f13(0) = 1.d0
         f14(0) = 1.d0
         f23(0) = 1.d0
         f24(0) = 1.d0
         f34(0) = 1.d0
c        vex3(1)--> AA
c        vex3(2)--> AB
c        vex3(3)--> BB
         aux12= r12(i)*dexp(-vex3(1)*r12(i))
         aux13= r13(i)*dexp(-vex3(2)*r13(i))
         aux14= r14(i)*dexp(-vex3(2)*r14(i))
         aux23= r23(i)*dexp(-vex3(2)*r23(i))
         aux24= r24(i)*dexp(-vex3(2)*r24(i))
         aux34= r34(i)*dexp(-vex3(3)*r34(i))
*     como el grado es mt, la potencia de orden mayor es mt-3
         do 2 k=1,mt-2
            f12(k) = aux12*f12(k-1)
            f13(k) = aux13*f13(k-1)
            f14(k) = aux14*f14(k-1)
            f23(k) = aux23*f23(k-1)
            f24(k) = aux24*f24(k-1)
            f34(k) = aux34*f34(k-1)
2        continue
         do 3 j=1,nfas43(mt)
            i1=ipa43(j,1)
            i2=ipa43(j,2)
            i3=ipa43(j,3)
            i4=ipa43(j,4)
            i5=ipa43(j,5)
            i6=ipa43(j,6)
            i7=ipa43(j,7)
            xx(j,i)=faabb(i1,i2,i3,i4,i5,i6,i7)
 3       continue
 1    continue
      end
**********************************************************************
      function faabb(k1,k2,k3,k4,k5,k6,jt)
      implicit real*8(a-h,o-z)
      parameter (igrad4=12)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
      if (jt.eq.1) then
         faabb=g(k1,k2,k3,k4,k5,k6)
         return
      endif
      if (jt.eq.2) then
         faabb=g(k1,k2,k3,k4,k5,k6)+g(k1,k4,k5,k2,k3,k6)
         return
      endif
      if (jt.eq.3) then
         faabb=g(k1,k2,k3,k4,k5,k6)+g(k1,k3,k2,k5,k4,k6)+
     *         g(k1,k4,k5,k2,k3,k6)+g(k1,k5,k4,k3,k2,k6)
         return
      endif
      end
****************************************************************
      subroutine saabb(nombre,n,mt,vex3,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie AABB
c     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*4 nombre
      include 'dimensions.inc'
      common/ipot43/ipa43(nexp43,7),nfas43(igrad4)
      common sm(nexp)
      dimension vex3(3)
      is=nfas43(mt)

      write(7,1) nombre,n,rms*627.51d0,emax*627.51d0,
     -is,mt-2,(vex3(i),i=1,3)
      do 2 k=1,is
         write(7,3) k,sm(k),k,ipa43(k,1),k,ipa43(k,2),k,ipa43(k,3),
     -         k,ipa43(k,4),k,ipa43(k,5),k,ipa43(k,6),k,ipa43(k,7)
2     continue
      write(7,10)
      write(7,11) 
      write(7,12) 
      write(7,13) mt-2
1     format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 6D PES      ',/,
     -'*     for the AABB system class fitted to ',i4,' points      ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_2=',i2,')                      ',/,
     -'      dimension vex3(3)                                      ',/,
     -'      dimension der(6)                                       ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),',
     ,'i7(is),cf(is)',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),   ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)    ',/,
     -'      data vex3(1)/',d23.15,'/                            ',/,
     -'      data vex3(2)/',d23.15,'/                            ',/,
     -'      data vex3(3)/',d23.15,'/                            ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/         ',/,
     -'      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/            ',/,
     -'      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            ')

3     format(
     -'      data cf(',i3,')/',d23.15,'/                           ',/,
     -'      data i1(',i3,')/',i2,'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,
     -'/,i4(',i3,')/',i2,'/                                        ',/,
     -'      data i5(',i3,')/',i2,'/,i6(',i3,')/',i2,'/,',
     -'i7(',i3,')/',i2,'/')
10    format(
     -'      ener=0.d0          ',/,
     -'      der12=0.d0         ',/,
     -'      der13=0.d0         ',/,
     -'      der14=0.d0         ',/,
     -'      der23=0.d0         ',/,
     -'      der24=0.d0         ',/,
     -'      der34=0.d0         ',/,
     -'      pux12= vex3(1)*r12 ',/,
     -'      pux13= vex3(2)*r13 ',/,
     -'      pux14= vex3(2)*r14 ',/,
     -'      pux23= vex3(2)*r23 ',/,
     -'      pux24= vex3(2)*r24 ',/,
     -'      pux34= vex3(3)*r34 ',/,
     -'      qux12= dexp(-pux12) ',/,
     -'      qux13= dexp(-pux13) ',/,
     -'      qux14= dexp(-pux14) ',/,
     -'      qux23= dexp(-pux23) ',/,
     -'      qux24= dexp(-pux24) ',/,
     -'      qux34= dexp(-pux34) ',/,
     -'      aux12= r12*qux12 ',/,
     -'      aux13= r13*qux13 ',/,
     -'      aux14= r14*qux14 ',/,
     -'      aux23= r23*qux23 ',/,
     -'      aux24= r24*qux24 ',/,
     -'      aux34= r34*qux34 ')

11    format(
     -'      do 1 i=1,mt_2               ',/,
     -'         f12(i) = aux12*f12(i-1)  ',/,
     -'         f13(i) = aux13*f13(i-1)  ',/,
     -'         f14(i) = aux14*f14(i-1)  ',/,
     -'         f23(i) = aux23*f23(i-1)  ',/,
     -'         f24(i) = aux24*f24(i-1)  ',/,
     -'         f34(i) = aux34*f34(i-1)  ',/,
     -'1        continue                 ')
12    format(
     -'      do 2 j=1,is     ',/,
     -'         k1=i1(j)     ',/,
     -'         k2=i2(j)     ',/,
     -'         k3=i3(j)     ',/,
     -'         k4=i4(j)     ',/,
     -'         k5=i5(j)     ',/,
     -'         k6=i6(j)     ',/,
     -'         k7=i7(j)     ',/,
     -'         cux1=f12(k1)*f34(k6)                      ',/,
     -'         cux2=f13(k2)*f14(k3)*f23(k4)*f24(k5)      ',/,
     -'         cux3=f13(k4)*f14(k5)*f23(k2)*f24(k3)      ',/,
     -'         cux4=f13(k3)*f14(k2)*f23(k5)*f24(k4)      ',/,
     -'         cux5=f13(k5)*f14(k4)*f23(k3)*f24(k2)      ',/,
     -'         if (k7.eq.1) faabb=cux1*cux2              ',/,
     -'         if (k7.eq.2) faabb=cux1*(cux2+cux3)       ',/,
     -'         if (k7.eq.3) faabb=cux1*(cux2+cux4+cux3+cux5) ',/,
     -'         ener=ener+cf(j)*faabb                         ',/,

     -'         if (iop.eq.1) then                        ',/,
     -'            d12=k1*faabb                       ',/,
     -'            d34=k6*faabb                       ',/,
     -'            if (k7.eq.1) then                      ',/,
     -'               d13=k2*faabb                       ',/,
     -'               d14=k3*faabb                       ',/,
     -'               d23=k4*faabb                       ',/,
     -'               d24=k5*faabb                       ',/,
     -'             else if (k7.eq.2) then               ',/,
     -'               d13=cux1*(k2*cux2+k4*cux3)         ',/,
     -'               d14=cux1*(k3*cux2+k5*cux3)         ',/,
     -'               d23=cux1*(k4*cux2+k2*cux3)         ',/,
     -'               d24=cux1*(k5*cux2+k3*cux3)         ',/,
     -'             else if (k7.eq.3) then                ',/,
     -'               d13=cux1*(k2*cux2+k4*cux3+k3*cux4+k5*cux5) ',/,
     -'               d14=cux1*(k3*cux2+k5*cux3+k2*cux4+k4*cux5) ',/,
     -'               d23=cux1*(k4*cux2+k2*cux3+k5*cux4+k3*cux5) ',/,
     -'               d24=cux1*(k5*cux2+k3*cux3+k4*cux4+k2*cux5) ',/,
     -'            endif                                        ',/,
     -'            d12=d12/f12(1)                       ',/,
     -'            d13=d13/f13(1)                       ',/,
     -'            d14=d14/f14(1)                       ',/,
     -'            d23=d23/f23(1)                       ',/,
     -'            d24=d24/f24(1)                       ',/,
     -'            d34=d34/f34(1)                       ',/,
     -'            der12=der12+cf(j)*d12                          ',/,
     -'            der13=der13+cf(j)*d13                          ',/,
     -'            der14=der14+cf(j)*d14                          ',/,
     -'            der23=der23+cf(j)*d23                          ',/,
     -'            der24=der24+cf(j)*d24                          ',/,
     -'            der34=der34+cf(j)*d34                          ',/,
     -'         endif                                             ',/,
     -'2     continue                                              ',/,
     -'      if (iop.eq.1) then                 ',/,
     -'         der(1)=der12*(1.d0-pux12)*qux12    ',/,
     -'         der(2)=der13*(1.d0-pux13)*qux13    ',/,
     -'         der(3)=der14*(1.d0-pux14)*qux14    ',/,
     -'         der(4)=der23*(1.d0-pux23)*qux23    ',/,
     -'         der(5)=der24*(1.d0-pux24)*qux24    ',/,
     -'         der(6)=der34*(1.d0-pux34)*qux34    ',/,
     -'      endif                                        ',/,
     -'      return                                                ',/,
     -'      end                                                   ')
**********************************************************************
13    format(
     -'*************************************************************',/,
     -'      function faabb(k1,k2,k3,k4,k5,k6,jt)                   ',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(mt_2=',i2,')                      ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
     -',/,
     -'      if (jt.eq.1) then                                      ',/,
     -'         faabb=g(k1,k2,k3,k4,k5,k6)                          ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      if (jt.eq.2) then                                      ',/,
     -'         faabb=g(k1,k2,k3,k4,k5,k6)+g(k1,k4,k5,k2,k3,k6)     ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      if (jt.eq.3) then                                      ',/,
     -'         faabb=g(k1,k2,k3,k4,k5,k6)+g(k1,k3,k2,k5,k4,k6)+    ',/,
     -'     *         g(k1,k4,k5,k2,k3,k6)+g(k1,k5,k4,k3,k2,k6)     ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      end                                                    ')
      end
**********************************************************************
      subroutine gfabcc(n,mt,vex4)
**********************************************************************
*     devuelve el valor de la fas j para cada punto i en xx(j,i)
**********************************************************************
      implicit real*8(a-h,o-z)
      include 'dimensions.inc'
      common/ipot44/ipa44(nexp44,7),nfas44(igrad4)
      common sm(nexp),xx(nexp,nmax),xener(nmax),sg(nexp*nexp),
     ,       r12(nmax),r13(nmax),r23(nmax),
     ,       r14(nmax),r24(nmax),r34(nmax)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      dimension vex4(4)
      do 1 i=1,n
         f12(0) = 1.d0
         f13(0) = 1.d0
         f14(0) = 1.d0
         f23(0) = 1.d0
         f24(0) = 1.d0
         f34(0) = 1.d0
c        vex4(1)--> AB
c        vex4(2)--> AC
c        vex4(3)--> BC
c        vex4(4)--> CC
         aux12= r12(i)*dexp(-vex4(1)*r12(i))
         aux13= r13(i)*dexp(-vex4(2)*r13(i))
         aux14= r14(i)*dexp(-vex4(2)*r14(i))
         aux23= r23(i)*dexp(-vex4(3)*r23(i))
         aux24= r24(i)*dexp(-vex4(3)*r24(i))
         aux34= r34(i)*dexp(-vex4(4)*r34(i))
         do 2 k=1,mt-2
            f12(k) = aux12*f12(k-1)
            f13(k) = aux13*f13(k-1)
            f14(k) = aux14*f14(k-1)
            f23(k) = aux23*f23(k-1)
            f24(k) = aux24*f24(k-1)
            f34(k) = aux34*f34(k-1)
2        continue
         do 3 j=1,nfas44(mt)
            i1=ipa44(j,1)
            i2=ipa44(j,2)
            i3=ipa44(j,3)
            i4=ipa44(j,4)
            i5=ipa44(j,5)
            i6=ipa44(j,6)
            i7=ipa44(j,7)
            xx(j,i)=fabcc(i1,i2,i3,i4,i5,i6,i7)
 3       continue
 1    continue
      end
**********************************************************************
      function fabcc(k1,k2,k3,k4,k5,k6,jt)
      implicit real*8(a-h,o-z)
      parameter (igrad4=12)
      common/basu4/f12(0:igrad4-2),f13(0:igrad4-2),f14(0:igrad4-2),
     ,             f23(0:igrad4-2),f24(0:igrad4-2),f34(0:igrad4-2)
      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
      if (jt.eq.1) then
         fabcc=g(k1,k2,k3,k4,k5,k6)
         return
      endif
      if (jt.eq.2) then
         fabcc=g(k1,k2,k3,k4,k5,k6)+g(k1,k3,k2,k5,k4,k6)
         return
      endif
      end
****************************************************************
      subroutine sabcc(nombre,n,mt,vex4,rms,emax)
****************************************************************
*     Escribe el programa fortran que calcula la superficie ABCC
c     necesita la matriz sm
****************************************************************
      implicit  real * 8 (a-h,o-z)
      character*4 nombre
      include 'dimensions.inc'
      common/ipot44/ipa44(nexp44,7),nfas44(igrad4)
      common sm(nexp)
      dimension vex4(4)
      is=nfas44(mt)

      write(7,1) nombre,n,rms*627.51d0,emax*627.51d0,
     -is,mt-2,(vex4(i),i=1,4)
      do 2 k=1,is
         write(7,3) k,sm(k),k,ipa44(k,1),k,ipa44(k,2),k,ipa44(k,3),
     -         k,ipa44(k,4),k,ipa44(k,5),k,ipa44(k,6),k,ipa44(k,7)
2     continue
      write(7,10)
      write(7,11) 
      write(7,12) 
      write(7,13) mt-2
1     format(
     -'*************************************************************',/,
     -'      subroutine ',a6,'(r12,r13,r14,r23,r24,r34,ener,der,iop)',/,
     -'*************************************************************',/,
     -'*     This subroutine computes the energies of a 6D PES      ',/,
     -'*     for the ABCC system class fitted to ',i4,' points      ',/,
     -'*     rms = ',f15.8,' kcal/mol                               ',/,
     -'*     emax= ',f15.8,' kcal/mol                               ',/,
     -'*************************************************************',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(is=',i5,', mt_2=',i2,')                      ',/,
     -'      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),',
     ,'i7(is),cf(is)',/,
     -'      dimension vex4(4)                                      ',/,
     -'      dimension der(6)                                       ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),   ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)    ',/,
     -'      data vex4(1)/',d23.15,'/                            ',/,
     -'      data vex4(2)/',d23.15,'/                            ',/,
     -'      data vex4(3)/',d23.15,'/                            ',/,
     -'      data vex4(4)/',d23.15,'/                            ',/,
     -'      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/         ',/,
     -'      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/         ',/,
     -'      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/            ',/,
     -'      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            ')
3     format(
     -'      data cf(',i3,')/',d23.15,'/                           ',/,
     -'      data i1(',i3,')/',i2,'/,i2(',i3,')/',i2,'/,i3(',i3,')/',i2,
     -'/,i4(',i3,')/',i2,'/                                        ',/,
     -'      data i5(',i3,')/',i2,'/,i6(',i3,')/',i2,'/,',
     -'i7(',i3,')/',i2,'/')
10    format(
     -'      ener=0.d0      ',/,
     -'      der12=0.d0         ',/,
     -'      der13=0.d0         ',/,
     -'      der14=0.d0         ',/,
     -'      der23=0.d0         ',/,
     -'      der24=0.d0         ',/,
     -'      der34=0.d0         ',/,
     -'      pux12= vex4(1)*r12 ',/,
     -'      pux13= vex4(2)*r13 ',/,
     -'      pux14= vex4(2)*r14 ',/,
     -'      pux23= vex4(3)*r23 ',/,
     -'      pux24= vex4(3)*r24 ',/,
     -'      pux34= vex4(4)*r34 ',/,
     -'      qux12= dexp(-pux12) ',/,
     -'      qux13= dexp(-pux13) ',/,
     -'      qux14= dexp(-pux14) ',/,
     -'      qux23= dexp(-pux23) ',/,
     -'      qux24= dexp(-pux24) ',/,
     -'      qux34= dexp(-pux34) ',/,
     -'      aux12= r12*qux12 ',/,
     -'      aux13= r13*qux13 ',/,
     -'      aux14= r14*qux14 ',/,
     -'      aux23= r23*qux23 ',/,
     -'      aux24= r24*qux24 ',/,
     -'      aux34= r34*qux34 ')

11    format(
     -'      do 1 i=1,mt_2               ',/,
     -'         f12(i) = aux12*f12(i-1)  ',/,
     -'         f13(i) = aux13*f13(i-1)  ',/,
     -'         f14(i) = aux14*f14(i-1)  ',/,
     -'         f23(i) = aux23*f23(i-1)  ',/,
     -'         f24(i) = aux24*f24(i-1)  ',/,
     -'         f34(i) = aux34*f34(i-1)  ',/,
     -'1        continue                 ')
12    format(
     -'      do 2 j=1,is     ',/,
     -'         k1=i1(j)     ',/,
     -'         k2=i2(j)     ',/,
     -'         k3=i3(j)     ',/,
     -'         k4=i4(j)     ',/,
     -'         k5=i5(j)     ',/,
     -'         k6=i6(j)     ',/,
     -'         k7=i7(j)     ',/,
     -'         cux1=f12(k1)*f34(k6)                      ',/,
     -'         cux2=f13(k2)*f14(k3)*f23(k4)*f24(k5)      ',/,
     -'         cux3=f13(k3)*f14(k2)*f23(k5)*f24(k4)      ',/,
     -'      if (k7.eq.1) then                            ',/,
     -'         fabcc=cux1*cux2                           ',/,
     -'       else                                        ',/,
     -'         fabcc=cux1*(cux2+cux3)                    ',/,
     -'      endif                                        ',/,
     -'         ener=ener+cf(j)*fabcc                     ',/,
     -'         if (iop.eq.1) then                        ',/,
     -'            d12=k1*fabcc                       ',/,
     -'            d34=k6*fabcc                       ',/,
     -'            if (k7.eq.1) then                     ',/,
     -'               d13=k2*fabcc                       ',/,
     -'               d14=k3*fabcc                       ',/,
     -'               d23=k4*fabcc                       ',/,
     -'               d24=k5*fabcc                       ',/,
     -'             else                                 ',/,
     -'               d13=cux1*(k2*cux2+k3*cux3)         ',/,
     -'               d14=cux1*(k3*cux2+k2*cux3)         ',/,
     -'               d23=cux1*(k4*cux2+k5*cux3)         ',/,
     -'               d24=cux1*(k5*cux2+k4*cux3)         ',/,
     -'            endif                                    ',/,
     -'            d12=d12/f12(1)                       ',/,
     -'            d13=d13/f13(1)                       ',/,
     -'            d14=d14/f14(1)                       ',/,
     -'            d23=d23/f23(1)                       ',/,
     -'            d24=d24/f24(1)                       ',/,
     -'            d34=d34/f34(1)                       ',/,
     -'            der12=der12+cf(j)*d12                          ',/,
     -'            der13=der13+cf(j)*d13                          ',/,
     -'            der14=der14+cf(j)*d14                          ',/,
     -'            der23=der23+cf(j)*d23                          ',/,
     -'            der24=der24+cf(j)*d24                          ',/,
     -'            der34=der34+cf(j)*d34                          ',/,
     -'         endif                                             ',/,
     -'2     continue                                              ',/,
     -'      if (iop.eq.1) then                 ',/,
     -'         der(1)=der12*(1.d0-pux12)*qux12    ',/,
     -'         der(2)=der13*(1.d0-pux13)*qux13    ',/,
     -'         der(3)=der14*(1.d0-pux14)*qux14    ',/,
     -'         der(4)=der23*(1.d0-pux23)*qux23    ',/,
     -'         der(5)=der24*(1.d0-pux24)*qux24    ',/,
     -'         der(6)=der34*(1.d0-pux34)*qux34    ',/,
     -'      endif                                        ',/,
     -'      return                                                ',/,
     -'      end                                                   ')
**********************************************************************
13    format(
     -'*************************************************************',/,
     -'      function fabcc(k1,k2,k3,k4,k5,k6,jt)                   ',/,
     -'      implicit real*8(a-h,o-z)                               ',/,
     -'      parameter(mt_2=',i2,')                      ',/,
     -'      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    ',/,
     -'     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     ',/,
     -'      g(i,j,k,l,m,n)=f12(i)*f13(j)*f14(k)*f23(l)*f24(m)*f34(n)
     -',/,
     -'      if (jt.eq.1) then                                      ',/,
     -'         fabcc=g(k1,k2,k3,k4,k5,k6)                          ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      if (jt.eq.2) then                                      ',/,
     -'         fabcc=g(k1,k2,k3,k4,k5,k6)+g(k1,k3,k2,k5,k4,k6)     ',/,
     -'         return                                              ',/,
     -'      endif                                                  ',/,
     -'      end                                                    ')
      end
      blockdata BB
      include 'dimensions.inc'
      common/ipot41/ipa41(nexp41,7),nfas41(igrad4)
      common/ipot42/ipa42(nexp42,7),nfas42(igrad4)
      common/ipot43/ipa43(nexp43,7),nfas43(igrad4)
      common/ipot44/ipa44(nexp44,7),nfas44(igrad4)
      common/ipot45/ipa45(nexp45,6),nfas45(igrad4)
*********** 0 0*******************************
*********** 1 0*******************************
*********** 2 0*******************************
      data (ipa41(   1,j),j=1,7)/  0, 0, 1, 0, 1, 1,24/
      data (ipa41(   2,j),j=1,7)/  0, 1, 0, 1, 1, 0,11/
*********** 3 2*******************************
      data (ipa41(   3,j),j=1,7)/  0, 0, 1, 0, 1, 2,24/
      data (ipa41(   4,j),j=1,7)/  0, 0, 1, 1, 0, 2,24/
      data (ipa41(   5,j),j=1,7)/  1, 1, 1, 1, 0, 0,13/
      data (ipa41(   6,j),j=1,7)/  0, 0, 1, 2, 0, 1,24/
      data (ipa41(   7,j),j=1,7)/  1, 1, 0, 0, 1, 1, 3/
*********** 4 5*******************************
      data (ipa41(   8,j),j=1,7)/  0, 0, 1, 0, 1, 3,24/
      data (ipa41(   9,j),j=1,7)/  0, 0, 1, 0, 2, 2,24/
      data (ipa41(  10,j),j=1,7)/  0, 0, 1, 1, 0, 3,24/
      data (ipa41(  11,j),j=1,7)/  0, 0, 1, 1, 1, 2,24/
      data (ipa41(  12,j),j=1,7)/  0, 0, 1, 2, 0, 2,24/
      data (ipa41(  13,j),j=1,7)/  0, 0, 1, 2, 1, 1,24/
      data (ipa41(  14,j),j=1,7)/  0, 0, 1, 3, 0, 1,24/
      data (ipa41(  15,j),j=1,7)/  0, 0, 2, 1, 1, 1,24/
      data (ipa41(  16,j),j=1,7)/  0, 0, 2, 2, 0, 1,24/
      data (ipa41(  17,j),j=1,7)/  1, 1, 0, 1, 1, 1, 6/
      data (ipa41(  18,j),j=1,7)/  0, 1, 1, 1, 2, 0,24/
*********** 5 11*******************************
      data (ipa41(  19,j),j=1,7)/  0, 0, 1, 0, 1, 4,24/
      data (ipa41(  20,j),j=1,7)/  0, 0, 1, 0, 2, 3,24/
      data (ipa41(  21,j),j=1,7)/  0, 0, 1, 1, 0, 4,24/
      data (ipa41(  22,j),j=1,7)/  0, 0, 1, 1, 1, 3,24/
      data (ipa41(  23,j),j=1,7)/  0, 0, 1, 1, 2, 2,24/
      data (ipa41(  24,j),j=1,7)/  0, 0, 1, 2, 0, 3,24/
      data (ipa41(  25,j),j=1,7)/  0, 0, 1, 2, 1, 2,24/
      data (ipa41(  26,j),j=1,7)/  0, 0, 1, 3, 0, 2,24/
      data (ipa41(  27,j),j=1,7)/  0, 0, 1, 3, 1, 1,24/
      data (ipa41(  28,j),j=1,7)/  0, 0, 1, 4, 0, 1,24/
      data (ipa41(  29,j),j=1,7)/  0, 0, 2, 0, 2, 2,24/
      data (ipa41(  30,j),j=1,7)/  0, 0, 2, 1, 1, 2,24/
      data (ipa41(  31,j),j=1,7)/  0, 2, 0, 2, 2, 0,11/
      data (ipa41(  32,j),j=1,7)/  0, 0, 2, 2, 1, 1,24/
      data (ipa41(  33,j),j=1,7)/  0, 0, 2, 3, 0, 1,24/
      data (ipa41(  34,j),j=1,7)/  0, 0, 3, 1, 1, 1,24/
      data (ipa41(  35,j),j=1,7)/  0, 1, 1, 1, 1, 2,24/
      data (ipa41(  36,j),j=1,7)/  0, 1, 1, 1, 2, 1,24/
      data (ipa41(  37,j),j=1,7)/  0, 1, 1, 1, 3, 0,24/
      data (ipa41(  38,j),j=1,7)/  0, 1, 1, 2, 2, 0,24/
      data (ipa41(  39,j),j=1,7)/  0, 1, 2, 2, 1, 0, 6/
      data (ipa41(  40,j),j=1,7)/  1, 1, 1, 1, 1, 1, 1/
*********** 6 22*******************************
      data (ipa41(  41,j),j=1,7)/  0, 0, 1, 0, 1, 5,24/
      data (ipa41(  42,j),j=1,7)/  0, 0, 1, 0, 2, 4,24/
      data (ipa41(  43,j),j=1,7)/  0, 0, 1, 0, 3, 3,24/
      data (ipa41(  44,j),j=1,7)/  0, 0, 1, 1, 0, 5,24/
      data (ipa41(  45,j),j=1,7)/  0, 0, 1, 1, 1, 4,24/
      data (ipa41(  46,j),j=1,7)/  0, 0, 1, 1, 2, 3,24/
      data (ipa41(  47,j),j=1,7)/  0, 0, 1, 2, 0, 4,24/
      data (ipa41(  48,j),j=1,7)/  0, 0, 1, 2, 1, 3,24/
      data (ipa41(  49,j),j=1,7)/  0, 0, 1, 2, 2, 2,24/
      data (ipa41(  50,j),j=1,7)/  0, 0, 1, 3, 0, 3,24/
      data (ipa41(  51,j),j=1,7)/  0, 0, 1, 3, 1, 2,24/
      data (ipa41(  52,j),j=1,7)/  0, 0, 1, 4, 0, 2,24/
      data (ipa41(  53,j),j=1,7)/  0, 0, 1, 4, 1, 1,24/
      data (ipa41(  54,j),j=1,7)/  0, 0, 1, 5, 0, 1,24/
      data (ipa41(  55,j),j=1,7)/  0, 0, 2, 0, 2, 3,24/
      data (ipa41(  56,j),j=1,7)/  0, 0, 2, 1, 1, 3,24/
      data (ipa41(  57,j),j=1,7)/  0, 0, 2, 1, 2, 2,24/
      data (ipa41(  58,j),j=1,7)/  0, 0, 2, 2, 0, 3,24/
      data (ipa41(  59,j),j=1,7)/  0, 0, 2, 2, 1, 2,24/
      data (ipa41(  60,j),j=1,7)/  0, 0, 2, 3, 0, 2,24/
      data (ipa41(  61,j),j=1,7)/  0, 0, 2, 3, 1, 1,13/
      data (ipa41(  62,j),j=1,7)/  0, 0, 2, 4, 0, 1,24/
      data (ipa41(  63,j),j=1,7)/  0, 0, 3, 1, 1, 2,24/
      data (ipa41(  64,j),j=1,7)/  0, 0, 3, 2, 1, 1,13/
      data (ipa41(  65,j),j=1,7)/  0, 0, 3, 3, 0, 1,24/
      data (ipa41(  66,j),j=1,7)/  0, 0, 4, 1, 1, 1,24/
      data (ipa41(  67,j),j=1,7)/  0, 1, 1, 1, 1, 3,24/
      data (ipa41(  68,j),j=1,7)/  0, 1, 1, 1, 2, 2,24/
      data (ipa41(  69,j),j=1,7)/  0, 1, 1, 1, 3, 1,24/
      data (ipa41(  70,j),j=1,7)/  0, 1, 1, 1, 4, 0,24/
      data (ipa41(  71,j),j=1,7)/  0, 1, 1, 2, 2, 1,24/
      data (ipa41(  72,j),j=1,7)/  0, 1, 1, 2, 3, 0,24/
      data (ipa41(  73,j),j=1,7)/  0, 1, 2, 1, 2, 1,24/
      data (ipa41(  74,j),j=1,7)/  0, 1, 2, 2, 1, 1,24/
      data (ipa41(  75,j),j=1,7)/  0, 1, 2, 2, 2, 0,24/
      data (ipa41(  76,j),j=1,7)/  0, 1, 2, 3, 1, 0,11/
      data (ipa41(  77,j),j=1,7)/  1, 1, 2, 1, 1, 1, 6/
*********** 7 37*******************************
      data (ipa41(  78,j),j=1,7)/  0, 0, 1, 0, 1, 6,24/
      data (ipa41(  79,j),j=1,7)/  0, 0, 1, 0, 2, 5,24/
      data (ipa41(  80,j),j=1,7)/  0, 0, 1, 0, 3, 4,24/
      data (ipa41(  81,j),j=1,7)/  0, 0, 1, 1, 0, 6,24/
      data (ipa41(  82,j),j=1,7)/  0, 0, 1, 1, 1, 5,24/
      data (ipa41(  83,j),j=1,7)/  0, 0, 1, 1, 2, 4,24/
      data (ipa41(  84,j),j=1,7)/  0, 0, 1, 1, 3, 3,24/
      data (ipa41(  85,j),j=1,7)/  0, 0, 1, 2, 0, 5,24/
      data (ipa41(  86,j),j=1,7)/  0, 0, 1, 2, 1, 4,24/
      data (ipa41(  87,j),j=1,7)/  0, 0, 1, 2, 2, 3,24/
      data (ipa41(  88,j),j=1,7)/  0, 0, 1, 3, 0, 4,24/
      data (ipa41(  89,j),j=1,7)/  0, 0, 1, 3, 1, 3,24/
      data (ipa41(  90,j),j=1,7)/  0, 0, 1, 3, 2, 2,13/
      data (ipa41(  91,j),j=1,7)/  0, 0, 1, 4, 0, 3,24/
      data (ipa41(  92,j),j=1,7)/  0, 0, 1, 4, 1, 2,24/
      data (ipa41(  93,j),j=1,7)/  0, 0, 1, 5, 0, 2,24/
      data (ipa41(  94,j),j=1,7)/  0, 0, 1, 5, 1, 1,24/
      data (ipa41(  95,j),j=1,7)/  0, 0, 1, 6, 0, 1,24/
      data (ipa41(  96,j),j=1,7)/  0, 0, 2, 0, 2, 4,24/
      data (ipa41(  97,j),j=1,7)/  0, 0, 2, 0, 3, 3,24/
      data (ipa41(  98,j),j=1,7)/  0, 0, 2, 1, 1, 4,24/
      data (ipa41(  99,j),j=1,7)/  0, 0, 2, 1, 2, 3,24/
      data (ipa41( 100,j),j=1,7)/  0, 0, 2, 2, 0, 4,24/
      data (ipa41( 101,j),j=1,7)/  0, 0, 2, 2, 1, 3,24/
      data (ipa41( 102,j),j=1,7)/  2, 2, 2, 2, 0, 0,13/
      data (ipa41( 103,j),j=1,7)/  0, 0, 2, 3, 0, 3,24/
      data (ipa41( 104,j),j=1,7)/  0, 0, 2, 3, 1, 2,24/
      data (ipa41( 105,j),j=1,7)/  0, 0, 2, 4, 0, 2,24/
      data (ipa41( 106,j),j=1,7)/  0, 0, 2, 4, 1, 1,13/
      data (ipa41( 107,j),j=1,7)/  0, 0, 2, 5, 0, 1,24/
      data (ipa41( 108,j),j=1,7)/  0, 0, 3, 1, 1, 3,24/
      data (ipa41( 109,j),j=1,7)/  0, 0, 3, 1, 2, 2,13/
      data (ipa41( 110,j),j=1,7)/  0, 0, 3, 2, 1, 2,24/
      data (ipa41( 111,j),j=1,7)/  0, 0, 3, 3, 0, 2,24/
      data (ipa41( 112,j),j=1,7)/  0, 0, 3, 3, 1, 1,24/
      data (ipa41( 113,j),j=1,7)/  0, 0, 3, 4, 0, 1,24/
      data (ipa41( 114,j),j=1,7)/  0, 0, 4, 1, 1, 2,24/
      data (ipa41( 115,j),j=1,7)/  0, 0, 4, 2, 1, 1,13/
      data (ipa41( 116,j),j=1,7)/  0, 0, 5, 1, 1, 1,24/
      data (ipa41( 117,j),j=1,7)/  0, 1, 1, 1, 1, 4,24/
      data (ipa41( 118,j),j=1,7)/  0, 1, 1, 1, 2, 3,24/
      data (ipa41( 119,j),j=1,7)/  0, 1, 1, 1, 3, 2,24/
      data (ipa41( 120,j),j=1,7)/  0, 1, 1, 1, 4, 1,24/
      data (ipa41( 121,j),j=1,7)/  0, 1, 1, 1, 5, 0,24/
      data (ipa41( 122,j),j=1,7)/  0, 1, 1, 2, 2, 2,24/
      data (ipa41( 123,j),j=1,7)/  0, 1, 1, 2, 3, 1,24/
      data (ipa41( 124,j),j=1,7)/  0, 1, 1, 2, 4, 0,24/
      data (ipa41( 125,j),j=1,7)/  0, 1, 1, 3, 3, 0,24/
      data (ipa41( 126,j),j=1,7)/  0, 1, 2, 1, 2, 2,24/
      data (ipa41( 127,j),j=1,7)/  0, 1, 2, 1, 3, 1,24/
      data (ipa41( 128,j),j=1,7)/  0, 1, 2, 2, 1, 2,24/
      data (ipa41( 129,j),j=1,7)/  0, 1, 2, 2, 2, 1,24/
      data (ipa41( 130,j),j=1,7)/  0, 2, 1, 3, 2, 0,11/
      data (ipa41( 131,j),j=1,7)/  0, 1, 2, 3, 1, 1,24/
      data (ipa41( 132,j),j=1,7)/  0, 1, 2, 3, 2, 0,24/
      data (ipa41( 133,j),j=1,7)/  0, 1, 2, 4, 1, 0,11/
      data (ipa41( 134,j),j=1,7)/  0, 1, 3, 3, 1, 0, 6/
      data (ipa41( 135,j),j=1,7)/  2, 2, 0, 0, 2, 2, 3/
      data (ipa41( 136,j),j=1,7)/  1, 1, 3, 1, 1, 1, 6/
      data (ipa41( 137,j),j=1,7)/  1, 1, 1, 1, 2, 2,13/
      data (ipa41( 138,j),j=1,7)/  1, 1, 2, 2, 1, 1, 3/
*********** 8 61*******************************
      data (ipa41( 139,j),j=1,7)/  0, 0, 1, 0, 1, 7,24/
      data (ipa41( 140,j),j=1,7)/  0, 0, 1, 0, 2, 6,24/
      data (ipa41( 141,j),j=1,7)/  0, 0, 1, 0, 3, 5,24/
      data (ipa41( 142,j),j=1,7)/  0, 0, 1, 0, 4, 4,24/
      data (ipa41( 143,j),j=1,7)/  0, 0, 1, 1, 0, 7,24/
      data (ipa41( 144,j),j=1,7)/  0, 0, 1, 1, 1, 6,24/
      data (ipa41( 145,j),j=1,7)/  0, 0, 1, 1, 2, 5,24/
      data (ipa41( 146,j),j=1,7)/  0, 0, 1, 1, 3, 4,24/
      data (ipa41( 147,j),j=1,7)/  0, 0, 1, 2, 0, 6,24/
      data (ipa41( 148,j),j=1,7)/  0, 0, 1, 2, 1, 5,24/
      data (ipa41( 149,j),j=1,7)/  0, 0, 1, 2, 2, 4,24/
      data (ipa41( 150,j),j=1,7)/  0, 0, 1, 2, 3, 3,13/
      data (ipa41( 151,j),j=1,7)/  0, 0, 1, 3, 0, 5,24/
      data (ipa41( 152,j),j=1,7)/  0, 0, 1, 3, 1, 4,24/
      data (ipa41( 153,j),j=1,7)/  0, 0, 1, 3, 2, 3,24/
      data (ipa41( 154,j),j=1,7)/  0, 0, 1, 4, 0, 4,24/
      data (ipa41( 155,j),j=1,7)/  0, 0, 1, 4, 1, 3,24/
      data (ipa41( 156,j),j=1,7)/  0, 0, 1, 4, 2, 2,13/
      data (ipa41( 157,j),j=1,7)/  0, 0, 1, 5, 0, 3,24/
      data (ipa41( 158,j),j=1,7)/  0, 0, 1, 5, 1, 2,24/
      data (ipa41( 159,j),j=1,7)/  0, 0, 1, 6, 0, 2,24/
      data (ipa41( 160,j),j=1,7)/  0, 0, 1, 6, 1, 1,24/
      data (ipa41( 161,j),j=1,7)/  0, 0, 1, 7, 0, 1,24/
      data (ipa41( 162,j),j=1,7)/  0, 0, 2, 0, 2, 5,24/
      data (ipa41( 163,j),j=1,7)/  0, 0, 2, 0, 3, 4,24/
      data (ipa41( 164,j),j=1,7)/  0, 0, 2, 1, 1, 5,24/
      data (ipa41( 165,j),j=1,7)/  0, 0, 2, 1, 2, 4,24/
      data (ipa41( 166,j),j=1,7)/  0, 0, 2, 1, 3, 3,13/
      data (ipa41( 167,j),j=1,7)/  0, 0, 2, 2, 0, 5,24/
      data (ipa41( 168,j),j=1,7)/  0, 0, 2, 2, 1, 4,24/
      data (ipa41( 169,j),j=1,7)/  0, 0, 2, 2, 2, 3,24/
      data (ipa41( 170,j),j=1,7)/  0, 0, 2, 3, 0, 4,24/
      data (ipa41( 171,j),j=1,7)/  0, 0, 2, 3, 1, 3,24/
      data (ipa41( 172,j),j=1,7)/  0, 0, 2, 3, 2, 2,24/
      data (ipa41( 173,j),j=1,7)/  0, 0, 2, 4, 0, 3,24/
      data (ipa41( 174,j),j=1,7)/  0, 0, 2, 4, 1, 2,24/
      data (ipa41( 175,j),j=1,7)/  0, 0, 2, 5, 0, 2,24/
      data (ipa41( 176,j),j=1,7)/  0, 0, 2, 5, 1, 1,13/
      data (ipa41( 177,j),j=1,7)/  0, 0, 2, 6, 0, 1,24/
      data (ipa41( 178,j),j=1,7)/  0, 0, 3, 0, 3, 3,24/
      data (ipa41( 179,j),j=1,7)/  0, 0, 3, 1, 1, 4,24/
      data (ipa41( 180,j),j=1,7)/  0, 0, 3, 1, 2, 3,24/
      data (ipa41( 181,j),j=1,7)/  0, 0, 3, 2, 1, 3,24/
      data (ipa41( 182,j),j=1,7)/  0, 0, 3, 2, 2, 2,24/
      data (ipa41( 183,j),j=1,7)/  0, 3, 0, 3, 3, 0,11/
      data (ipa41( 184,j),j=1,7)/  0, 0, 3, 3, 1, 2,24/
      data (ipa41( 185,j),j=1,7)/  0, 0, 3, 4, 0, 2,24/
      data (ipa41( 186,j),j=1,7)/  0, 0, 3, 4, 1, 1,13/
      data (ipa41( 187,j),j=1,7)/  0, 0, 3, 5, 0, 1,24/
      data (ipa41( 188,j),j=1,7)/  0, 0, 4, 1, 1, 3,24/
      data (ipa41( 189,j),j=1,7)/  0, 0, 4, 1, 2, 2,13/
      data (ipa41( 190,j),j=1,7)/  0, 0, 4, 2, 1, 2,24/
      data (ipa41( 191,j),j=1,7)/  0, 0, 4, 3, 1, 1,13/
      data (ipa41( 192,j),j=1,7)/  0, 0, 4, 4, 0, 1,24/
      data (ipa41( 193,j),j=1,7)/  0, 0, 5, 1, 1, 2,24/
      data (ipa41( 194,j),j=1,7)/  0, 0, 5, 2, 1, 1,13/
      data (ipa41( 195,j),j=1,7)/  0, 0, 6, 1, 1, 1,24/
      data (ipa41( 196,j),j=1,7)/  0, 1, 1, 1, 1, 5,24/
      data (ipa41( 197,j),j=1,7)/  0, 1, 1, 1, 2, 4,24/
      data (ipa41( 198,j),j=1,7)/  0, 1, 1, 1, 3, 3,24/
      data (ipa41( 199,j),j=1,7)/  0, 1, 1, 1, 4, 2,24/
      data (ipa41( 200,j),j=1,7)/  0, 1, 1, 1, 5, 1,24/
      data (ipa41( 201,j),j=1,7)/  0, 1, 1, 1, 6, 0,24/
      data (ipa41( 202,j),j=1,7)/  1, 1, 0, 3, 2, 2,13/
      data (ipa41( 203,j),j=1,7)/  0, 1, 1, 2, 3, 2,24/
      data (ipa41( 204,j),j=1,7)/  0, 1, 1, 2, 4, 1,24/
      data (ipa41( 205,j),j=1,7)/  0, 1, 1, 2, 5, 0,24/
      data (ipa41( 206,j),j=1,7)/  0, 1, 1, 3, 3, 1,24/
      data (ipa41( 207,j),j=1,7)/  0, 1, 1, 3, 4, 0,24/
      data (ipa41( 208,j),j=1,7)/  0, 1, 2, 1, 2, 3,24/
      data (ipa41( 209,j),j=1,7)/  0, 1, 2, 1, 3, 2,24/
      data (ipa41( 210,j),j=1,7)/  0, 1, 2, 1, 4, 1,24/
      data (ipa41( 211,j),j=1,7)/  1, 2, 0, 3, 2, 1,11/
      data (ipa41( 212,j),j=1,7)/  0, 1, 2, 2, 2, 2,24/
      data (ipa41( 213,j),j=1,7)/  0, 1, 2, 2, 3, 1,24/
      data (ipa41( 214,j),j=1,7)/  0, 2, 1, 4, 2, 0,11/
      data (ipa41( 215,j),j=1,7)/  0, 1, 2, 3, 1, 2,24/
      data (ipa41( 216,j),j=1,7)/  0, 1, 2, 3, 2, 1,24/
      data (ipa41( 217,j),j=1,7)/  0, 1, 2, 3, 3, 0,24/
      data (ipa41( 218,j),j=1,7)/  0, 1, 2, 4, 1, 1,24/
      data (ipa41( 219,j),j=1,7)/  0, 1, 2, 4, 2, 0,24/
      data (ipa41( 220,j),j=1,7)/  0, 1, 2, 5, 1, 0,11/
      data (ipa41( 221,j),j=1,7)/  0, 1, 3, 1, 3, 1,24/
      data (ipa41( 222,j),j=1,7)/  0, 1, 3, 2, 2, 1,24/
      data (ipa41( 223,j),j=1,7)/  0, 1, 3, 3, 1, 1,24/
      data (ipa41( 224,j),j=1,7)/  0, 3, 1, 2, 3, 0,11/
      data (ipa41( 225,j),j=1,7)/  0, 1, 3, 4, 1, 0,11/
      data (ipa41( 226,j),j=1,7)/  0, 2, 2, 2, 2, 1,24/
      data (ipa41( 227,j),j=1,7)/  0, 2, 2, 2, 3, 0,24/
      data (ipa41( 228,j),j=1,7)/  1, 1, 4, 1, 1, 1, 6/
      data (ipa41( 229,j),j=1,7)/  1, 1, 1, 1, 2, 3,24/
      data (ipa41( 230,j),j=1,7)/  1, 1, 1, 2, 2, 2, 4/
      data (ipa41( 231,j),j=1,7)/  1, 1, 2, 1, 2, 2,24/
      data (ipa41( 232,j),j=1,7)/  1, 2, 1, 2, 2, 1,11/
      data (ipa41( 233,j),j=1,7)/  1, 1, 2, 3, 1, 1,24/
*********** 9 95*******************************
      data (ipa41( 234,j),j=1,7)/  0, 0, 1, 0, 1, 8,24/
      data (ipa41( 235,j),j=1,7)/  0, 0, 1, 0, 2, 7,24/
      data (ipa41( 236,j),j=1,7)/  0, 0, 1, 0, 3, 6,24/
      data (ipa41( 237,j),j=1,7)/  0, 0, 1, 0, 4, 5,24/
      data (ipa41( 238,j),j=1,7)/  0, 0, 1, 1, 0, 8,24/
      data (ipa41( 239,j),j=1,7)/  0, 0, 1, 1, 1, 7,24/
      data (ipa41( 240,j),j=1,7)/  0, 0, 1, 1, 2, 6,24/
      data (ipa41( 241,j),j=1,7)/  0, 0, 1, 1, 3, 5,24/
      data (ipa41( 242,j),j=1,7)/  0, 0, 1, 1, 4, 4,24/
      data (ipa41( 243,j),j=1,7)/  0, 0, 1, 2, 0, 7,24/
      data (ipa41( 244,j),j=1,7)/  0, 0, 1, 2, 1, 6,24/
      data (ipa41( 245,j),j=1,7)/  0, 0, 1, 2, 2, 5,24/
      data (ipa41( 246,j),j=1,7)/  0, 0, 1, 2, 3, 4,24/
      data (ipa41( 247,j),j=1,7)/  0, 0, 1, 3, 0, 6,24/
      data (ipa41( 248,j),j=1,7)/  0, 0, 1, 3, 1, 5,24/
      data (ipa41( 249,j),j=1,7)/  0, 0, 1, 3, 2, 4,24/
      data (ipa41( 250,j),j=1,7)/  0, 0, 1, 3, 3, 3,24/
      data (ipa41( 251,j),j=1,7)/  0, 0, 1, 4, 0, 5,24/
      data (ipa41( 252,j),j=1,7)/  0, 0, 1, 4, 1, 4,24/
      data (ipa41( 253,j),j=1,7)/  0, 0, 1, 4, 2, 3,24/
      data (ipa41( 254,j),j=1,7)/  0, 0, 1, 5, 0, 4,24/
      data (ipa41( 255,j),j=1,7)/  0, 0, 1, 5, 1, 3,24/
      data (ipa41( 256,j),j=1,7)/  0, 0, 1, 5, 2, 2,13/
      data (ipa41( 257,j),j=1,7)/  0, 0, 1, 6, 0, 3,24/
      data (ipa41( 258,j),j=1,7)/  0, 0, 1, 6, 1, 2,24/
      data (ipa41( 259,j),j=1,7)/  0, 0, 1, 7, 0, 2,24/
      data (ipa41( 260,j),j=1,7)/  0, 0, 1, 7, 1, 1,24/
      data (ipa41( 261,j),j=1,7)/  0, 0, 1, 8, 0, 1,24/
      data (ipa41( 262,j),j=1,7)/  0, 0, 2, 0, 2, 6,24/
      data (ipa41( 263,j),j=1,7)/  0, 0, 2, 0, 3, 5,24/
      data (ipa41( 264,j),j=1,7)/  0, 0, 2, 0, 4, 4,24/
      data (ipa41( 265,j),j=1,7)/  0, 0, 2, 1, 1, 6,24/
      data (ipa41( 266,j),j=1,7)/  0, 0, 2, 1, 2, 5,24/
      data (ipa41( 267,j),j=1,7)/  0, 0, 2, 1, 3, 4,24/
      data (ipa41( 268,j),j=1,7)/  0, 0, 2, 2, 0, 6,24/
      data (ipa41( 269,j),j=1,7)/  0, 0, 2, 2, 1, 5,24/
      data (ipa41( 270,j),j=1,7)/  0, 0, 2, 2, 2, 4,24/
      data (ipa41( 271,j),j=1,7)/  0, 0, 2, 2, 3, 3,24/
      data (ipa41( 272,j),j=1,7)/  0, 0, 2, 3, 0, 5,24/
      data (ipa41( 273,j),j=1,7)/  0, 0, 2, 3, 1, 4,24/
      data (ipa41( 274,j),j=1,7)/  0, 0, 2, 3, 2, 3,24/
      data (ipa41( 275,j),j=1,7)/  0, 0, 2, 4, 0, 4,24/
      data (ipa41( 276,j),j=1,7)/  0, 0, 2, 4, 1, 3,24/
      data (ipa41( 277,j),j=1,7)/  0, 0, 2, 4, 2, 2,24/
      data (ipa41( 278,j),j=1,7)/  0, 0, 2, 5, 0, 3,24/
      data (ipa41( 279,j),j=1,7)/  0, 0, 2, 5, 1, 2,24/
      data (ipa41( 280,j),j=1,7)/  0, 0, 2, 6, 0, 2,24/
      data (ipa41( 281,j),j=1,7)/  0, 0, 2, 6, 1, 1,13/
      data (ipa41( 282,j),j=1,7)/  0, 0, 2, 7, 0, 1,24/
      data (ipa41( 283,j),j=1,7)/  0, 0, 3, 0, 3, 4,24/
      data (ipa41( 284,j),j=1,7)/  0, 0, 3, 1, 1, 5,24/
      data (ipa41( 285,j),j=1,7)/  0, 0, 3, 1, 2, 4,24/
      data (ipa41( 286,j),j=1,7)/  0, 0, 3, 1, 3, 3,24/
      data (ipa41( 287,j),j=1,7)/  0, 0, 3, 2, 1, 4,24/
      data (ipa41( 288,j),j=1,7)/  0, 0, 3, 2, 2, 3,24/
      data (ipa41( 289,j),j=1,7)/  0, 0, 3, 3, 0, 4,24/
      data (ipa41( 290,j),j=1,7)/  0, 0, 3, 3, 1, 3,24/
      data (ipa41( 291,j),j=1,7)/  0, 0, 3, 3, 2, 2,24/
      data (ipa41( 292,j),j=1,7)/  0, 0, 3, 4, 0, 3,24/
      data (ipa41( 293,j),j=1,7)/  0, 0, 3, 4, 1, 2,24/
      data (ipa41( 294,j),j=1,7)/  0, 0, 3, 5, 0, 2,24/
      data (ipa41( 295,j),j=1,7)/  0, 0, 3, 5, 1, 1,13/
      data (ipa41( 296,j),j=1,7)/  0, 0, 3, 6, 0, 1,24/
      data (ipa41( 297,j),j=1,7)/  0, 0, 4, 1, 1, 4,24/
      data (ipa41( 298,j),j=1,7)/  0, 0, 4, 1, 2, 3,24/
      data (ipa41( 299,j),j=1,7)/  0, 0, 4, 2, 1, 3,24/
      data (ipa41( 300,j),j=1,7)/  0, 0, 4, 2, 2, 2,24/
      data (ipa41( 301,j),j=1,7)/  0, 0, 4, 3, 1, 2,24/
      data (ipa41( 302,j),j=1,7)/  0, 0, 4, 4, 0, 2,24/
      data (ipa41( 303,j),j=1,7)/  0, 0, 4, 4, 1, 1,24/
      data (ipa41( 304,j),j=1,7)/  0, 0, 4, 5, 0, 1,24/
      data (ipa41( 305,j),j=1,7)/  0, 0, 5, 1, 1, 3,24/
      data (ipa41( 306,j),j=1,7)/  0, 0, 5, 1, 2, 2,13/
      data (ipa41( 307,j),j=1,7)/  0, 0, 5, 2, 1, 2,24/
      data (ipa41( 308,j),j=1,7)/  0, 0, 5, 3, 1, 1,13/
      data (ipa41( 309,j),j=1,7)/  0, 0, 6, 1, 1, 2,24/
      data (ipa41( 310,j),j=1,7)/  0, 0, 6, 2, 1, 1,13/
      data (ipa41( 311,j),j=1,7)/  0, 0, 7, 1, 1, 1,24/
      data (ipa41( 312,j),j=1,7)/  0, 1, 1, 1, 1, 6,24/
      data (ipa41( 313,j),j=1,7)/  0, 1, 1, 1, 2, 5,24/
      data (ipa41( 314,j),j=1,7)/  0, 1, 1, 1, 3, 4,24/
      data (ipa41( 315,j),j=1,7)/  0, 1, 1, 1, 4, 3,24/
      data (ipa41( 316,j),j=1,7)/  0, 1, 1, 1, 5, 2,24/
      data (ipa41( 317,j),j=1,7)/  0, 1, 1, 1, 6, 1,24/
      data (ipa41( 318,j),j=1,7)/  0, 1, 1, 1, 7, 0,24/
      data (ipa41( 319,j),j=1,7)/  1, 1, 0, 4, 2, 2,13/
      data (ipa41( 320,j),j=1,7)/  0, 1, 1, 2, 3, 3,24/
      data (ipa41( 321,j),j=1,7)/  0, 1, 1, 2, 4, 2,24/
      data (ipa41( 322,j),j=1,7)/  0, 1, 1, 2, 5, 1,24/
      data (ipa41( 323,j),j=1,7)/  0, 1, 1, 2, 6, 0,24/
      data (ipa41( 324,j),j=1,7)/  1, 1, 0, 2, 3, 3,13/
      data (ipa41( 325,j),j=1,7)/  0, 1, 1, 3, 4, 1,24/
      data (ipa41( 326,j),j=1,7)/  0, 1, 1, 3, 5, 0,24/
      data (ipa41( 327,j),j=1,7)/  0, 1, 1, 4, 4, 0,24/
      data (ipa41( 328,j),j=1,7)/  0, 1, 2, 1, 2, 4,24/
      data (ipa41( 329,j),j=1,7)/  0, 1, 2, 1, 3, 3,24/
      data (ipa41( 330,j),j=1,7)/  0, 1, 2, 1, 4, 2,24/
      data (ipa41( 331,j),j=1,7)/  0, 1, 2, 1, 5, 1,24/
      data (ipa41( 332,j),j=1,7)/  1, 2, 0, 4, 2, 1,11/
      data (ipa41( 333,j),j=1,7)/  0, 1, 2, 2, 2, 3,24/
      data (ipa41( 334,j),j=1,7)/  0, 1, 2, 2, 3, 2,24/
      data (ipa41( 335,j),j=1,7)/  0, 1, 2, 2, 4, 1,24/
      data (ipa41( 336,j),j=1,7)/  0, 2, 1, 5, 2, 0,11/
      data (ipa41( 337,j),j=1,7)/  0, 1, 2, 3, 1, 3,24/
      data (ipa41( 338,j),j=1,7)/  0, 1, 2, 3, 2, 2,24/
      data (ipa41( 339,j),j=1,7)/  0, 1, 2, 3, 3, 1,24/
      data (ipa41( 340,j),j=1,7)/  0, 1, 2, 3, 4, 0,24/
      data (ipa41( 341,j),j=1,7)/  0, 1, 2, 4, 1, 2,24/
      data (ipa41( 342,j),j=1,7)/  0, 1, 2, 4, 2, 1,24/
      data (ipa41( 343,j),j=1,7)/  0, 1, 2, 4, 3, 0,24/
      data (ipa41( 344,j),j=1,7)/  0, 1, 2, 5, 1, 1,24/
      data (ipa41( 345,j),j=1,7)/  0, 1, 2, 5, 2, 0,24/
      data (ipa41( 346,j),j=1,7)/  0, 1, 2, 6, 1, 0,11/
      data (ipa41( 347,j),j=1,7)/  0, 1, 3, 1, 3, 2,24/
      data (ipa41( 348,j),j=1,7)/  0, 1, 3, 1, 4, 1,24/
      data (ipa41( 349,j),j=1,7)/  0, 1, 3, 2, 2, 2,24/
      data (ipa41( 350,j),j=1,7)/  0, 1, 3, 2, 3, 1,24/
      data (ipa41( 351,j),j=1,7)/  1, 3, 0, 2, 3, 1,11/
      data (ipa41( 352,j),j=1,7)/  0, 1, 3, 3, 2, 1,24/
      data (ipa41( 353,j),j=1,7)/  0, 1, 3, 3, 3, 0,24/
      data (ipa41( 354,j),j=1,7)/  0, 1, 3, 4, 1, 1,24/
      data (ipa41( 355,j),j=1,7)/  0, 1, 3, 4, 2, 0,24/
      data (ipa41( 356,j),j=1,7)/  0, 1, 3, 5, 1, 0,11/
      data (ipa41( 357,j),j=1,7)/  0, 1, 4, 2, 2, 1,24/
      data (ipa41( 358,j),j=1,7)/  0, 1, 4, 4, 1, 0, 6/
      data (ipa41( 359,j),j=1,7)/  2, 2, 0, 2, 2, 2, 6/
      data (ipa41( 360,j),j=1,7)/  0, 2, 2, 2, 3, 1,24/
      data (ipa41( 361,j),j=1,7)/  0, 2, 2, 2, 4, 0,24/
      data (ipa41( 362,j),j=1,7)/  0, 2, 2, 3, 3, 0,24/
      data (ipa41( 363,j),j=1,7)/  0, 2, 3, 3, 2, 0, 6/
      data (ipa41( 364,j),j=1,7)/  1, 1, 5, 1, 1, 1, 6/
      data (ipa41( 365,j),j=1,7)/  1, 1, 1, 1, 2, 4,24/
      data (ipa41( 366,j),j=1,7)/  1, 1, 1, 1, 3, 3,13/
      data (ipa41( 367,j),j=1,7)/  1, 1, 1, 2, 2, 3,24/
      data (ipa41( 368,j),j=1,7)/  1, 1, 2, 1, 2, 3,24/
      data (ipa41( 369,j),j=1,7)/  1, 1, 2, 2, 1, 3,24/
      data (ipa41( 370,j),j=1,7)/  2, 2, 2, 2, 1, 1,13/
      data (ipa41( 371,j),j=1,7)/  1, 1, 2, 3, 1, 2,24/
      data (ipa41( 372,j),j=1,7)/  1, 1, 2, 4, 1, 1,24/
      data (ipa41( 373,j),j=1,7)/  1, 1, 3, 3, 1, 1, 3/
      data (ipa41( 374,j),j=1,7)/  2, 2, 1, 1, 2, 2, 3/
*********** 10 141*******************************
      data (ipa41( 375,j),j=1,7)/  0, 0, 1, 0, 1, 9,24/
      data (ipa41( 376,j),j=1,7)/  0, 0, 1, 0, 2, 8,24/
      data (ipa41( 377,j),j=1,7)/  0, 0, 1, 0, 3, 7,24/
      data (ipa41( 378,j),j=1,7)/  0, 0, 1, 0, 4, 6,24/
      data (ipa41( 379,j),j=1,7)/  0, 0, 1, 0, 5, 5,24/
      data (ipa41( 380,j),j=1,7)/  0, 0, 1, 1, 0, 9,24/
      data (ipa41( 381,j),j=1,7)/  0, 0, 1, 1, 1, 8,24/
      data (ipa41( 382,j),j=1,7)/  0, 0, 1, 1, 2, 7,24/
      data (ipa41( 383,j),j=1,7)/  0, 0, 1, 1, 3, 6,24/
      data (ipa41( 384,j),j=1,7)/  0, 0, 1, 1, 4, 5,24/
      data (ipa41( 385,j),j=1,7)/  0, 0, 1, 2, 0, 8,24/
      data (ipa41( 386,j),j=1,7)/  0, 0, 1, 2, 1, 7,24/
      data (ipa41( 387,j),j=1,7)/  0, 0, 1, 2, 2, 6,24/
      data (ipa41( 388,j),j=1,7)/  0, 0, 1, 2, 3, 5,24/
      data (ipa41( 389,j),j=1,7)/  0, 0, 1, 2, 4, 4,13/
      data (ipa41( 390,j),j=1,7)/  0, 0, 1, 3, 0, 7,24/
      data (ipa41( 391,j),j=1,7)/  0, 0, 1, 3, 1, 6,24/
      data (ipa41( 392,j),j=1,7)/  0, 0, 1, 3, 2, 5,24/
      data (ipa41( 393,j),j=1,7)/  0, 0, 1, 3, 3, 4,24/
      data (ipa41( 394,j),j=1,7)/  0, 0, 1, 4, 0, 6,24/
      data (ipa41( 395,j),j=1,7)/  0, 0, 1, 4, 1, 5,24/
      data (ipa41( 396,j),j=1,7)/  0, 0, 1, 4, 2, 4,24/
      data (ipa41( 397,j),j=1,7)/  0, 0, 1, 4, 3, 3,13/
      data (ipa41( 398,j),j=1,7)/  0, 0, 1, 5, 0, 5,24/
      data (ipa41( 399,j),j=1,7)/  0, 0, 1, 5, 1, 4,24/
      data (ipa41( 400,j),j=1,7)/  0, 0, 1, 5, 2, 3,24/
      data (ipa41( 401,j),j=1,7)/  0, 0, 1, 6, 0, 4,24/
      data (ipa41( 402,j),j=1,7)/  0, 0, 1, 6, 1, 3,24/
      data (ipa41( 403,j),j=1,7)/  0, 0, 1, 6, 2, 2,13/
      data (ipa41( 404,j),j=1,7)/  0, 0, 1, 7, 0, 3,24/
      data (ipa41( 405,j),j=1,7)/  0, 0, 1, 7, 1, 2,24/
      data (ipa41( 406,j),j=1,7)/  0, 0, 1, 8, 0, 2,24/
      data (ipa41( 407,j),j=1,7)/  0, 0, 1, 8, 1, 1,24/
      data (ipa41( 408,j),j=1,7)/  0, 0, 1, 9, 0, 1,24/
      data (ipa41( 409,j),j=1,7)/  0, 0, 2, 0, 2, 7,24/
      data (ipa41( 410,j),j=1,7)/  0, 0, 2, 0, 3, 6,24/
      data (ipa41( 411,j),j=1,7)/  0, 0, 2, 0, 4, 5,24/
      data (ipa41( 412,j),j=1,7)/  0, 0, 2, 1, 1, 7,24/
      data (ipa41( 413,j),j=1,7)/  0, 0, 2, 1, 2, 6,24/
      data (ipa41( 414,j),j=1,7)/  0, 0, 2, 1, 3, 5,24/
      data (ipa41( 415,j),j=1,7)/  0, 0, 2, 1, 4, 4,13/
      data (ipa41( 416,j),j=1,7)/  0, 0, 2, 2, 0, 7,24/
      data (ipa41( 417,j),j=1,7)/  0, 0, 2, 2, 1, 6,24/
      data (ipa41( 418,j),j=1,7)/  0, 0, 2, 2, 2, 5,24/
      data (ipa41( 419,j),j=1,7)/  0, 0, 2, 2, 3, 4,24/
      data (ipa41( 420,j),j=1,7)/  0, 0, 2, 3, 0, 6,24/
      data (ipa41( 421,j),j=1,7)/  0, 0, 2, 3, 1, 5,24/
      data (ipa41( 422,j),j=1,7)/  0, 0, 2, 3, 2, 4,24/
      data (ipa41( 423,j),j=1,7)/  0, 0, 2, 3, 3, 3,24/
      data (ipa41( 424,j),j=1,7)/  0, 0, 2, 4, 0, 5,24/
      data (ipa41( 425,j),j=1,7)/  0, 0, 2, 4, 1, 4,24/
      data (ipa41( 426,j),j=1,7)/  0, 0, 2, 4, 2, 3,24/
      data (ipa41( 427,j),j=1,7)/  0, 0, 2, 5, 0, 4,24/
      data (ipa41( 428,j),j=1,7)/  0, 0, 2, 5, 1, 3,24/
      data (ipa41( 429,j),j=1,7)/  0, 0, 2, 5, 2, 2,24/
      data (ipa41( 430,j),j=1,7)/  0, 0, 2, 6, 0, 3,24/
      data (ipa41( 431,j),j=1,7)/  0, 0, 2, 6, 1, 2,24/
      data (ipa41( 432,j),j=1,7)/  0, 0, 2, 7, 0, 2,24/
      data (ipa41( 433,j),j=1,7)/  0, 0, 2, 7, 1, 1,13/
      data (ipa41( 434,j),j=1,7)/  0, 0, 2, 8, 0, 1,24/
      data (ipa41( 435,j),j=1,7)/  0, 0, 3, 0, 3, 5,24/
      data (ipa41( 436,j),j=1,7)/  0, 0, 3, 0, 4, 4,24/
      data (ipa41( 437,j),j=1,7)/  0, 0, 3, 1, 1, 6,24/
      data (ipa41( 438,j),j=1,7)/  0, 0, 3, 1, 2, 5,24/
      data (ipa41( 439,j),j=1,7)/  0, 0, 3, 1, 3, 4,24/
      data (ipa41( 440,j),j=1,7)/  0, 0, 3, 2, 1, 5,24/
      data (ipa41( 441,j),j=1,7)/  0, 0, 3, 2, 2, 4,24/
      data (ipa41( 442,j),j=1,7)/  0, 0, 3, 2, 3, 3,24/
      data (ipa41( 443,j),j=1,7)/  0, 0, 3, 3, 0, 5,24/
      data (ipa41( 444,j),j=1,7)/  0, 0, 3, 3, 1, 4,24/
      data (ipa41( 445,j),j=1,7)/  0, 0, 3, 3, 2, 3,24/
      data (ipa41( 446,j),j=1,7)/  0, 0, 3, 4, 0, 4,24/
      data (ipa41( 447,j),j=1,7)/  0, 0, 3, 4, 1, 3,24/
      data (ipa41( 448,j),j=1,7)/  0, 0, 3, 4, 2, 2,13/
      data (ipa41( 449,j),j=1,7)/  0, 0, 3, 5, 0, 3,24/
      data (ipa41( 450,j),j=1,7)/  0, 0, 3, 5, 1, 2,24/
      data (ipa41( 451,j),j=1,7)/  0, 0, 3, 6, 0, 2,24/
      data (ipa41( 452,j),j=1,7)/  0, 0, 3, 6, 1, 1,13/
      data (ipa41( 453,j),j=1,7)/  0, 0, 3, 7, 0, 1,24/
      data (ipa41( 454,j),j=1,7)/  0, 0, 4, 1, 1, 5,24/
      data (ipa41( 455,j),j=1,7)/  0, 0, 4, 1, 2, 4,24/
      data (ipa41( 456,j),j=1,7)/  0, 0, 4, 1, 3, 3,13/
      data (ipa41( 457,j),j=1,7)/  0, 0, 4, 2, 1, 4,24/
      data (ipa41( 458,j),j=1,7)/  0, 0, 4, 2, 2, 3,24/
      data (ipa41( 459,j),j=1,7)/  0, 0, 4, 3, 1, 3,24/
      data (ipa41( 460,j),j=1,7)/  0, 0, 4, 3, 2, 2,13/
      data (ipa41( 461,j),j=1,7)/  0, 0, 4, 4, 0, 3,24/
      data (ipa41( 462,j),j=1,7)/  0, 0, 4, 4, 1, 2,24/
      data (ipa41( 463,j),j=1,7)/  0, 0, 4, 5, 0, 2,24/
      data (ipa41( 464,j),j=1,7)/  0, 0, 4, 5, 1, 1,13/
      data (ipa41( 465,j),j=1,7)/  0, 0, 4, 6, 0, 1,24/
      data (ipa41( 466,j),j=1,7)/  0, 0, 5, 1, 1, 4,24/
      data (ipa41( 467,j),j=1,7)/  0, 0, 5, 1, 2, 3,24/
      data (ipa41( 468,j),j=1,7)/  0, 0, 5, 2, 1, 3,24/
      data (ipa41( 469,j),j=1,7)/  0, 0, 5, 2, 2, 2,24/
      data (ipa41( 470,j),j=1,7)/  0, 0, 5, 3, 1, 2,24/
      data (ipa41( 471,j),j=1,7)/  0, 0, 5, 4, 1, 1,13/
      data (ipa41( 472,j),j=1,7)/  0, 0, 5, 5, 0, 1,24/
      data (ipa41( 473,j),j=1,7)/  0, 0, 6, 1, 1, 3,24/
      data (ipa41( 474,j),j=1,7)/  0, 0, 6, 1, 2, 2,13/
      data (ipa41( 475,j),j=1,7)/  0, 0, 6, 2, 1, 2,24/
      data (ipa41( 476,j),j=1,7)/  0, 0, 6, 3, 1, 1,13/
      data (ipa41( 477,j),j=1,7)/  0, 0, 7, 1, 1, 2,24/
      data (ipa41( 478,j),j=1,7)/  0, 0, 7, 2, 1, 1,13/
      data (ipa41( 479,j),j=1,7)/  0, 0, 8, 1, 1, 1,24/
      data (ipa41( 480,j),j=1,7)/  0, 1, 1, 1, 1, 7,24/
      data (ipa41( 481,j),j=1,7)/  0, 1, 1, 1, 2, 6,24/
      data (ipa41( 482,j),j=1,7)/  0, 1, 1, 1, 3, 5,24/
      data (ipa41( 483,j),j=1,7)/  0, 1, 1, 1, 4, 4,24/
      data (ipa41( 484,j),j=1,7)/  0, 1, 1, 1, 5, 3,24/
      data (ipa41( 485,j),j=1,7)/  0, 1, 1, 1, 6, 2,24/
      data (ipa41( 486,j),j=1,7)/  0, 1, 1, 1, 7, 1,24/
      data (ipa41( 487,j),j=1,7)/  0, 1, 1, 1, 8, 0,24/
      data (ipa41( 488,j),j=1,7)/  1, 1, 0, 5, 2, 2,13/
      data (ipa41( 489,j),j=1,7)/  0, 1, 1, 2, 3, 4,24/
      data (ipa41( 490,j),j=1,7)/  0, 1, 1, 2, 4, 3,24/
      data (ipa41( 491,j),j=1,7)/  0, 1, 1, 2, 5, 2,24/
      data (ipa41( 492,j),j=1,7)/  0, 1, 1, 2, 6, 1,24/
      data (ipa41( 493,j),j=1,7)/  0, 1, 1, 2, 7, 0,24/
      data (ipa41( 494,j),j=1,7)/  0, 1, 1, 3, 3, 3,24/
      data (ipa41( 495,j),j=1,7)/  0, 1, 1, 3, 4, 2,24/
      data (ipa41( 496,j),j=1,7)/  0, 1, 1, 3, 5, 1,24/
      data (ipa41( 497,j),j=1,7)/  0, 1, 1, 3, 6, 0,24/
      data (ipa41( 498,j),j=1,7)/  0, 1, 1, 4, 4, 1,24/
      data (ipa41( 499,j),j=1,7)/  0, 1, 1, 4, 5, 0,24/
      data (ipa41( 500,j),j=1,7)/  0, 1, 2, 1, 2, 5,24/
      data (ipa41( 501,j),j=1,7)/  0, 1, 2, 1, 3, 4,24/
      data (ipa41( 502,j),j=1,7)/  0, 1, 2, 1, 4, 3,24/
      data (ipa41( 503,j),j=1,7)/  0, 1, 2, 1, 5, 2,24/
      data (ipa41( 504,j),j=1,7)/  0, 1, 2, 1, 6, 1,24/
      data (ipa41( 505,j),j=1,7)/  1, 2, 0, 5, 2, 1,11/
      data (ipa41( 506,j),j=1,7)/  0, 1, 2, 2, 2, 4,24/
      data (ipa41( 507,j),j=1,7)/  0, 1, 2, 2, 3, 3,24/
      data (ipa41( 508,j),j=1,7)/  0, 1, 2, 2, 4, 2,24/
      data (ipa41( 509,j),j=1,7)/  0, 1, 2, 2, 5, 1,24/
      data (ipa41( 510,j),j=1,7)/  0, 2, 1, 6, 2, 0,11/
      data (ipa41( 511,j),j=1,7)/  0, 1, 2, 3, 1, 4,24/
      data (ipa41( 512,j),j=1,7)/  0, 1, 2, 3, 2, 3,24/
      data (ipa41( 513,j),j=1,7)/  0, 1, 2, 3, 3, 2,24/
      data (ipa41( 514,j),j=1,7)/  0, 1, 2, 3, 4, 1,24/
      data (ipa41( 515,j),j=1,7)/  0, 1, 2, 3, 5, 0,24/
      data (ipa41( 516,j),j=1,7)/  0, 1, 2, 4, 1, 3,24/
      data (ipa41( 517,j),j=1,7)/  0, 1, 2, 4, 2, 2,24/
      data (ipa41( 518,j),j=1,7)/  0, 1, 2, 4, 3, 1,24/
      data (ipa41( 519,j),j=1,7)/  0, 1, 2, 4, 4, 0,24/
      data (ipa41( 520,j),j=1,7)/  0, 1, 2, 5, 1, 2,24/
      data (ipa41( 521,j),j=1,7)/  0, 1, 2, 5, 2, 1,24/
      data (ipa41( 522,j),j=1,7)/  0, 1, 2, 5, 3, 0,24/
      data (ipa41( 523,j),j=1,7)/  0, 1, 2, 6, 1, 1,24/
      data (ipa41( 524,j),j=1,7)/  0, 1, 2, 6, 2, 0,24/
      data (ipa41( 525,j),j=1,7)/  0, 1, 2, 7, 1, 0,11/
      data (ipa41( 526,j),j=1,7)/  0, 1, 3, 1, 3, 3,24/
      data (ipa41( 527,j),j=1,7)/  0, 1, 3, 1, 4, 2,24/
      data (ipa41( 528,j),j=1,7)/  0, 1, 3, 1, 5, 1,24/
      data (ipa41( 529,j),j=1,7)/  0, 1, 3, 2, 2, 3,24/
      data (ipa41( 530,j),j=1,7)/  0, 1, 3, 2, 3, 2,24/
      data (ipa41( 531,j),j=1,7)/  0, 1, 3, 2, 4, 1,24/
      data (ipa41( 532,j),j=1,7)/  0, 1, 3, 3, 1, 3,24/
      data (ipa41( 533,j),j=1,7)/  0, 1, 3, 3, 2, 2,24/
      data (ipa41( 534,j),j=1,7)/  0, 1, 3, 3, 3, 1,24/
      data (ipa41( 535,j),j=1,7)/  0, 3, 1, 4, 3, 0,11/
      data (ipa41( 536,j),j=1,7)/  0, 1, 3, 4, 1, 2,24/
      data (ipa41( 537,j),j=1,7)/  0, 1, 3, 4, 2, 1,24/
      data (ipa41( 538,j),j=1,7)/  0, 1, 3, 4, 3, 0,24/
      data (ipa41( 539,j),j=1,7)/  0, 1, 3, 5, 1, 1,24/
      data (ipa41( 540,j),j=1,7)/  0, 1, 3, 5, 2, 0,24/
      data (ipa41( 541,j),j=1,7)/  0, 1, 3, 6, 1, 0,11/
      data (ipa41( 542,j),j=1,7)/  0, 1, 4, 1, 4, 1,24/
      data (ipa41( 543,j),j=1,7)/  0, 1, 4, 2, 2, 2,24/
      data (ipa41( 544,j),j=1,7)/  0, 1, 4, 2, 3, 1,24/
      data (ipa41( 545,j),j=1,7)/  0, 1, 4, 3, 2, 1,24/
      data (ipa41( 546,j),j=1,7)/  0, 1, 4, 4, 1, 1,24/
      data (ipa41( 547,j),j=1,7)/  0, 4, 1, 2, 4, 0,11/
      data (ipa41( 548,j),j=1,7)/  0, 1, 4, 5, 1, 0,11/
      data (ipa41( 549,j),j=1,7)/  0, 1, 5, 2, 2, 1,24/
      data (ipa41( 550,j),j=1,7)/  0, 2, 2, 2, 2, 3,24/
      data (ipa41( 551,j),j=1,7)/  0, 2, 2, 2, 3, 2,24/
      data (ipa41( 552,j),j=1,7)/  0, 2, 2, 2, 4, 1,24/
      data (ipa41( 553,j),j=1,7)/  0, 2, 2, 2, 5, 0,24/
      data (ipa41( 554,j),j=1,7)/  2, 2, 0, 1, 3, 3,13/
      data (ipa41( 555,j),j=1,7)/  0, 2, 2, 3, 4, 0,24/
      data (ipa41( 556,j),j=1,7)/  0, 2, 3, 2, 3, 1,24/
      data (ipa41( 557,j),j=1,7)/  2, 3, 0, 1, 3, 2,11/
      data (ipa41( 558,j),j=1,7)/  0, 2, 3, 3, 3, 0,24/
      data (ipa41( 559,j),j=1,7)/  0, 2, 3, 4, 2, 0,11/
      data (ipa41( 560,j),j=1,7)/  1, 1, 6, 1, 1, 1, 6/
      data (ipa41( 561,j),j=1,7)/  1, 1, 1, 1, 2, 5,24/
      data (ipa41( 562,j),j=1,7)/  1, 1, 1, 1, 3, 4,24/
      data (ipa41( 563,j),j=1,7)/  1, 1, 1, 2, 2, 4,24/
      data (ipa41( 564,j),j=1,7)/  1, 1, 1, 2, 3, 3,24/
      data (ipa41( 565,j),j=1,7)/  1, 1, 2, 1, 2, 4,24/
      data (ipa41( 566,j),j=1,7)/  1, 1, 2, 1, 3, 3,24/
      data (ipa41( 567,j),j=1,7)/  1, 1, 2, 2, 1, 4,24/
      data (ipa41( 568,j),j=1,7)/  1, 1, 2, 2, 2, 3,24/
      data (ipa41( 569,j),j=1,7)/  1, 1, 2, 3, 1, 3,24/
      data (ipa41( 570,j),j=1,7)/  1, 1, 2, 3, 2, 2,24/
      data (ipa41( 571,j),j=1,7)/  1, 1, 2, 4, 1, 2,24/
      data (ipa41( 572,j),j=1,7)/  1, 1, 2, 5, 1, 1,24/
      data (ipa41( 573,j),j=1,7)/  1, 1, 3, 2, 2, 2,24/
      data (ipa41( 574,j),j=1,7)/  1, 1, 3, 3, 1, 2,24/
      data (ipa41( 575,j),j=1,7)/  1, 1, 3, 4, 1, 1,24/
      data (ipa41( 576,j),j=1,7)/  2, 2, 1, 2, 2, 2, 6/
      data (ipa41( 577,j),j=1,7)/  1, 2, 2, 2, 3, 1,24/
*********** 11 203*******************************
      data (ipa41( 578,j),j=1,7)/  0, 0, 1, 0, 1,10,24/
      data (ipa41( 579,j),j=1,7)/  0, 0, 1, 0, 2, 9,24/
      data (ipa41( 580,j),j=1,7)/  0, 0, 1, 0, 3, 8,24/
      data (ipa41( 581,j),j=1,7)/  0, 0, 1, 0, 4, 7,24/
      data (ipa41( 582,j),j=1,7)/  0, 0, 1, 0, 5, 6,24/
      data (ipa41( 583,j),j=1,7)/  0, 0, 1, 1, 0,10,24/
      data (ipa41( 584,j),j=1,7)/  0, 0, 1, 1, 1, 9,24/
      data (ipa41( 585,j),j=1,7)/  0, 0, 1, 1, 2, 8,24/
      data (ipa41( 586,j),j=1,7)/  0, 0, 1, 1, 3, 7,24/
      data (ipa41( 587,j),j=1,7)/  0, 0, 1, 1, 4, 6,24/
      data (ipa41( 588,j),j=1,7)/  0, 0, 1, 1, 5, 5,24/
      data (ipa41( 589,j),j=1,7)/  0, 0, 1, 2, 0, 9,24/
      data (ipa41( 590,j),j=1,7)/  0, 0, 1, 2, 1, 8,24/
      data (ipa41( 591,j),j=1,7)/  0, 0, 1, 2, 2, 7,24/
      data (ipa41( 592,j),j=1,7)/  0, 0, 1, 2, 3, 6,24/
      data (ipa41( 593,j),j=1,7)/  0, 0, 1, 2, 4, 5,24/
      data (ipa41( 594,j),j=1,7)/  0, 0, 1, 3, 0, 8,24/
      data (ipa41( 595,j),j=1,7)/  0, 0, 1, 3, 1, 7,24/
      data (ipa41( 596,j),j=1,7)/  0, 0, 1, 3, 2, 6,24/
      data (ipa41( 597,j),j=1,7)/  0, 0, 1, 3, 3, 5,24/
      data (ipa41( 598,j),j=1,7)/  0, 0, 1, 3, 4, 4,13/
      data (ipa41( 599,j),j=1,7)/  0, 0, 1, 4, 0, 7,24/
      data (ipa41( 600,j),j=1,7)/  0, 0, 1, 4, 1, 6,24/
      data (ipa41( 601,j),j=1,7)/  0, 0, 1, 4, 2, 5,24/
      data (ipa41( 602,j),j=1,7)/  0, 0, 1, 4, 3, 4,24/
      data (ipa41( 603,j),j=1,7)/  0, 0, 1, 5, 0, 6,24/
      data (ipa41( 604,j),j=1,7)/  0, 0, 1, 5, 1, 5,24/
      data (ipa41( 605,j),j=1,7)/  0, 0, 1, 5, 2, 4,24/
      data (ipa41( 606,j),j=1,7)/  0, 0, 1, 5, 3, 3,13/
      data (ipa41( 607,j),j=1,7)/  0, 0, 1, 6, 0, 5,24/
      data (ipa41( 608,j),j=1,7)/  0, 0, 1, 6, 1, 4,24/
      data (ipa41( 609,j),j=1,7)/  0, 0, 1, 6, 2, 3,24/
      data (ipa41( 610,j),j=1,7)/  0, 0, 1, 7, 0, 4,24/
      data (ipa41( 611,j),j=1,7)/  0, 0, 1, 7, 1, 3,24/
      data (ipa41( 612,j),j=1,7)/  0, 0, 1, 7, 2, 2,13/
      data (ipa41( 613,j),j=1,7)/  0, 0, 1, 8, 0, 3,24/
      data (ipa41( 614,j),j=1,7)/  0, 0, 1, 8, 1, 2,24/
      data (ipa41( 615,j),j=1,7)/  0, 0, 1, 9, 0, 2,24/
      data (ipa41( 616,j),j=1,7)/  0, 0, 1, 9, 1, 1,24/
      data (ipa41( 617,j),j=1,7)/  0, 0, 1,10, 0, 1,24/
      data (ipa41( 618,j),j=1,7)/  0, 0, 2, 0, 2, 8,24/
      data (ipa41( 619,j),j=1,7)/  0, 0, 2, 0, 3, 7,24/
      data (ipa41( 620,j),j=1,7)/  0, 0, 2, 0, 4, 6,24/
      data (ipa41( 621,j),j=1,7)/  0, 0, 2, 0, 5, 5,24/
      data (ipa41( 622,j),j=1,7)/  0, 0, 2, 1, 1, 8,24/
      data (ipa41( 623,j),j=1,7)/  0, 0, 2, 1, 2, 7,24/
      data (ipa41( 624,j),j=1,7)/  0, 0, 2, 1, 3, 6,24/
      data (ipa41( 625,j),j=1,7)/  0, 0, 2, 1, 4, 5,24/
      data (ipa41( 626,j),j=1,7)/  0, 0, 2, 2, 0, 8,24/
      data (ipa41( 627,j),j=1,7)/  0, 0, 2, 2, 1, 7,24/
      data (ipa41( 628,j),j=1,7)/  0, 0, 2, 2, 2, 6,24/
      data (ipa41( 629,j),j=1,7)/  0, 0, 2, 2, 3, 5,24/
      data (ipa41( 630,j),j=1,7)/  0, 0, 2, 2, 4, 4,24/
      data (ipa41( 631,j),j=1,7)/  0, 0, 2, 3, 0, 7,24/
      data (ipa41( 632,j),j=1,7)/  0, 0, 2, 3, 1, 6,24/
      data (ipa41( 633,j),j=1,7)/  0, 0, 2, 3, 2, 5,24/
      data (ipa41( 634,j),j=1,7)/  0, 0, 2, 3, 3, 4,24/
      data (ipa41( 635,j),j=1,7)/  0, 0, 2, 4, 0, 6,24/
      data (ipa41( 636,j),j=1,7)/  0, 0, 2, 4, 1, 5,24/
      data (ipa41( 637,j),j=1,7)/  0, 0, 2, 4, 2, 4,24/
      data (ipa41( 638,j),j=1,7)/  0, 0, 2, 4, 3, 3,13/
      data (ipa41( 639,j),j=1,7)/  0, 0, 2, 5, 0, 5,24/
      data (ipa41( 640,j),j=1,7)/  0, 0, 2, 5, 1, 4,24/
      data (ipa41( 641,j),j=1,7)/  0, 0, 2, 5, 2, 3,24/
      data (ipa41( 642,j),j=1,7)/  0, 0, 2, 6, 0, 4,24/
      data (ipa41( 643,j),j=1,7)/  0, 0, 2, 6, 1, 3,24/
      data (ipa41( 644,j),j=1,7)/  0, 0, 2, 6, 2, 2,24/
      data (ipa41( 645,j),j=1,7)/  0, 0, 2, 7, 0, 3,24/
      data (ipa41( 646,j),j=1,7)/  0, 0, 2, 7, 1, 2,24/
      data (ipa41( 647,j),j=1,7)/  0, 0, 2, 8, 0, 2,24/
      data (ipa41( 648,j),j=1,7)/  0, 0, 2, 8, 1, 1,13/
      data (ipa41( 649,j),j=1,7)/  0, 0, 2, 9, 0, 1,24/
      data (ipa41( 650,j),j=1,7)/  0, 0, 3, 0, 3, 6,24/
      data (ipa41( 651,j),j=1,7)/  0, 0, 3, 0, 4, 5,24/
      data (ipa41( 652,j),j=1,7)/  0, 0, 3, 1, 1, 7,24/
      data (ipa41( 653,j),j=1,7)/  0, 0, 3, 1, 2, 6,24/
      data (ipa41( 654,j),j=1,7)/  0, 0, 3, 1, 3, 5,24/
      data (ipa41( 655,j),j=1,7)/  0, 0, 3, 1, 4, 4,13/
      data (ipa41( 656,j),j=1,7)/  0, 0, 3, 2, 1, 6,24/
      data (ipa41( 657,j),j=1,7)/  0, 0, 3, 2, 2, 5,24/
      data (ipa41( 658,j),j=1,7)/  0, 0, 3, 2, 3, 4,24/
      data (ipa41( 659,j),j=1,7)/  0, 0, 3, 3, 0, 6,24/
      data (ipa41( 660,j),j=1,7)/  0, 0, 3, 3, 1, 5,24/
      data (ipa41( 661,j),j=1,7)/  0, 0, 3, 3, 2, 4,24/
      data (ipa41( 662,j),j=1,7)/  3, 3, 3, 3, 0, 0,13/
      data (ipa41( 663,j),j=1,7)/  0, 0, 3, 4, 0, 5,24/
      data (ipa41( 664,j),j=1,7)/  0, 0, 3, 4, 1, 4,24/
      data (ipa41( 665,j),j=1,7)/  0, 0, 3, 4, 2, 3,24/
      data (ipa41( 666,j),j=1,7)/  0, 0, 3, 5, 0, 4,24/
      data (ipa41( 667,j),j=1,7)/  0, 0, 3, 5, 1, 3,24/
      data (ipa41( 668,j),j=1,7)/  0, 0, 3, 5, 2, 2,13/
      data (ipa41( 669,j),j=1,7)/  0, 0, 3, 6, 0, 3,24/
      data (ipa41( 670,j),j=1,7)/  0, 0, 3, 6, 1, 2,24/
      data (ipa41( 671,j),j=1,7)/  0, 0, 3, 7, 0, 2,24/
      data (ipa41( 672,j),j=1,7)/  0, 0, 3, 7, 1, 1,13/
      data (ipa41( 673,j),j=1,7)/  0, 0, 3, 8, 0, 1,24/
      data (ipa41( 674,j),j=1,7)/  0, 0, 4, 0, 4, 4,24/
      data (ipa41( 675,j),j=1,7)/  0, 0, 4, 1, 1, 6,24/
      data (ipa41( 676,j),j=1,7)/  0, 0, 4, 1, 2, 5,24/
      data (ipa41( 677,j),j=1,7)/  0, 0, 4, 1, 3, 4,24/
      data (ipa41( 678,j),j=1,7)/  0, 0, 4, 2, 1, 5,24/
      data (ipa41( 679,j),j=1,7)/  0, 0, 4, 2, 2, 4,24/
      data (ipa41( 680,j),j=1,7)/  0, 0, 4, 2, 3, 3,13/
      data (ipa41( 681,j),j=1,7)/  0, 0, 4, 3, 1, 4,24/
      data (ipa41( 682,j),j=1,7)/  0, 0, 4, 3, 2, 3,24/
      data (ipa41( 683,j),j=1,7)/  0, 4, 0, 4, 4, 0,11/
      data (ipa41( 684,j),j=1,7)/  0, 0, 4, 4, 1, 3,24/
      data (ipa41( 685,j),j=1,7)/  0, 0, 4, 4, 2, 2,24/
      data (ipa41( 686,j),j=1,7)/  0, 0, 4, 5, 0, 3,24/
      data (ipa41( 687,j),j=1,7)/  0, 0, 4, 5, 1, 2,24/
      data (ipa41( 688,j),j=1,7)/  0, 0, 4, 6, 0, 2,24/
      data (ipa41( 689,j),j=1,7)/  0, 0, 4, 6, 1, 1,13/
      data (ipa41( 690,j),j=1,7)/  0, 0, 4, 7, 0, 1,24/
      data (ipa41( 691,j),j=1,7)/  0, 0, 5, 1, 1, 5,24/
      data (ipa41( 692,j),j=1,7)/  0, 0, 5, 1, 2, 4,24/
      data (ipa41( 693,j),j=1,7)/  0, 0, 5, 1, 3, 3,13/
      data (ipa41( 694,j),j=1,7)/  0, 0, 5, 2, 1, 4,24/
      data (ipa41( 695,j),j=1,7)/  0, 0, 5, 2, 2, 3,24/
      data (ipa41( 696,j),j=1,7)/  0, 0, 5, 3, 1, 3,24/
      data (ipa41( 697,j),j=1,7)/  0, 0, 5, 3, 2, 2,13/
      data (ipa41( 698,j),j=1,7)/  0, 0, 5, 4, 1, 2,24/
      data (ipa41( 699,j),j=1,7)/  0, 0, 5, 5, 0, 2,24/
      data (ipa41( 700,j),j=1,7)/  0, 0, 5, 5, 1, 1,24/
      data (ipa41( 701,j),j=1,7)/  0, 0, 5, 6, 0, 1,24/
      data (ipa41( 702,j),j=1,7)/  0, 0, 6, 1, 1, 4,24/
      data (ipa41( 703,j),j=1,7)/  0, 0, 6, 1, 2, 3,24/
      data (ipa41( 704,j),j=1,7)/  0, 0, 6, 2, 1, 3,24/
      data (ipa41( 705,j),j=1,7)/  0, 0, 6, 2, 2, 2,24/
      data (ipa41( 706,j),j=1,7)/  0, 0, 6, 3, 1, 2,24/
      data (ipa41( 707,j),j=1,7)/  0, 0, 6, 4, 1, 1,13/
      data (ipa41( 708,j),j=1,7)/  0, 0, 7, 1, 1, 3,24/
      data (ipa41( 709,j),j=1,7)/  0, 0, 7, 1, 2, 2,13/
      data (ipa41( 710,j),j=1,7)/  0, 0, 7, 2, 1, 2,24/
      data (ipa41( 711,j),j=1,7)/  0, 0, 7, 3, 1, 1,13/
      data (ipa41( 712,j),j=1,7)/  0, 0, 8, 1, 1, 2,24/
      data (ipa41( 713,j),j=1,7)/  0, 0, 8, 2, 1, 1,13/
      data (ipa41( 714,j),j=1,7)/  0, 0, 9, 1, 1, 1,24/
      data (ipa41( 715,j),j=1,7)/  0, 1, 1, 1, 1, 8,24/
      data (ipa41( 716,j),j=1,7)/  0, 1, 1, 1, 2, 7,24/
      data (ipa41( 717,j),j=1,7)/  0, 1, 1, 1, 3, 6,24/
      data (ipa41( 718,j),j=1,7)/  0, 1, 1, 1, 4, 5,24/
      data (ipa41( 719,j),j=1,7)/  0, 1, 1, 1, 5, 4,24/
      data (ipa41( 720,j),j=1,7)/  0, 1, 1, 1, 6, 3,24/
      data (ipa41( 721,j),j=1,7)/  0, 1, 1, 1, 7, 2,24/
      data (ipa41( 722,j),j=1,7)/  0, 1, 1, 1, 8, 1,24/
      data (ipa41( 723,j),j=1,7)/  0, 1, 1, 1, 9, 0,24/
      data (ipa41( 724,j),j=1,7)/  1, 1, 0, 6, 2, 2,13/
      data (ipa41( 725,j),j=1,7)/  0, 1, 1, 2, 3, 5,24/
      data (ipa41( 726,j),j=1,7)/  0, 1, 1, 2, 4, 4,24/
      data (ipa41( 727,j),j=1,7)/  0, 1, 1, 2, 5, 3,24/
      data (ipa41( 728,j),j=1,7)/  0, 1, 1, 2, 6, 2,24/
      data (ipa41( 729,j),j=1,7)/  0, 1, 1, 2, 7, 1,24/
      data (ipa41( 730,j),j=1,7)/  0, 1, 1, 2, 8, 0,24/
      data (ipa41( 731,j),j=1,7)/  1, 1, 0, 4, 3, 3,13/
      data (ipa41( 732,j),j=1,7)/  0, 1, 1, 3, 4, 3,24/
      data (ipa41( 733,j),j=1,7)/  0, 1, 1, 3, 5, 2,24/
      data (ipa41( 734,j),j=1,7)/  0, 1, 1, 3, 6, 1,24/
      data (ipa41( 735,j),j=1,7)/  0, 1, 1, 3, 7, 0,24/
      data (ipa41( 736,j),j=1,7)/  1, 1, 0, 2, 4, 4,13/
      data (ipa41( 737,j),j=1,7)/  0, 1, 1, 4, 5, 1,24/
      data (ipa41( 738,j),j=1,7)/  0, 1, 1, 4, 6, 0,24/
      data (ipa41( 739,j),j=1,7)/  0, 1, 1, 5, 5, 0,24/
      data (ipa41( 740,j),j=1,7)/  0, 1, 2, 1, 2, 6,24/
      data (ipa41( 741,j),j=1,7)/  0, 1, 2, 1, 3, 5,24/
      data (ipa41( 742,j),j=1,7)/  0, 1, 2, 1, 4, 4,24/
      data (ipa41( 743,j),j=1,7)/  0, 1, 2, 1, 5, 3,24/
      data (ipa41( 744,j),j=1,7)/  0, 1, 2, 1, 6, 2,24/
      data (ipa41( 745,j),j=1,7)/  0, 1, 2, 1, 7, 1,24/
      data (ipa41( 746,j),j=1,7)/  1, 2, 0, 6, 2, 1,11/
      data (ipa41( 747,j),j=1,7)/  0, 1, 2, 2, 2, 5,24/
      data (ipa41( 748,j),j=1,7)/  0, 1, 2, 2, 3, 4,24/
      data (ipa41( 749,j),j=1,7)/  0, 1, 2, 2, 4, 3,24/
      data (ipa41( 750,j),j=1,7)/  0, 1, 2, 2, 5, 2,24/
      data (ipa41( 751,j),j=1,7)/  0, 1, 2, 2, 6, 1,24/
      data (ipa41( 752,j),j=1,7)/  0, 2, 1, 7, 2, 0,11/
      data (ipa41( 753,j),j=1,7)/  0, 1, 2, 3, 1, 5,24/
      data (ipa41( 754,j),j=1,7)/  0, 1, 2, 3, 2, 4,24/
      data (ipa41( 755,j),j=1,7)/  0, 1, 2, 3, 3, 3,24/
      data (ipa41( 756,j),j=1,7)/  0, 1, 2, 3, 4, 2,24/
      data (ipa41( 757,j),j=1,7)/  0, 1, 2, 3, 5, 1,24/
      data (ipa41( 758,j),j=1,7)/  0, 1, 2, 3, 6, 0,24/
      data (ipa41( 759,j),j=1,7)/  0, 1, 2, 4, 1, 4,24/
      data (ipa41( 760,j),j=1,7)/  0, 1, 2, 4, 2, 3,24/
      data (ipa41( 761,j),j=1,7)/  0, 1, 2, 4, 3, 2,24/
      data (ipa41( 762,j),j=1,7)/  0, 1, 2, 4, 4, 1,24/
      data (ipa41( 763,j),j=1,7)/  0, 1, 2, 4, 5, 0,24/
      data (ipa41( 764,j),j=1,7)/  0, 1, 2, 5, 1, 3,24/
      data (ipa41( 765,j),j=1,7)/  0, 1, 2, 5, 2, 2,24/
      data (ipa41( 766,j),j=1,7)/  0, 1, 2, 5, 3, 1,24/
      data (ipa41( 767,j),j=1,7)/  0, 1, 2, 5, 4, 0,24/
      data (ipa41( 768,j),j=1,7)/  0, 1, 2, 6, 1, 2,24/
      data (ipa41( 769,j),j=1,7)/  0, 1, 2, 6, 2, 1,24/
      data (ipa41( 770,j),j=1,7)/  0, 1, 2, 6, 3, 0,24/
      data (ipa41( 771,j),j=1,7)/  0, 1, 2, 7, 1, 1,24/
      data (ipa41( 772,j),j=1,7)/  0, 1, 2, 7, 2, 0,24/
      data (ipa41( 773,j),j=1,7)/  0, 1, 2, 8, 1, 0,11/
      data (ipa41( 774,j),j=1,7)/  0, 1, 3, 1, 3, 4,24/
      data (ipa41( 775,j),j=1,7)/  0, 1, 3, 1, 4, 3,24/
      data (ipa41( 776,j),j=1,7)/  0, 1, 3, 1, 5, 2,24/
      data (ipa41( 777,j),j=1,7)/  0, 1, 3, 1, 6, 1,24/
      data (ipa41( 778,j),j=1,7)/  0, 1, 3, 2, 2, 4,24/
      data (ipa41( 779,j),j=1,7)/  0, 1, 3, 2, 3, 3,24/
      data (ipa41( 780,j),j=1,7)/  0, 1, 3, 2, 4, 2,24/
      data (ipa41( 781,j),j=1,7)/  0, 1, 3, 2, 5, 1,24/
      data (ipa41( 782,j),j=1,7)/  1, 3, 0, 4, 3, 1,11/
      data (ipa41( 783,j),j=1,7)/  0, 1, 3, 3, 2, 3,24/
      data (ipa41( 784,j),j=1,7)/  0, 1, 3, 3, 3, 2,24/
      data (ipa41( 785,j),j=1,7)/  0, 1, 3, 3, 4, 1,24/
      data (ipa41( 786,j),j=1,7)/  0, 3, 1, 5, 3, 0,11/
      data (ipa41( 787,j),j=1,7)/  0, 1, 3, 4, 1, 3,24/
      data (ipa41( 788,j),j=1,7)/  0, 1, 3, 4, 2, 2,24/
      data (ipa41( 789,j),j=1,7)/  0, 1, 3, 4, 3, 1,24/
      data (ipa41( 790,j),j=1,7)/  0, 1, 3, 4, 4, 0,24/
      data (ipa41( 791,j),j=1,7)/  0, 1, 3, 5, 1, 2,24/
      data (ipa41( 792,j),j=1,7)/  0, 1, 3, 5, 2, 1,24/
      data (ipa41( 793,j),j=1,7)/  0, 1, 3, 5, 3, 0,24/
      data (ipa41( 794,j),j=1,7)/  0, 1, 3, 6, 1, 1,24/
      data (ipa41( 795,j),j=1,7)/  0, 1, 3, 6, 2, 0,24/
      data (ipa41( 796,j),j=1,7)/  0, 1, 3, 7, 1, 0,11/
      data (ipa41( 797,j),j=1,7)/  0, 1, 4, 1, 4, 2,24/
      data (ipa41( 798,j),j=1,7)/  0, 1, 4, 1, 5, 1,24/
      data (ipa41( 799,j),j=1,7)/  0, 1, 4, 2, 2, 3,24/
      data (ipa41( 800,j),j=1,7)/  0, 1, 4, 2, 3, 2,24/
      data (ipa41( 801,j),j=1,7)/  0, 1, 4, 2, 4, 1,24/
      data (ipa41( 802,j),j=1,7)/  0, 1, 4, 3, 2, 2,24/
      data (ipa41( 803,j),j=1,7)/  0, 1, 4, 3, 3, 1,24/
      data (ipa41( 804,j),j=1,7)/  1, 4, 0, 2, 4, 1,11/
      data (ipa41( 805,j),j=1,7)/  0, 1, 4, 4, 2, 1,24/
      data (ipa41( 806,j),j=1,7)/  0, 4, 1, 3, 4, 0,11/
      data (ipa41( 807,j),j=1,7)/  0, 1, 4, 5, 1, 1,24/
      data (ipa41( 808,j),j=1,7)/  0, 1, 4, 5, 2, 0,24/
      data (ipa41( 809,j),j=1,7)/  0, 1, 4, 6, 1, 0,11/
      data (ipa41( 810,j),j=1,7)/  0, 1, 5, 2, 2, 2,24/
      data (ipa41( 811,j),j=1,7)/  0, 1, 5, 2, 3, 1,24/
      data (ipa41( 812,j),j=1,7)/  0, 1, 5, 3, 2, 1,24/
      data (ipa41( 813,j),j=1,7)/  0, 1, 5, 5, 1, 0, 6/
      data (ipa41( 814,j),j=1,7)/  0, 1, 6, 2, 2, 1,24/
      data (ipa41( 815,j),j=1,7)/  0, 2, 2, 2, 2, 4,24/
      data (ipa41( 816,j),j=1,7)/  0, 2, 2, 2, 3, 3,24/
      data (ipa41( 817,j),j=1,7)/  0, 2, 2, 2, 4, 2,24/
      data (ipa41( 818,j),j=1,7)/  0, 2, 2, 2, 5, 1,24/
      data (ipa41( 819,j),j=1,7)/  0, 2, 2, 2, 6, 0,24/
      data (ipa41( 820,j),j=1,7)/  0, 2, 2, 3, 3, 2,24/
      data (ipa41( 821,j),j=1,7)/  0, 2, 2, 3, 4, 1,24/
      data (ipa41( 822,j),j=1,7)/  0, 2, 2, 3, 5, 0,24/
      data (ipa41( 823,j),j=1,7)/  0, 2, 2, 4, 4, 0,24/
      data (ipa41( 824,j),j=1,7)/  0, 2, 3, 2, 3, 2,24/
      data (ipa41( 825,j),j=1,7)/  0, 2, 3, 2, 4, 1,24/
      data (ipa41( 826,j),j=1,7)/  0, 2, 3, 3, 2, 2,24/
      data (ipa41( 827,j),j=1,7)/  0, 2, 3, 3, 3, 1,24/
      data (ipa41( 828,j),j=1,7)/  0, 3, 2, 4, 3, 0,11/
      data (ipa41( 829,j),j=1,7)/  0, 2, 3, 4, 2, 1,24/
      data (ipa41( 830,j),j=1,7)/  0, 2, 3, 4, 3, 0,24/
      data (ipa41( 831,j),j=1,7)/  0, 2, 3, 5, 2, 0,11/
      data (ipa41( 832,j),j=1,7)/  0, 2, 4, 4, 2, 0, 6/
      data (ipa41( 833,j),j=1,7)/  3, 3, 0, 0, 3, 3, 3/
      data (ipa41( 834,j),j=1,7)/  1, 1, 7, 1, 1, 1, 6/
      data (ipa41( 835,j),j=1,7)/  1, 1, 1, 1, 2, 6,24/
      data (ipa41( 836,j),j=1,7)/  1, 1, 1, 1, 3, 5,24/
      data (ipa41( 837,j),j=1,7)/  1, 1, 1, 1, 4, 4,13/
      data (ipa41( 838,j),j=1,7)/  1, 1, 1, 2, 2, 5,24/
      data (ipa41( 839,j),j=1,7)/  1, 1, 1, 2, 3, 4,24/
      data (ipa41( 840,j),j=1,7)/  1, 1, 1, 3, 3, 3, 4/
      data (ipa41( 841,j),j=1,7)/  1, 1, 2, 1, 2, 5,24/
      data (ipa41( 842,j),j=1,7)/  1, 1, 2, 1, 3, 4,24/
      data (ipa41( 843,j),j=1,7)/  1, 1, 2, 2, 1, 5,24/
      data (ipa41( 844,j),j=1,7)/  1, 1, 2, 2, 2, 4,24/
      data (ipa41( 845,j),j=1,7)/  1, 1, 2, 2, 3, 3,24/
      data (ipa41( 846,j),j=1,7)/  1, 1, 2, 3, 1, 4,24/
      data (ipa41( 847,j),j=1,7)/  1, 1, 2, 3, 2, 3,24/
      data (ipa41( 848,j),j=1,7)/  1, 1, 2, 4, 1, 3,24/
      data (ipa41( 849,j),j=1,7)/  1, 1, 2, 4, 2, 2,24/
      data (ipa41( 850,j),j=1,7)/  1, 1, 2, 5, 1, 2,24/
      data (ipa41( 851,j),j=1,7)/  1, 1, 2, 6, 1, 1,24/
      data (ipa41( 852,j),j=1,7)/  1, 1, 3, 1, 3, 3,24/
      data (ipa41( 853,j),j=1,7)/  1, 1, 3, 2, 2, 3,24/
      data (ipa41( 854,j),j=1,7)/  1, 3, 1, 3, 3, 1,11/
      data (ipa41( 855,j),j=1,7)/  1, 1, 3, 3, 2, 2,24/
      data (ipa41( 856,j),j=1,7)/  1, 1, 3, 4, 1, 2,24/
      data (ipa41( 857,j),j=1,7)/  1, 1, 3, 5, 1, 1,24/
      data (ipa41( 858,j),j=1,7)/  1, 1, 4, 2, 2, 2,24/
      data (ipa41( 859,j),j=1,7)/  1, 1, 4, 4, 1, 1, 3/
      data (ipa41( 860,j),j=1,7)/  1, 2, 2, 2, 2, 3,24/
      data (ipa41( 861,j),j=1,7)/  1, 2, 2, 2, 3, 2,24/
      data (ipa41( 862,j),j=1,7)/  1, 2, 2, 2, 4, 1,24/
      data (ipa41( 863,j),j=1,7)/  1, 2, 2, 3, 3, 1,24/
      data (ipa41( 864,j),j=1,7)/  1, 2, 3, 3, 2, 1, 6/
      data (ipa41( 865,j),j=1,7)/  2, 2, 2, 2, 2, 2, 1/
*********** 12 288*******************************
      data (nfas41(i),i=1,12)/0,0,2,7,18,40,77,138,233,374,577,865/
*********** 0 0*******************************
*********** 1 0*******************************
*********** 2 0*******************************
      data (ipa42(   1,j),j=1,7)/  0, 0, 1, 0, 1, 1, 2/
      data (ipa42(   2,j),j=1,7)/  0, 0, 1, 1, 0, 1, 3/
      data (ipa42(   3,j),j=1,7)/  0, 1, 1, 0, 1, 0, 3/
      data (ipa42(   4,j),j=1,7)/  1, 1, 1, 0, 0, 0, 1/
*********** 3 4*******************************
      data (ipa42(   5,j),j=1,7)/  0, 0, 1, 0, 1, 2, 3/
      data (ipa42(   6,j),j=1,7)/  0, 0, 1, 1, 0, 2, 3/
      data (ipa42(   7,j),j=1,7)/  0, 0, 1, 1, 1, 1, 2/
      data (ipa42(   8,j),j=1,7)/  0, 0, 1, 2, 0, 1, 3/
      data (ipa42(   9,j),j=1,7)/  0, 0, 2, 0, 1, 1, 2/
      data (ipa42(  10,j),j=1,7)/  0, 0, 2, 1, 0, 1, 3/
      data (ipa42(  11,j),j=1,7)/  0, 1, 1, 0, 1, 1, 3/
      data (ipa42(  12,j),j=1,7)/  0, 1, 1, 0, 2, 0, 3/
      data (ipa42(  13,j),j=1,7)/  1, 1, 0, 0, 1, 1, 2/
      data (ipa42(  14,j),j=1,7)/  0, 1, 2, 0, 1, 0, 3/
      data (ipa42(  15,j),j=1,7)/  0, 1, 2, 1, 0, 0, 3/
      data (ipa42(  16,j),j=1,7)/  1, 1, 1, 1, 0, 0, 2/
      data (ipa42(  17,j),j=1,7)/  1, 1, 2, 0, 0, 0, 2/
*********** 4 13*******************************
      data (ipa42(  18,j),j=1,7)/  0, 0, 1, 0, 1, 3, 3/
      data (ipa42(  19,j),j=1,7)/  0, 0, 1, 0, 2, 2, 2/
      data (ipa42(  20,j),j=1,7)/  0, 0, 1, 1, 0, 3, 3/
      data (ipa42(  21,j),j=1,7)/  0, 0, 1, 1, 1, 2, 3/
      data (ipa42(  22,j),j=1,7)/  0, 0, 1, 2, 0, 2, 3/
      data (ipa42(  23,j),j=1,7)/  0, 0, 1, 2, 1, 1, 2/
      data (ipa42(  24,j),j=1,7)/  0, 0, 1, 3, 0, 1, 3/
      data (ipa42(  25,j),j=1,7)/  0, 0, 2, 0, 1, 2, 3/
      data (ipa42(  26,j),j=1,7)/  0, 0, 2, 1, 0, 2, 3/
      data (ipa42(  27,j),j=1,7)/  0, 0, 2, 1, 1, 1, 2/
      data (ipa42(  28,j),j=1,7)/  0, 0, 2, 2, 0, 1, 3/
      data (ipa42(  29,j),j=1,7)/  0, 0, 3, 0, 1, 1, 2/
      data (ipa42(  30,j),j=1,7)/  0, 0, 3, 1, 0, 1, 3/
      data (ipa42(  31,j),j=1,7)/  0, 1, 1, 0, 1, 2, 3/
      data (ipa42(  32,j),j=1,7)/  0, 1, 1, 0, 2, 1, 3/
      data (ipa42(  33,j),j=1,7)/  0, 1, 1, 0, 3, 0, 3/
      data (ipa42(  34,j),j=1,7)/  1, 1, 0, 1, 1, 1, 2/
      data (ipa42(  35,j),j=1,7)/  0, 1, 1, 1, 2, 0, 3/
      data (ipa42(  36,j),j=1,7)/  0, 1, 2, 0, 1, 1, 3/
      data (ipa42(  37,j),j=1,7)/  0, 1, 2, 0, 2, 0, 3/
      data (ipa42(  38,j),j=1,7)/  0, 1, 2, 1, 0, 1, 3/
      data (ipa42(  39,j),j=1,7)/  0, 1, 2, 1, 1, 0, 3/
      data (ipa42(  40,j),j=1,7)/  0, 1, 2, 2, 0, 0, 3/
      data (ipa42(  41,j),j=1,7)/  0, 1, 3, 0, 1, 0, 3/
      data (ipa42(  42,j),j=1,7)/  0, 1, 3, 1, 0, 0, 3/
      data (ipa42(  43,j),j=1,7)/  0, 2, 2, 0, 1, 0, 3/
      data (ipa42(  44,j),j=1,7)/  1, 1, 1, 2, 0, 0, 2/
      data (ipa42(  45,j),j=1,7)/  1, 1, 1, 0, 1, 1, 2/
      data (ipa42(  46,j),j=1,7)/  1, 1, 2, 0, 0, 1, 3/
      data (ipa42(  47,j),j=1,7)/  1, 1, 2, 1, 0, 0, 2/
      data (ipa42(  48,j),j=1,7)/  1, 1, 3, 0, 0, 0, 2/
      data (ipa42(  49,j),j=1,7)/  2, 2, 1, 0, 0, 0, 2/
*********** 5 32*******************************
      data (ipa42(  50,j),j=1,7)/  0, 0, 1, 0, 1, 4, 3/
      data (ipa42(  51,j),j=1,7)/  0, 0, 1, 0, 2, 3, 3/
      data (ipa42(  52,j),j=1,7)/  0, 0, 1, 1, 0, 4, 3/
      data (ipa42(  53,j),j=1,7)/  0, 0, 1, 1, 1, 3, 3/
      data (ipa42(  54,j),j=1,7)/  0, 0, 1, 1, 2, 2, 2/
      data (ipa42(  55,j),j=1,7)/  0, 0, 1, 2, 0, 3, 3/
      data (ipa42(  56,j),j=1,7)/  0, 0, 1, 2, 1, 2, 3/
      data (ipa42(  57,j),j=1,7)/  0, 0, 1, 3, 0, 2, 3/
      data (ipa42(  58,j),j=1,7)/  0, 0, 1, 3, 1, 1, 2/
      data (ipa42(  59,j),j=1,7)/  0, 0, 1, 4, 0, 1, 3/
      data (ipa42(  60,j),j=1,7)/  0, 0, 2, 0, 1, 3, 3/
      data (ipa42(  61,j),j=1,7)/  0, 0, 2, 0, 2, 2, 2/
      data (ipa42(  62,j),j=1,7)/  0, 0, 2, 1, 0, 3, 3/
      data (ipa42(  63,j),j=1,7)/  0, 0, 2, 1, 1, 2, 3/
      data (ipa42(  64,j),j=1,7)/  0, 0, 2, 2, 0, 2, 3/
      data (ipa42(  65,j),j=1,7)/  0, 0, 2, 2, 1, 1, 2/
      data (ipa42(  66,j),j=1,7)/  0, 0, 2, 3, 0, 1, 3/
      data (ipa42(  67,j),j=1,7)/  0, 0, 3, 0, 1, 2, 3/
      data (ipa42(  68,j),j=1,7)/  0, 0, 3, 1, 0, 2, 3/
      data (ipa42(  69,j),j=1,7)/  0, 0, 3, 1, 1, 1, 2/
      data (ipa42(  70,j),j=1,7)/  0, 0, 3, 2, 0, 1, 3/
      data (ipa42(  71,j),j=1,7)/  0, 0, 4, 0, 1, 1, 2/
      data (ipa42(  72,j),j=1,7)/  0, 0, 4, 1, 0, 1, 3/
      data (ipa42(  73,j),j=1,7)/  0, 1, 1, 0, 1, 3, 3/
      data (ipa42(  74,j),j=1,7)/  0, 1, 1, 0, 2, 2, 3/
      data (ipa42(  75,j),j=1,7)/  0, 1, 1, 0, 3, 1, 3/
      data (ipa42(  76,j),j=1,7)/  0, 1, 1, 0, 4, 0, 3/
      data (ipa42(  77,j),j=1,7)/  1, 1, 0, 2, 1, 1, 2/
      data (ipa42(  78,j),j=1,7)/  0, 1, 1, 1, 2, 1, 3/
      data (ipa42(  79,j),j=1,7)/  0, 1, 1, 1, 3, 0, 3/
      data (ipa42(  80,j),j=1,7)/  1, 1, 0, 0, 2, 2, 2/
      data (ipa42(  81,j),j=1,7)/  0, 1, 2, 0, 1, 2, 3/
      data (ipa42(  82,j),j=1,7)/  0, 1, 2, 0, 2, 1, 3/
      data (ipa42(  83,j),j=1,7)/  0, 1, 2, 0, 3, 0, 3/
      data (ipa42(  84,j),j=1,7)/  0, 1, 2, 1, 0, 2, 3/
      data (ipa42(  85,j),j=1,7)/  0, 1, 2, 1, 1, 1, 3/
      data (ipa42(  86,j),j=1,7)/  0, 1, 2, 1, 2, 0, 3/
      data (ipa42(  87,j),j=1,7)/  0, 1, 2, 2, 0, 1, 3/
      data (ipa42(  88,j),j=1,7)/  0, 1, 2, 2, 1, 0, 3/
      data (ipa42(  89,j),j=1,7)/  0, 1, 2, 3, 0, 0, 3/
      data (ipa42(  90,j),j=1,7)/  0, 1, 3, 0, 1, 1, 3/
      data (ipa42(  91,j),j=1,7)/  0, 1, 3, 0, 2, 0, 3/
      data (ipa42(  92,j),j=1,7)/  0, 1, 3, 1, 0, 1, 3/
      data (ipa42(  93,j),j=1,7)/  0, 1, 3, 1, 1, 0, 3/
      data (ipa42(  94,j),j=1,7)/  0, 1, 3, 2, 0, 0, 3/
      data (ipa42(  95,j),j=1,7)/  0, 1, 4, 0, 1, 0, 3/
      data (ipa42(  96,j),j=1,7)/  0, 1, 4, 1, 0, 0, 3/
      data (ipa42(  97,j),j=1,7)/  0, 2, 2, 0, 1, 1, 3/
      data (ipa42(  98,j),j=1,7)/  0, 2, 2, 0, 2, 0, 3/
      data (ipa42(  99,j),j=1,7)/  2, 2, 0, 0, 1, 1, 2/
      data (ipa42( 100,j),j=1,7)/  0, 2, 3, 0, 1, 0, 3/
      data (ipa42( 101,j),j=1,7)/  0, 2, 3, 1, 0, 0, 3/
      data (ipa42( 102,j),j=1,7)/  1, 1, 1, 3, 0, 0, 2/
      data (ipa42( 103,j),j=1,7)/  1, 1, 1, 0, 1, 2, 3/
      data (ipa42( 104,j),j=1,7)/  1, 1, 1, 1, 1, 1, 1/
      data (ipa42( 105,j),j=1,7)/  1, 1, 2, 0, 0, 2, 3/
      data (ipa42( 106,j),j=1,7)/  1, 1, 2, 0, 1, 1, 2/
      data (ipa42( 107,j),j=1,7)/  1, 1, 2, 1, 0, 1, 3/
      data (ipa42( 108,j),j=1,7)/  1, 1, 2, 2, 0, 0, 2/
      data (ipa42( 109,j),j=1,7)/  1, 1, 3, 0, 0, 1, 3/
      data (ipa42( 110,j),j=1,7)/  1, 1, 3, 1, 0, 0, 2/
      data (ipa42( 111,j),j=1,7)/  1, 1, 4, 0, 0, 0, 2/
      data (ipa42( 112,j),j=1,7)/  2, 2, 1, 1, 0, 0, 2/
      data (ipa42( 113,j),j=1,7)/  1, 2, 2, 0, 1, 0, 3/
      data (ipa42( 114,j),j=1,7)/  1, 2, 3, 0, 0, 0, 3/
      data (ipa42( 115,j),j=1,7)/  2, 2, 2, 0, 0, 0, 1/
*********** 6 66*******************************
      data (ipa42( 116,j),j=1,7)/  0, 0, 1, 0, 1, 5, 3/
      data (ipa42( 117,j),j=1,7)/  0, 0, 1, 0, 2, 4, 3/
      data (ipa42( 118,j),j=1,7)/  0, 0, 1, 0, 3, 3, 2/
      data (ipa42( 119,j),j=1,7)/  0, 0, 1, 1, 0, 5, 3/
      data (ipa42( 120,j),j=1,7)/  0, 0, 1, 1, 1, 4, 3/
      data (ipa42( 121,j),j=1,7)/  0, 0, 1, 1, 2, 3, 3/
      data (ipa42( 122,j),j=1,7)/  0, 0, 1, 2, 0, 4, 3/
      data (ipa42( 123,j),j=1,7)/  0, 0, 1, 2, 1, 3, 3/
      data (ipa42( 124,j),j=1,7)/  0, 0, 1, 2, 2, 2, 2/
      data (ipa42( 125,j),j=1,7)/  0, 0, 1, 3, 0, 3, 3/
      data (ipa42( 126,j),j=1,7)/  0, 0, 1, 3, 1, 2, 3/
      data (ipa42( 127,j),j=1,7)/  0, 0, 1, 4, 0, 2, 3/
      data (ipa42( 128,j),j=1,7)/  0, 0, 1, 4, 1, 1, 2/
      data (ipa42( 129,j),j=1,7)/  0, 0, 1, 5, 0, 1, 3/
      data (ipa42( 130,j),j=1,7)/  0, 0, 2, 0, 1, 4, 3/
      data (ipa42( 131,j),j=1,7)/  0, 0, 2, 0, 2, 3, 3/
      data (ipa42( 132,j),j=1,7)/  0, 0, 2, 1, 0, 4, 3/
      data (ipa42( 133,j),j=1,7)/  0, 0, 2, 1, 1, 3, 3/
      data (ipa42( 134,j),j=1,7)/  0, 0, 2, 1, 2, 2, 2/
      data (ipa42( 135,j),j=1,7)/  0, 0, 2, 2, 0, 3, 3/
      data (ipa42( 136,j),j=1,7)/  0, 0, 2, 2, 1, 2, 3/
      data (ipa42( 137,j),j=1,7)/  0, 0, 2, 3, 0, 2, 3/
      data (ipa42( 138,j),j=1,7)/  0, 0, 2, 3, 1, 1, 2/
      data (ipa42( 139,j),j=1,7)/  0, 0, 2, 4, 0, 1, 3/
      data (ipa42( 140,j),j=1,7)/  0, 0, 3, 0, 1, 3, 3/
      data (ipa42( 141,j),j=1,7)/  0, 0, 3, 0, 2, 2, 2/
      data (ipa42( 142,j),j=1,7)/  0, 0, 3, 1, 0, 3, 3/
      data (ipa42( 143,j),j=1,7)/  0, 0, 3, 1, 1, 2, 3/
      data (ipa42( 144,j),j=1,7)/  0, 0, 3, 2, 0, 2, 3/
      data (ipa42( 145,j),j=1,7)/  0, 0, 3, 2, 1, 1, 2/
      data (ipa42( 146,j),j=1,7)/  0, 0, 3, 3, 0, 1, 3/
      data (ipa42( 147,j),j=1,7)/  0, 0, 4, 0, 1, 2, 3/
      data (ipa42( 148,j),j=1,7)/  0, 0, 4, 1, 0, 2, 3/
      data (ipa42( 149,j),j=1,7)/  0, 0, 4, 1, 1, 1, 2/
      data (ipa42( 150,j),j=1,7)/  0, 0, 4, 2, 0, 1, 3/
      data (ipa42( 151,j),j=1,7)/  0, 0, 5, 0, 1, 1, 2/
      data (ipa42( 152,j),j=1,7)/  0, 0, 5, 1, 0, 1, 3/
      data (ipa42( 153,j),j=1,7)/  0, 1, 1, 0, 1, 4, 3/
      data (ipa42( 154,j),j=1,7)/  0, 1, 1, 0, 2, 3, 3/
      data (ipa42( 155,j),j=1,7)/  0, 1, 1, 0, 3, 2, 3/
      data (ipa42( 156,j),j=1,7)/  0, 1, 1, 0, 4, 1, 3/
      data (ipa42( 157,j),j=1,7)/  0, 1, 1, 0, 5, 0, 3/
      data (ipa42( 158,j),j=1,7)/  1, 1, 0, 3, 1, 1, 2/
      data (ipa42( 159,j),j=1,7)/  0, 1, 1, 1, 2, 2, 3/
      data (ipa42( 160,j),j=1,7)/  0, 1, 1, 1, 3, 1, 3/
      data (ipa42( 161,j),j=1,7)/  0, 1, 1, 1, 4, 0, 3/
      data (ipa42( 162,j),j=1,7)/  1, 1, 0, 1, 2, 2, 2/
      data (ipa42( 163,j),j=1,7)/  0, 1, 1, 2, 3, 0, 3/
      data (ipa42( 164,j),j=1,7)/  0, 1, 2, 0, 1, 3, 3/
      data (ipa42( 165,j),j=1,7)/  0, 1, 2, 0, 2, 2, 3/
      data (ipa42( 166,j),j=1,7)/  0, 1, 2, 0, 3, 1, 3/
      data (ipa42( 167,j),j=1,7)/  0, 1, 2, 0, 4, 0, 3/
      data (ipa42( 168,j),j=1,7)/  0, 1, 2, 1, 0, 3, 3/
      data (ipa42( 169,j),j=1,7)/  0, 1, 2, 1, 1, 2, 3/
      data (ipa42( 170,j),j=1,7)/  0, 1, 2, 1, 2, 1, 3/
      data (ipa42( 171,j),j=1,7)/  0, 1, 2, 1, 3, 0, 3/
      data (ipa42( 172,j),j=1,7)/  0, 1, 2, 2, 0, 2, 3/
      data (ipa42( 173,j),j=1,7)/  0, 1, 2, 2, 1, 1, 3/
      data (ipa42( 174,j),j=1,7)/  0, 1, 2, 2, 2, 0, 3/
      data (ipa42( 175,j),j=1,7)/  0, 1, 2, 3, 0, 1, 3/
      data (ipa42( 176,j),j=1,7)/  0, 1, 2, 3, 1, 0, 3/
      data (ipa42( 177,j),j=1,7)/  0, 1, 2, 4, 0, 0, 3/
      data (ipa42( 178,j),j=1,7)/  0, 1, 3, 0, 1, 2, 3/
      data (ipa42( 179,j),j=1,7)/  0, 1, 3, 0, 2, 1, 3/
      data (ipa42( 180,j),j=1,7)/  0, 1, 3, 0, 3, 0, 3/
      data (ipa42( 181,j),j=1,7)/  0, 1, 3, 1, 0, 2, 3/
      data (ipa42( 182,j),j=1,7)/  0, 1, 3, 1, 1, 1, 3/
      data (ipa42( 183,j),j=1,7)/  0, 1, 3, 1, 2, 0, 3/
      data (ipa42( 184,j),j=1,7)/  0, 1, 3, 2, 0, 1, 3/
      data (ipa42( 185,j),j=1,7)/  0, 1, 3, 2, 1, 0, 3/
      data (ipa42( 186,j),j=1,7)/  0, 1, 3, 3, 0, 0, 3/
      data (ipa42( 187,j),j=1,7)/  0, 1, 4, 0, 1, 1, 3/
      data (ipa42( 188,j),j=1,7)/  0, 1, 4, 0, 2, 0, 3/
      data (ipa42( 189,j),j=1,7)/  0, 1, 4, 1, 0, 1, 3/
      data (ipa42( 190,j),j=1,7)/  0, 1, 4, 1, 1, 0, 3/
      data (ipa42( 191,j),j=1,7)/  0, 1, 4, 2, 0, 0, 3/
      data (ipa42( 192,j),j=1,7)/  0, 1, 5, 0, 1, 0, 3/
      data (ipa42( 193,j),j=1,7)/  0, 1, 5, 1, 0, 0, 3/
      data (ipa42( 194,j),j=1,7)/  0, 2, 2, 0, 1, 2, 3/
      data (ipa42( 195,j),j=1,7)/  0, 2, 2, 0, 2, 1, 3/
      data (ipa42( 196,j),j=1,7)/  0, 2, 2, 0, 3, 0, 3/
      data (ipa42( 197,j),j=1,7)/  2, 2, 0, 1, 1, 1, 2/
      data (ipa42( 198,j),j=1,7)/  0, 2, 2, 1, 2, 0, 3/
      data (ipa42( 199,j),j=1,7)/  0, 2, 3, 0, 1, 1, 3/
      data (ipa42( 200,j),j=1,7)/  0, 2, 3, 0, 2, 0, 3/
      data (ipa42( 201,j),j=1,7)/  0, 2, 3, 1, 0, 1, 3/
      data (ipa42( 202,j),j=1,7)/  0, 2, 3, 1, 1, 0, 3/
      data (ipa42( 203,j),j=1,7)/  0, 2, 3, 2, 0, 0, 3/
      data (ipa42( 204,j),j=1,7)/  0, 2, 4, 0, 1, 0, 3/
      data (ipa42( 205,j),j=1,7)/  0, 2, 4, 1, 0, 0, 3/
      data (ipa42( 206,j),j=1,7)/  0, 3, 3, 0, 1, 0, 3/
      data (ipa42( 207,j),j=1,7)/  1, 1, 1, 4, 0, 0, 2/
      data (ipa42( 208,j),j=1,7)/  1, 1, 1, 0, 1, 3, 3/
      data (ipa42( 209,j),j=1,7)/  1, 1, 1, 0, 2, 2, 2/
      data (ipa42( 210,j),j=1,7)/  1, 1, 1, 2, 1, 1, 2/
      data (ipa42( 211,j),j=1,7)/  1, 1, 2, 0, 0, 3, 3/
      data (ipa42( 212,j),j=1,7)/  1, 1, 2, 0, 1, 2, 3/
      data (ipa42( 213,j),j=1,7)/  1, 1, 2, 1, 0, 2, 3/
      data (ipa42( 214,j),j=1,7)/  1, 1, 2, 1, 1, 1, 2/
      data (ipa42( 215,j),j=1,7)/  1, 1, 2, 2, 0, 1, 3/
      data (ipa42( 216,j),j=1,7)/  1, 1, 2, 3, 0, 0, 2/
      data (ipa42( 217,j),j=1,7)/  1, 1, 3, 0, 0, 2, 3/
      data (ipa42( 218,j),j=1,7)/  1, 1, 3, 0, 1, 1, 2/
      data (ipa42( 219,j),j=1,7)/  1, 1, 3, 1, 0, 1, 3/
      data (ipa42( 220,j),j=1,7)/  1, 1, 3, 2, 0, 0, 2/
      data (ipa42( 221,j),j=1,7)/  1, 1, 4, 0, 0, 1, 3/
      data (ipa42( 222,j),j=1,7)/  1, 1, 4, 1, 0, 0, 2/
      data (ipa42( 223,j),j=1,7)/  1, 1, 5, 0, 0, 0, 2/
      data (ipa42( 224,j),j=1,7)/  2, 2, 1, 2, 0, 0, 2/
      data (ipa42( 225,j),j=1,7)/  1, 2, 2, 0, 1, 1, 3/
      data (ipa42( 226,j),j=1,7)/  1, 2, 2, 0, 2, 0, 3/
      data (ipa42( 227,j),j=1,7)/  2, 2, 1, 0, 1, 1, 2/
      data (ipa42( 228,j),j=1,7)/  1, 2, 3, 0, 0, 1, 3/
      data (ipa42( 229,j),j=1,7)/  1, 2, 3, 0, 1, 0, 3/
      data (ipa42( 230,j),j=1,7)/  1, 2, 3, 1, 0, 0, 3/
      data (ipa42( 231,j),j=1,7)/  1, 2, 4, 0, 0, 0, 3/
      data (ipa42( 232,j),j=1,7)/  3, 3, 1, 0, 0, 0, 2/
      data (ipa42( 233,j),j=1,7)/  2, 2, 2, 1, 0, 0, 2/
      data (ipa42( 234,j),j=1,7)/  2, 2, 3, 0, 0, 0, 2/
*********** 7 119*******************************
      data (ipa42( 235,j),j=1,7)/  0, 0, 1, 0, 1, 6, 3/
      data (ipa42( 236,j),j=1,7)/  0, 0, 1, 0, 2, 5, 3/
      data (ipa42( 237,j),j=1,7)/  0, 0, 1, 0, 3, 4, 3/
      data (ipa42( 238,j),j=1,7)/  0, 0, 1, 1, 0, 6, 3/
      data (ipa42( 239,j),j=1,7)/  0, 0, 1, 1, 1, 5, 3/
      data (ipa42( 240,j),j=1,7)/  0, 0, 1, 1, 2, 4, 3/
      data (ipa42( 241,j),j=1,7)/  0, 0, 1, 1, 3, 3, 2/
      data (ipa42( 242,j),j=1,7)/  0, 0, 1, 2, 0, 5, 3/
      data (ipa42( 243,j),j=1,7)/  0, 0, 1, 2, 1, 4, 3/
      data (ipa42( 244,j),j=1,7)/  0, 0, 1, 2, 2, 3, 3/
      data (ipa42( 245,j),j=1,7)/  0, 0, 1, 3, 0, 4, 3/
      data (ipa42( 246,j),j=1,7)/  0, 0, 1, 3, 1, 3, 3/
      data (ipa42( 247,j),j=1,7)/  0, 0, 1, 3, 2, 2, 2/
      data (ipa42( 248,j),j=1,7)/  0, 0, 1, 4, 0, 3, 3/
      data (ipa42( 249,j),j=1,7)/  0, 0, 1, 4, 1, 2, 3/
      data (ipa42( 250,j),j=1,7)/  0, 0, 1, 5, 0, 2, 3/
      data (ipa42( 251,j),j=1,7)/  0, 0, 1, 5, 1, 1, 2/
      data (ipa42( 252,j),j=1,7)/  0, 0, 1, 6, 0, 1, 3/
      data (ipa42( 253,j),j=1,7)/  0, 0, 2, 0, 1, 5, 3/
      data (ipa42( 254,j),j=1,7)/  0, 0, 2, 0, 2, 4, 3/
      data (ipa42( 255,j),j=1,7)/  0, 0, 2, 0, 3, 3, 2/
      data (ipa42( 256,j),j=1,7)/  0, 0, 2, 1, 0, 5, 3/
      data (ipa42( 257,j),j=1,7)/  0, 0, 2, 1, 1, 4, 3/
      data (ipa42( 258,j),j=1,7)/  0, 0, 2, 1, 2, 3, 3/
      data (ipa42( 259,j),j=1,7)/  0, 0, 2, 2, 0, 4, 3/
      data (ipa42( 260,j),j=1,7)/  0, 0, 2, 2, 1, 3, 3/
      data (ipa42( 261,j),j=1,7)/  0, 0, 2, 2, 2, 2, 2/
      data (ipa42( 262,j),j=1,7)/  0, 0, 2, 3, 0, 3, 3/
      data (ipa42( 263,j),j=1,7)/  0, 0, 2, 3, 1, 2, 3/
      data (ipa42( 264,j),j=1,7)/  0, 0, 2, 4, 0, 2, 3/
      data (ipa42( 265,j),j=1,7)/  0, 0, 2, 4, 1, 1, 2/
      data (ipa42( 266,j),j=1,7)/  0, 0, 2, 5, 0, 1, 3/
      data (ipa42( 267,j),j=1,7)/  0, 0, 3, 0, 1, 4, 3/
      data (ipa42( 268,j),j=1,7)/  0, 0, 3, 0, 2, 3, 3/
      data (ipa42( 269,j),j=1,7)/  0, 0, 3, 1, 0, 4, 3/
      data (ipa42( 270,j),j=1,7)/  0, 0, 3, 1, 1, 3, 3/
      data (ipa42( 271,j),j=1,7)/  0, 0, 3, 1, 2, 2, 2/
      data (ipa42( 272,j),j=1,7)/  0, 0, 3, 2, 0, 3, 3/
      data (ipa42( 273,j),j=1,7)/  0, 0, 3, 2, 1, 2, 3/
      data (ipa42( 274,j),j=1,7)/  0, 0, 3, 3, 0, 2, 3/
      data (ipa42( 275,j),j=1,7)/  0, 0, 3, 3, 1, 1, 2/
      data (ipa42( 276,j),j=1,7)/  0, 0, 3, 4, 0, 1, 3/
      data (ipa42( 277,j),j=1,7)/  0, 0, 4, 0, 1, 3, 3/
      data (ipa42( 278,j),j=1,7)/  0, 0, 4, 0, 2, 2, 2/
      data (ipa42( 279,j),j=1,7)/  0, 0, 4, 1, 0, 3, 3/
      data (ipa42( 280,j),j=1,7)/  0, 0, 4, 1, 1, 2, 3/
      data (ipa42( 281,j),j=1,7)/  0, 0, 4, 2, 0, 2, 3/
      data (ipa42( 282,j),j=1,7)/  0, 0, 4, 2, 1, 1, 2/
      data (ipa42( 283,j),j=1,7)/  0, 0, 4, 3, 0, 1, 3/
      data (ipa42( 284,j),j=1,7)/  0, 0, 5, 0, 1, 2, 3/
      data (ipa42( 285,j),j=1,7)/  0, 0, 5, 1, 0, 2, 3/
      data (ipa42( 286,j),j=1,7)/  0, 0, 5, 1, 1, 1, 2/
      data (ipa42( 287,j),j=1,7)/  0, 0, 5, 2, 0, 1, 3/
      data (ipa42( 288,j),j=1,7)/  0, 0, 6, 0, 1, 1, 2/
      data (ipa42( 289,j),j=1,7)/  0, 0, 6, 1, 0, 1, 3/
      data (ipa42( 290,j),j=1,7)/  0, 1, 1, 0, 1, 5, 3/
      data (ipa42( 291,j),j=1,7)/  0, 1, 1, 0, 2, 4, 3/
      data (ipa42( 292,j),j=1,7)/  0, 1, 1, 0, 3, 3, 3/
      data (ipa42( 293,j),j=1,7)/  0, 1, 1, 0, 4, 2, 3/
      data (ipa42( 294,j),j=1,7)/  0, 1, 1, 0, 5, 1, 3/
      data (ipa42( 295,j),j=1,7)/  0, 1, 1, 0, 6, 0, 3/
      data (ipa42( 296,j),j=1,7)/  1, 1, 0, 4, 1, 1, 2/
      data (ipa42( 297,j),j=1,7)/  0, 1, 1, 1, 2, 3, 3/
      data (ipa42( 298,j),j=1,7)/  0, 1, 1, 1, 3, 2, 3/
      data (ipa42( 299,j),j=1,7)/  0, 1, 1, 1, 4, 1, 3/
      data (ipa42( 300,j),j=1,7)/  0, 1, 1, 1, 5, 0, 3/
      data (ipa42( 301,j),j=1,7)/  1, 1, 0, 2, 2, 2, 2/
      data (ipa42( 302,j),j=1,7)/  0, 1, 1, 2, 3, 1, 3/
      data (ipa42( 303,j),j=1,7)/  0, 1, 1, 2, 4, 0, 3/
      data (ipa42( 304,j),j=1,7)/  1, 1, 0, 0, 3, 3, 2/
      data (ipa42( 305,j),j=1,7)/  0, 1, 2, 0, 1, 4, 3/
      data (ipa42( 306,j),j=1,7)/  0, 1, 2, 0, 2, 3, 3/
      data (ipa42( 307,j),j=1,7)/  0, 1, 2, 0, 3, 2, 3/
      data (ipa42( 308,j),j=1,7)/  0, 1, 2, 0, 4, 1, 3/
      data (ipa42( 309,j),j=1,7)/  0, 1, 2, 0, 5, 0, 3/
      data (ipa42( 310,j),j=1,7)/  0, 1, 2, 1, 0, 4, 3/
      data (ipa42( 311,j),j=1,7)/  0, 1, 2, 1, 1, 3, 3/
      data (ipa42( 312,j),j=1,7)/  0, 1, 2, 1, 2, 2, 3/
      data (ipa42( 313,j),j=1,7)/  0, 1, 2, 1, 3, 1, 3/
      data (ipa42( 314,j),j=1,7)/  0, 1, 2, 1, 4, 0, 3/
      data (ipa42( 315,j),j=1,7)/  0, 1, 2, 2, 0, 3, 3/
      data (ipa42( 316,j),j=1,7)/  0, 1, 2, 2, 1, 2, 3/
      data (ipa42( 317,j),j=1,7)/  0, 1, 2, 2, 2, 1, 3/
      data (ipa42( 318,j),j=1,7)/  0, 1, 2, 2, 3, 0, 3/
      data (ipa42( 319,j),j=1,7)/  0, 1, 2, 3, 0, 2, 3/
      data (ipa42( 320,j),j=1,7)/  0, 1, 2, 3, 1, 1, 3/
      data (ipa42( 321,j),j=1,7)/  0, 1, 2, 3, 2, 0, 3/
      data (ipa42( 322,j),j=1,7)/  0, 1, 2, 4, 0, 1, 3/
      data (ipa42( 323,j),j=1,7)/  0, 1, 2, 4, 1, 0, 3/
      data (ipa42( 324,j),j=1,7)/  0, 1, 2, 5, 0, 0, 3/
      data (ipa42( 325,j),j=1,7)/  0, 1, 3, 0, 1, 3, 3/
      data (ipa42( 326,j),j=1,7)/  0, 1, 3, 0, 2, 2, 3/
      data (ipa42( 327,j),j=1,7)/  0, 1, 3, 0, 3, 1, 3/
      data (ipa42( 328,j),j=1,7)/  0, 1, 3, 0, 4, 0, 3/
      data (ipa42( 329,j),j=1,7)/  0, 1, 3, 1, 0, 3, 3/
      data (ipa42( 330,j),j=1,7)/  0, 1, 3, 1, 1, 2, 3/
      data (ipa42( 331,j),j=1,7)/  0, 1, 3, 1, 2, 1, 3/
      data (ipa42( 332,j),j=1,7)/  0, 1, 3, 1, 3, 0, 3/
      data (ipa42( 333,j),j=1,7)/  0, 1, 3, 2, 0, 2, 3/
      data (ipa42( 334,j),j=1,7)/  0, 1, 3, 2, 1, 1, 3/
      data (ipa42( 335,j),j=1,7)/  0, 1, 3, 2, 2, 0, 3/
      data (ipa42( 336,j),j=1,7)/  0, 1, 3, 3, 0, 1, 3/
      data (ipa42( 337,j),j=1,7)/  0, 1, 3, 3, 1, 0, 3/
      data (ipa42( 338,j),j=1,7)/  0, 1, 3, 4, 0, 0, 3/
      data (ipa42( 339,j),j=1,7)/  0, 1, 4, 0, 1, 2, 3/
      data (ipa42( 340,j),j=1,7)/  0, 1, 4, 0, 2, 1, 3/
      data (ipa42( 341,j),j=1,7)/  0, 1, 4, 0, 3, 0, 3/
      data (ipa42( 342,j),j=1,7)/  0, 1, 4, 1, 0, 2, 3/
      data (ipa42( 343,j),j=1,7)/  0, 1, 4, 1, 1, 1, 3/
      data (ipa42( 344,j),j=1,7)/  0, 1, 4, 1, 2, 0, 3/
      data (ipa42( 345,j),j=1,7)/  0, 1, 4, 2, 0, 1, 3/
      data (ipa42( 346,j),j=1,7)/  0, 1, 4, 2, 1, 0, 3/
      data (ipa42( 347,j),j=1,7)/  0, 1, 4, 3, 0, 0, 3/
      data (ipa42( 348,j),j=1,7)/  0, 1, 5, 0, 1, 1, 3/
      data (ipa42( 349,j),j=1,7)/  0, 1, 5, 0, 2, 0, 3/
      data (ipa42( 350,j),j=1,7)/  0, 1, 5, 1, 0, 1, 3/
      data (ipa42( 351,j),j=1,7)/  0, 1, 5, 1, 1, 0, 3/
      data (ipa42( 352,j),j=1,7)/  0, 1, 5, 2, 0, 0, 3/
      data (ipa42( 353,j),j=1,7)/  0, 1, 6, 0, 1, 0, 3/
      data (ipa42( 354,j),j=1,7)/  0, 1, 6, 1, 0, 0, 3/
      data (ipa42( 355,j),j=1,7)/  0, 2, 2, 0, 1, 3, 3/
      data (ipa42( 356,j),j=1,7)/  0, 2, 2, 0, 2, 2, 3/
      data (ipa42( 357,j),j=1,7)/  0, 2, 2, 0, 3, 1, 3/
      data (ipa42( 358,j),j=1,7)/  0, 2, 2, 0, 4, 0, 3/
      data (ipa42( 359,j),j=1,7)/  2, 2, 0, 2, 1, 1, 2/
      data (ipa42( 360,j),j=1,7)/  0, 2, 2, 1, 2, 1, 3/
      data (ipa42( 361,j),j=1,7)/  0, 2, 2, 1, 3, 0, 3/
      data (ipa42( 362,j),j=1,7)/  2, 2, 0, 0, 2, 2, 2/
      data (ipa42( 363,j),j=1,7)/  0, 2, 3, 0, 1, 2, 3/
      data (ipa42( 364,j),j=1,7)/  0, 2, 3, 0, 2, 1, 3/
      data (ipa42( 365,j),j=1,7)/  0, 2, 3, 0, 3, 0, 3/
      data (ipa42( 366,j),j=1,7)/  0, 2, 3, 1, 0, 2, 3/
      data (ipa42( 367,j),j=1,7)/  0, 2, 3, 1, 1, 1, 3/
      data (ipa42( 368,j),j=1,7)/  0, 2, 3, 1, 2, 0, 3/
      data (ipa42( 369,j),j=1,7)/  0, 2, 3, 2, 0, 1, 3/
      data (ipa42( 370,j),j=1,7)/  0, 2, 3, 2, 1, 0, 3/
      data (ipa42( 371,j),j=1,7)/  0, 2, 3, 3, 0, 0, 3/
      data (ipa42( 372,j),j=1,7)/  0, 2, 4, 0, 1, 1, 3/
      data (ipa42( 373,j),j=1,7)/  0, 2, 4, 0, 2, 0, 3/
      data (ipa42( 374,j),j=1,7)/  0, 2, 4, 1, 0, 1, 3/
      data (ipa42( 375,j),j=1,7)/  0, 2, 4, 1, 1, 0, 3/
      data (ipa42( 376,j),j=1,7)/  0, 2, 4, 2, 0, 0, 3/
      data (ipa42( 377,j),j=1,7)/  0, 2, 5, 0, 1, 0, 3/
      data (ipa42( 378,j),j=1,7)/  0, 2, 5, 1, 0, 0, 3/
      data (ipa42( 379,j),j=1,7)/  0, 3, 3, 0, 1, 1, 3/
      data (ipa42( 380,j),j=1,7)/  0, 3, 3, 0, 2, 0, 3/
      data (ipa42( 381,j),j=1,7)/  3, 3, 0, 0, 1, 1, 2/
      data (ipa42( 382,j),j=1,7)/  0, 3, 4, 0, 1, 0, 3/
      data (ipa42( 383,j),j=1,7)/  0, 3, 4, 1, 0, 0, 3/
      data (ipa42( 384,j),j=1,7)/  1, 1, 1, 5, 0, 0, 2/
      data (ipa42( 385,j),j=1,7)/  1, 1, 1, 0, 1, 4, 3/
      data (ipa42( 386,j),j=1,7)/  1, 1, 1, 0, 2, 3, 3/
      data (ipa42( 387,j),j=1,7)/  1, 1, 1, 3, 1, 1, 2/
      data (ipa42( 388,j),j=1,7)/  1, 1, 1, 1, 2, 2, 2/
      data (ipa42( 389,j),j=1,7)/  1, 1, 2, 0, 0, 4, 3/
      data (ipa42( 390,j),j=1,7)/  1, 1, 2, 0, 1, 3, 3/
      data (ipa42( 391,j),j=1,7)/  1, 1, 2, 0, 2, 2, 2/
      data (ipa42( 392,j),j=1,7)/  1, 1, 2, 1, 0, 3, 3/
      data (ipa42( 393,j),j=1,7)/  1, 1, 2, 1, 1, 2, 3/
      data (ipa42( 394,j),j=1,7)/  1, 1, 2, 2, 0, 2, 3/
      data (ipa42( 395,j),j=1,7)/  1, 1, 2, 2, 1, 1, 2/
      data (ipa42( 396,j),j=1,7)/  1, 1, 2, 3, 0, 1, 3/
      data (ipa42( 397,j),j=1,7)/  1, 1, 2, 4, 0, 0, 2/
      data (ipa42( 398,j),j=1,7)/  1, 1, 3, 0, 0, 3, 3/
      data (ipa42( 399,j),j=1,7)/  1, 1, 3, 0, 1, 2, 3/
      data (ipa42( 400,j),j=1,7)/  1, 1, 3, 1, 0, 2, 3/
      data (ipa42( 401,j),j=1,7)/  1, 1, 3, 1, 1, 1, 2/
      data (ipa42( 402,j),j=1,7)/  1, 1, 3, 2, 0, 1, 3/
      data (ipa42( 403,j),j=1,7)/  1, 1, 3, 3, 0, 0, 2/
      data (ipa42( 404,j),j=1,7)/  1, 1, 4, 0, 0, 2, 3/
      data (ipa42( 405,j),j=1,7)/  1, 1, 4, 0, 1, 1, 2/
      data (ipa42( 406,j),j=1,7)/  1, 1, 4, 1, 0, 1, 3/
      data (ipa42( 407,j),j=1,7)/  1, 1, 4, 2, 0, 0, 2/
      data (ipa42( 408,j),j=1,7)/  1, 1, 5, 0, 0, 1, 3/
      data (ipa42( 409,j),j=1,7)/  1, 1, 5, 1, 0, 0, 2/
      data (ipa42( 410,j),j=1,7)/  1, 1, 6, 0, 0, 0, 2/
      data (ipa42( 411,j),j=1,7)/  2, 2, 1, 3, 0, 0, 2/
      data (ipa42( 412,j),j=1,7)/  1, 2, 2, 0, 1, 2, 3/
      data (ipa42( 413,j),j=1,7)/  1, 2, 2, 0, 2, 1, 3/
      data (ipa42( 414,j),j=1,7)/  1, 2, 2, 0, 3, 0, 3/
      data (ipa42( 415,j),j=1,7)/  2, 2, 1, 1, 1, 1, 2/
      data (ipa42( 416,j),j=1,7)/  1, 2, 2, 1, 2, 0, 3/
      data (ipa42( 417,j),j=1,7)/  1, 2, 3, 0, 0, 2, 3/
      data (ipa42( 418,j),j=1,7)/  1, 2, 3, 0, 1, 1, 3/
      data (ipa42( 419,j),j=1,7)/  1, 2, 3, 0, 2, 0, 3/
      data (ipa42( 420,j),j=1,7)/  1, 2, 3, 1, 0, 1, 3/
      data (ipa42( 421,j),j=1,7)/  1, 2, 3, 1, 1, 0, 3/
      data (ipa42( 422,j),j=1,7)/  1, 2, 3, 2, 0, 0, 3/
      data (ipa42( 423,j),j=1,7)/  1, 2, 4, 0, 0, 1, 3/
      data (ipa42( 424,j),j=1,7)/  1, 2, 4, 0, 1, 0, 3/
      data (ipa42( 425,j),j=1,7)/  1, 2, 4, 1, 0, 0, 3/
      data (ipa42( 426,j),j=1,7)/  1, 2, 5, 0, 0, 0, 3/
      data (ipa42( 427,j),j=1,7)/  3, 3, 1, 1, 0, 0, 2/
      data (ipa42( 428,j),j=1,7)/  1, 3, 3, 0, 1, 0, 3/
      data (ipa42( 429,j),j=1,7)/  1, 3, 4, 0, 0, 0, 3/
      data (ipa42( 430,j),j=1,7)/  2, 2, 2, 2, 0, 0, 2/
      data (ipa42( 431,j),j=1,7)/  2, 2, 2, 0, 1, 1, 2/
      data (ipa42( 432,j),j=1,7)/  2, 2, 3, 0, 0, 1, 3/
      data (ipa42( 433,j),j=1,7)/  2, 2, 3, 1, 0, 0, 2/
      data (ipa42( 434,j),j=1,7)/  2, 2, 4, 0, 0, 0, 2/
      data (ipa42( 435,j),j=1,7)/  3, 3, 2, 0, 0, 0, 2/
*********** 8 201*******************************
      data (ipa42( 436,j),j=1,7)/  0, 0, 1, 0, 1, 7, 3/
      data (ipa42( 437,j),j=1,7)/  0, 0, 1, 0, 2, 6, 3/
      data (ipa42( 438,j),j=1,7)/  0, 0, 1, 0, 3, 5, 3/
      data (ipa42( 439,j),j=1,7)/  0, 0, 1, 0, 4, 4, 2/
      data (ipa42( 440,j),j=1,7)/  0, 0, 1, 1, 0, 7, 3/
      data (ipa42( 441,j),j=1,7)/  0, 0, 1, 1, 1, 6, 3/
      data (ipa42( 442,j),j=1,7)/  0, 0, 1, 1, 2, 5, 3/
      data (ipa42( 443,j),j=1,7)/  0, 0, 1, 1, 3, 4, 3/
      data (ipa42( 444,j),j=1,7)/  0, 0, 1, 2, 0, 6, 3/
      data (ipa42( 445,j),j=1,7)/  0, 0, 1, 2, 1, 5, 3/
      data (ipa42( 446,j),j=1,7)/  0, 0, 1, 2, 2, 4, 3/
      data (ipa42( 447,j),j=1,7)/  0, 0, 1, 2, 3, 3, 2/
      data (ipa42( 448,j),j=1,7)/  0, 0, 1, 3, 0, 5, 3/
      data (ipa42( 449,j),j=1,7)/  0, 0, 1, 3, 1, 4, 3/
      data (ipa42( 450,j),j=1,7)/  0, 0, 1, 3, 2, 3, 3/
      data (ipa42( 451,j),j=1,7)/  0, 0, 1, 4, 0, 4, 3/
      data (ipa42( 452,j),j=1,7)/  0, 0, 1, 4, 1, 3, 3/
      data (ipa42( 453,j),j=1,7)/  0, 0, 1, 4, 2, 2, 2/
      data (ipa42( 454,j),j=1,7)/  0, 0, 1, 5, 0, 3, 3/
      data (ipa42( 455,j),j=1,7)/  0, 0, 1, 5, 1, 2, 3/
      data (ipa42( 456,j),j=1,7)/  0, 0, 1, 6, 0, 2, 3/
      data (ipa42( 457,j),j=1,7)/  0, 0, 1, 6, 1, 1, 2/
      data (ipa42( 458,j),j=1,7)/  0, 0, 1, 7, 0, 1, 3/
      data (ipa42( 459,j),j=1,7)/  0, 0, 2, 0, 1, 6, 3/
      data (ipa42( 460,j),j=1,7)/  0, 0, 2, 0, 2, 5, 3/
      data (ipa42( 461,j),j=1,7)/  0, 0, 2, 0, 3, 4, 3/
      data (ipa42( 462,j),j=1,7)/  0, 0, 2, 1, 0, 6, 3/
      data (ipa42( 463,j),j=1,7)/  0, 0, 2, 1, 1, 5, 3/
      data (ipa42( 464,j),j=1,7)/  0, 0, 2, 1, 2, 4, 3/
      data (ipa42( 465,j),j=1,7)/  0, 0, 2, 1, 3, 3, 2/
      data (ipa42( 466,j),j=1,7)/  0, 0, 2, 2, 0, 5, 3/
      data (ipa42( 467,j),j=1,7)/  0, 0, 2, 2, 1, 4, 3/
      data (ipa42( 468,j),j=1,7)/  0, 0, 2, 2, 2, 3, 3/
      data (ipa42( 469,j),j=1,7)/  0, 0, 2, 3, 0, 4, 3/
      data (ipa42( 470,j),j=1,7)/  0, 0, 2, 3, 1, 3, 3/
      data (ipa42( 471,j),j=1,7)/  0, 0, 2, 3, 2, 2, 2/
      data (ipa42( 472,j),j=1,7)/  0, 0, 2, 4, 0, 3, 3/
      data (ipa42( 473,j),j=1,7)/  0, 0, 2, 4, 1, 2, 3/
      data (ipa42( 474,j),j=1,7)/  0, 0, 2, 5, 0, 2, 3/
      data (ipa42( 475,j),j=1,7)/  0, 0, 2, 5, 1, 1, 2/
      data (ipa42( 476,j),j=1,7)/  0, 0, 2, 6, 0, 1, 3/
      data (ipa42( 477,j),j=1,7)/  0, 0, 3, 0, 1, 5, 3/
      data (ipa42( 478,j),j=1,7)/  0, 0, 3, 0, 2, 4, 3/
      data (ipa42( 479,j),j=1,7)/  0, 0, 3, 0, 3, 3, 2/
      data (ipa42( 480,j),j=1,7)/  0, 0, 3, 1, 0, 5, 3/
      data (ipa42( 481,j),j=1,7)/  0, 0, 3, 1, 1, 4, 3/
      data (ipa42( 482,j),j=1,7)/  0, 0, 3, 1, 2, 3, 3/
      data (ipa42( 483,j),j=1,7)/  0, 0, 3, 2, 0, 4, 3/
      data (ipa42( 484,j),j=1,7)/  0, 0, 3, 2, 1, 3, 3/
      data (ipa42( 485,j),j=1,7)/  0, 0, 3, 2, 2, 2, 2/
      data (ipa42( 486,j),j=1,7)/  0, 0, 3, 3, 0, 3, 3/
      data (ipa42( 487,j),j=1,7)/  0, 0, 3, 3, 1, 2, 3/
      data (ipa42( 488,j),j=1,7)/  0, 0, 3, 4, 0, 2, 3/
      data (ipa42( 489,j),j=1,7)/  0, 0, 3, 4, 1, 1, 2/
      data (ipa42( 490,j),j=1,7)/  0, 0, 3, 5, 0, 1, 3/
      data (ipa42( 491,j),j=1,7)/  0, 0, 4, 0, 1, 4, 3/
      data (ipa42( 492,j),j=1,7)/  0, 0, 4, 0, 2, 3, 3/
      data (ipa42( 493,j),j=1,7)/  0, 0, 4, 1, 0, 4, 3/
      data (ipa42( 494,j),j=1,7)/  0, 0, 4, 1, 1, 3, 3/
      data (ipa42( 495,j),j=1,7)/  0, 0, 4, 1, 2, 2, 2/
      data (ipa42( 496,j),j=1,7)/  0, 0, 4, 2, 0, 3, 3/
      data (ipa42( 497,j),j=1,7)/  0, 0, 4, 2, 1, 2, 3/
      data (ipa42( 498,j),j=1,7)/  0, 0, 4, 3, 0, 2, 3/
      data (ipa42( 499,j),j=1,7)/  0, 0, 4, 3, 1, 1, 2/
      data (ipa42( 500,j),j=1,7)/  0, 0, 4, 4, 0, 1, 3/
      data (ipa42( 501,j),j=1,7)/  0, 0, 5, 0, 1, 3, 3/
      data (ipa42( 502,j),j=1,7)/  0, 0, 5, 0, 2, 2, 2/
      data (ipa42( 503,j),j=1,7)/  0, 0, 5, 1, 0, 3, 3/
      data (ipa42( 504,j),j=1,7)/  0, 0, 5, 1, 1, 2, 3/
      data (ipa42( 505,j),j=1,7)/  0, 0, 5, 2, 0, 2, 3/
      data (ipa42( 506,j),j=1,7)/  0, 0, 5, 2, 1, 1, 2/
      data (ipa42( 507,j),j=1,7)/  0, 0, 5, 3, 0, 1, 3/
      data (ipa42( 508,j),j=1,7)/  0, 0, 6, 0, 1, 2, 3/
      data (ipa42( 509,j),j=1,7)/  0, 0, 6, 1, 0, 2, 3/
      data (ipa42( 510,j),j=1,7)/  0, 0, 6, 1, 1, 1, 2/
      data (ipa42( 511,j),j=1,7)/  0, 0, 6, 2, 0, 1, 3/
      data (ipa42( 512,j),j=1,7)/  0, 0, 7, 0, 1, 1, 2/
      data (ipa42( 513,j),j=1,7)/  0, 0, 7, 1, 0, 1, 3/
      data (ipa42( 514,j),j=1,7)/  0, 1, 1, 0, 1, 6, 3/
      data (ipa42( 515,j),j=1,7)/  0, 1, 1, 0, 2, 5, 3/
      data (ipa42( 516,j),j=1,7)/  0, 1, 1, 0, 3, 4, 3/
      data (ipa42( 517,j),j=1,7)/  0, 1, 1, 0, 4, 3, 3/
      data (ipa42( 518,j),j=1,7)/  0, 1, 1, 0, 5, 2, 3/
      data (ipa42( 519,j),j=1,7)/  0, 1, 1, 0, 6, 1, 3/
      data (ipa42( 520,j),j=1,7)/  0, 1, 1, 0, 7, 0, 3/
      data (ipa42( 521,j),j=1,7)/  1, 1, 0, 5, 1, 1, 2/
      data (ipa42( 522,j),j=1,7)/  0, 1, 1, 1, 2, 4, 3/
      data (ipa42( 523,j),j=1,7)/  0, 1, 1, 1, 3, 3, 3/
      data (ipa42( 524,j),j=1,7)/  0, 1, 1, 1, 4, 2, 3/
      data (ipa42( 525,j),j=1,7)/  0, 1, 1, 1, 5, 1, 3/
      data (ipa42( 526,j),j=1,7)/  0, 1, 1, 1, 6, 0, 3/
      data (ipa42( 527,j),j=1,7)/  1, 1, 0, 3, 2, 2, 2/
      data (ipa42( 528,j),j=1,7)/  0, 1, 1, 2, 3, 2, 3/
      data (ipa42( 529,j),j=1,7)/  0, 1, 1, 2, 4, 1, 3/
      data (ipa42( 530,j),j=1,7)/  0, 1, 1, 2, 5, 0, 3/
      data (ipa42( 531,j),j=1,7)/  1, 1, 0, 1, 3, 3, 2/
      data (ipa42( 532,j),j=1,7)/  0, 1, 1, 3, 4, 0, 3/
      data (ipa42( 533,j),j=1,7)/  0, 1, 2, 0, 1, 5, 3/
      data (ipa42( 534,j),j=1,7)/  0, 1, 2, 0, 2, 4, 3/
      data (ipa42( 535,j),j=1,7)/  0, 1, 2, 0, 3, 3, 3/
      data (ipa42( 536,j),j=1,7)/  0, 1, 2, 0, 4, 2, 3/
      data (ipa42( 537,j),j=1,7)/  0, 1, 2, 0, 5, 1, 3/
      data (ipa42( 538,j),j=1,7)/  0, 1, 2, 0, 6, 0, 3/
      data (ipa42( 539,j),j=1,7)/  0, 1, 2, 1, 0, 5, 3/
      data (ipa42( 540,j),j=1,7)/  0, 1, 2, 1, 1, 4, 3/
      data (ipa42( 541,j),j=1,7)/  0, 1, 2, 1, 2, 3, 3/
      data (ipa42( 542,j),j=1,7)/  0, 1, 2, 1, 3, 2, 3/
      data (ipa42( 543,j),j=1,7)/  0, 1, 2, 1, 4, 1, 3/
      data (ipa42( 544,j),j=1,7)/  0, 1, 2, 1, 5, 0, 3/
      data (ipa42( 545,j),j=1,7)/  0, 1, 2, 2, 0, 4, 3/
      data (ipa42( 546,j),j=1,7)/  0, 1, 2, 2, 1, 3, 3/
      data (ipa42( 547,j),j=1,7)/  0, 1, 2, 2, 2, 2, 3/
      data (ipa42( 548,j),j=1,7)/  0, 1, 2, 2, 3, 1, 3/
      data (ipa42( 549,j),j=1,7)/  0, 1, 2, 2, 4, 0, 3/
      data (ipa42( 550,j),j=1,7)/  0, 1, 2, 3, 0, 3, 3/
      data (ipa42( 551,j),j=1,7)/  0, 1, 2, 3, 1, 2, 3/
      data (ipa42( 552,j),j=1,7)/  0, 1, 2, 3, 2, 1, 3/
      data (ipa42( 553,j),j=1,7)/  0, 1, 2, 3, 3, 0, 3/
      data (ipa42( 554,j),j=1,7)/  0, 1, 2, 4, 0, 2, 3/
      data (ipa42( 555,j),j=1,7)/  0, 1, 2, 4, 1, 1, 3/
      data (ipa42( 556,j),j=1,7)/  0, 1, 2, 4, 2, 0, 3/
      data (ipa42( 557,j),j=1,7)/  0, 1, 2, 5, 0, 1, 3/
      data (ipa42( 558,j),j=1,7)/  0, 1, 2, 5, 1, 0, 3/
      data (ipa42( 559,j),j=1,7)/  0, 1, 2, 6, 0, 0, 3/
      data (ipa42( 560,j),j=1,7)/  0, 1, 3, 0, 1, 4, 3/
      data (ipa42( 561,j),j=1,7)/  0, 1, 3, 0, 2, 3, 3/
      data (ipa42( 562,j),j=1,7)/  0, 1, 3, 0, 3, 2, 3/
      data (ipa42( 563,j),j=1,7)/  0, 1, 3, 0, 4, 1, 3/
      data (ipa42( 564,j),j=1,7)/  0, 1, 3, 0, 5, 0, 3/
      data (ipa42( 565,j),j=1,7)/  0, 1, 3, 1, 0, 4, 3/
      data (ipa42( 566,j),j=1,7)/  0, 1, 3, 1, 1, 3, 3/
      data (ipa42( 567,j),j=1,7)/  0, 1, 3, 1, 2, 2, 3/
      data (ipa42( 568,j),j=1,7)/  0, 1, 3, 1, 3, 1, 3/
      data (ipa42( 569,j),j=1,7)/  0, 1, 3, 1, 4, 0, 3/
      data (ipa42( 570,j),j=1,7)/  0, 1, 3, 2, 0, 3, 3/
      data (ipa42( 571,j),j=1,7)/  0, 1, 3, 2, 1, 2, 3/
      data (ipa42( 572,j),j=1,7)/  0, 1, 3, 2, 2, 1, 3/
      data (ipa42( 573,j),j=1,7)/  0, 1, 3, 2, 3, 0, 3/
      data (ipa42( 574,j),j=1,7)/  0, 1, 3, 3, 0, 2, 3/
      data (ipa42( 575,j),j=1,7)/  0, 1, 3, 3, 1, 1, 3/
      data (ipa42( 576,j),j=1,7)/  0, 1, 3, 3, 2, 0, 3/
      data (ipa42( 577,j),j=1,7)/  0, 1, 3, 4, 0, 1, 3/
      data (ipa42( 578,j),j=1,7)/  0, 1, 3, 4, 1, 0, 3/
      data (ipa42( 579,j),j=1,7)/  0, 1, 3, 5, 0, 0, 3/
      data (ipa42( 580,j),j=1,7)/  0, 1, 4, 0, 1, 3, 3/
      data (ipa42( 581,j),j=1,7)/  0, 1, 4, 0, 2, 2, 3/
      data (ipa42( 582,j),j=1,7)/  0, 1, 4, 0, 3, 1, 3/
      data (ipa42( 583,j),j=1,7)/  0, 1, 4, 0, 4, 0, 3/
      data (ipa42( 584,j),j=1,7)/  0, 1, 4, 1, 0, 3, 3/
      data (ipa42( 585,j),j=1,7)/  0, 1, 4, 1, 1, 2, 3/
      data (ipa42( 586,j),j=1,7)/  0, 1, 4, 1, 2, 1, 3/
      data (ipa42( 587,j),j=1,7)/  0, 1, 4, 1, 3, 0, 3/
      data (ipa42( 588,j),j=1,7)/  0, 1, 4, 2, 0, 2, 3/
      data (ipa42( 589,j),j=1,7)/  0, 1, 4, 2, 1, 1, 3/
      data (ipa42( 590,j),j=1,7)/  0, 1, 4, 2, 2, 0, 3/
      data (ipa42( 591,j),j=1,7)/  0, 1, 4, 3, 0, 1, 3/
      data (ipa42( 592,j),j=1,7)/  0, 1, 4, 3, 1, 0, 3/
      data (ipa42( 593,j),j=1,7)/  0, 1, 4, 4, 0, 0, 3/
      data (ipa42( 594,j),j=1,7)/  0, 1, 5, 0, 1, 2, 3/
      data (ipa42( 595,j),j=1,7)/  0, 1, 5, 0, 2, 1, 3/
      data (ipa42( 596,j),j=1,7)/  0, 1, 5, 0, 3, 0, 3/
      data (ipa42( 597,j),j=1,7)/  0, 1, 5, 1, 0, 2, 3/
      data (ipa42( 598,j),j=1,7)/  0, 1, 5, 1, 1, 1, 3/
      data (ipa42( 599,j),j=1,7)/  0, 1, 5, 1, 2, 0, 3/
      data (ipa42( 600,j),j=1,7)/  0, 1, 5, 2, 0, 1, 3/
      data (ipa42( 601,j),j=1,7)/  0, 1, 5, 2, 1, 0, 3/
      data (ipa42( 602,j),j=1,7)/  0, 1, 5, 3, 0, 0, 3/
      data (ipa42( 603,j),j=1,7)/  0, 1, 6, 0, 1, 1, 3/
      data (ipa42( 604,j),j=1,7)/  0, 1, 6, 0, 2, 0, 3/
      data (ipa42( 605,j),j=1,7)/  0, 1, 6, 1, 0, 1, 3/
      data (ipa42( 606,j),j=1,7)/  0, 1, 6, 1, 1, 0, 3/
      data (ipa42( 607,j),j=1,7)/  0, 1, 6, 2, 0, 0, 3/
      data (ipa42( 608,j),j=1,7)/  0, 1, 7, 0, 1, 0, 3/
      data (ipa42( 609,j),j=1,7)/  0, 1, 7, 1, 0, 0, 3/
      data (ipa42( 610,j),j=1,7)/  0, 2, 2, 0, 1, 4, 3/
      data (ipa42( 611,j),j=1,7)/  0, 2, 2, 0, 2, 3, 3/
      data (ipa42( 612,j),j=1,7)/  0, 2, 2, 0, 3, 2, 3/
      data (ipa42( 613,j),j=1,7)/  0, 2, 2, 0, 4, 1, 3/
      data (ipa42( 614,j),j=1,7)/  0, 2, 2, 0, 5, 0, 3/
      data (ipa42( 615,j),j=1,7)/  2, 2, 0, 3, 1, 1, 2/
      data (ipa42( 616,j),j=1,7)/  0, 2, 2, 1, 2, 2, 3/
      data (ipa42( 617,j),j=1,7)/  0, 2, 2, 1, 3, 1, 3/
      data (ipa42( 618,j),j=1,7)/  0, 2, 2, 1, 4, 0, 3/
      data (ipa42( 619,j),j=1,7)/  2, 2, 0, 1, 2, 2, 2/
      data (ipa42( 620,j),j=1,7)/  0, 2, 2, 2, 3, 0, 3/
      data (ipa42( 621,j),j=1,7)/  0, 2, 3, 0, 1, 3, 3/
      data (ipa42( 622,j),j=1,7)/  0, 2, 3, 0, 2, 2, 3/
      data (ipa42( 623,j),j=1,7)/  0, 2, 3, 0, 3, 1, 3/
      data (ipa42( 624,j),j=1,7)/  0, 2, 3, 0, 4, 0, 3/
      data (ipa42( 625,j),j=1,7)/  0, 2, 3, 1, 0, 3, 3/
      data (ipa42( 626,j),j=1,7)/  0, 2, 3, 1, 1, 2, 3/
      data (ipa42( 627,j),j=1,7)/  0, 2, 3, 1, 2, 1, 3/
      data (ipa42( 628,j),j=1,7)/  0, 2, 3, 1, 3, 0, 3/
      data (ipa42( 629,j),j=1,7)/  0, 2, 3, 2, 0, 2, 3/
      data (ipa42( 630,j),j=1,7)/  0, 2, 3, 2, 1, 1, 3/
      data (ipa42( 631,j),j=1,7)/  0, 2, 3, 2, 2, 0, 3/
      data (ipa42( 632,j),j=1,7)/  0, 2, 3, 3, 0, 1, 3/
      data (ipa42( 633,j),j=1,7)/  0, 2, 3, 3, 1, 0, 3/
      data (ipa42( 634,j),j=1,7)/  0, 2, 3, 4, 0, 0, 3/
      data (ipa42( 635,j),j=1,7)/  0, 2, 4, 0, 1, 2, 3/
      data (ipa42( 636,j),j=1,7)/  0, 2, 4, 0, 2, 1, 3/
      data (ipa42( 637,j),j=1,7)/  0, 2, 4, 0, 3, 0, 3/
      data (ipa42( 638,j),j=1,7)/  0, 2, 4, 1, 0, 2, 3/
      data (ipa42( 639,j),j=1,7)/  0, 2, 4, 1, 1, 1, 3/
      data (ipa42( 640,j),j=1,7)/  0, 2, 4, 1, 2, 0, 3/
      data (ipa42( 641,j),j=1,7)/  0, 2, 4, 2, 0, 1, 3/
      data (ipa42( 642,j),j=1,7)/  0, 2, 4, 2, 1, 0, 3/
      data (ipa42( 643,j),j=1,7)/  0, 2, 4, 3, 0, 0, 3/
      data (ipa42( 644,j),j=1,7)/  0, 2, 5, 0, 1, 1, 3/
      data (ipa42( 645,j),j=1,7)/  0, 2, 5, 0, 2, 0, 3/
      data (ipa42( 646,j),j=1,7)/  0, 2, 5, 1, 0, 1, 3/
      data (ipa42( 647,j),j=1,7)/  0, 2, 5, 1, 1, 0, 3/
      data (ipa42( 648,j),j=1,7)/  0, 2, 5, 2, 0, 0, 3/
      data (ipa42( 649,j),j=1,7)/  0, 2, 6, 0, 1, 0, 3/
      data (ipa42( 650,j),j=1,7)/  0, 2, 6, 1, 0, 0, 3/
      data (ipa42( 651,j),j=1,7)/  0, 3, 3, 0, 1, 2, 3/
      data (ipa42( 652,j),j=1,7)/  0, 3, 3, 0, 2, 1, 3/
      data (ipa42( 653,j),j=1,7)/  0, 3, 3, 0, 3, 0, 3/
      data (ipa42( 654,j),j=1,7)/  3, 3, 0, 1, 1, 1, 2/
      data (ipa42( 655,j),j=1,7)/  0, 3, 3, 1, 2, 0, 3/
      data (ipa42( 656,j),j=1,7)/  0, 3, 4, 0, 1, 1, 3/
      data (ipa42( 657,j),j=1,7)/  0, 3, 4, 0, 2, 0, 3/
      data (ipa42( 658,j),j=1,7)/  0, 3, 4, 1, 0, 1, 3/
      data (ipa42( 659,j),j=1,7)/  0, 3, 4, 1, 1, 0, 3/
      data (ipa42( 660,j),j=1,7)/  0, 3, 4, 2, 0, 0, 3/
      data (ipa42( 661,j),j=1,7)/  0, 3, 5, 0, 1, 0, 3/
      data (ipa42( 662,j),j=1,7)/  0, 3, 5, 1, 0, 0, 3/
      data (ipa42( 663,j),j=1,7)/  0, 4, 4, 0, 1, 0, 3/
      data (ipa42( 664,j),j=1,7)/  1, 1, 1, 6, 0, 0, 2/
      data (ipa42( 665,j),j=1,7)/  1, 1, 1, 0, 1, 5, 3/
      data (ipa42( 666,j),j=1,7)/  1, 1, 1, 0, 2, 4, 3/
      data (ipa42( 667,j),j=1,7)/  1, 1, 1, 0, 3, 3, 2/
      data (ipa42( 668,j),j=1,7)/  1, 1, 1, 4, 1, 1, 2/
      data (ipa42( 669,j),j=1,7)/  1, 1, 1, 1, 2, 3, 3/
      data (ipa42( 670,j),j=1,7)/  1, 1, 1, 2, 2, 2, 1/
      data (ipa42( 671,j),j=1,7)/  1, 1, 2, 0, 0, 5, 3/
      data (ipa42( 672,j),j=1,7)/  1, 1, 2, 0, 1, 4, 3/
      data (ipa42( 673,j),j=1,7)/  1, 1, 2, 0, 2, 3, 3/
      data (ipa42( 674,j),j=1,7)/  1, 1, 2, 1, 0, 4, 3/
      data (ipa42( 675,j),j=1,7)/  1, 1, 2, 1, 1, 3, 3/
      data (ipa42( 676,j),j=1,7)/  1, 1, 2, 1, 2, 2, 2/
      data (ipa42( 677,j),j=1,7)/  1, 1, 2, 2, 0, 3, 3/
      data (ipa42( 678,j),j=1,7)/  1, 1, 2, 2, 1, 2, 3/
      data (ipa42( 679,j),j=1,7)/  1, 1, 2, 3, 0, 2, 3/
      data (ipa42( 680,j),j=1,7)/  1, 1, 2, 3, 1, 1, 2/
      data (ipa42( 681,j),j=1,7)/  1, 1, 2, 4, 0, 1, 3/
      data (ipa42( 682,j),j=1,7)/  1, 1, 2, 5, 0, 0, 2/
      data (ipa42( 683,j),j=1,7)/  1, 1, 3, 0, 0, 4, 3/
      data (ipa42( 684,j),j=1,7)/  1, 1, 3, 0, 1, 3, 3/
      data (ipa42( 685,j),j=1,7)/  1, 1, 3, 0, 2, 2, 2/
      data (ipa42( 686,j),j=1,7)/  1, 1, 3, 1, 0, 3, 3/
      data (ipa42( 687,j),j=1,7)/  1, 1, 3, 1, 1, 2, 3/
      data (ipa42( 688,j),j=1,7)/  1, 1, 3, 2, 0, 2, 3/
      data (ipa42( 689,j),j=1,7)/  1, 1, 3, 2, 1, 1, 2/
      data (ipa42( 690,j),j=1,7)/  1, 1, 3, 3, 0, 1, 3/
      data (ipa42( 691,j),j=1,7)/  1, 1, 3, 4, 0, 0, 2/
      data (ipa42( 692,j),j=1,7)/  1, 1, 4, 0, 0, 3, 3/
      data (ipa42( 693,j),j=1,7)/  1, 1, 4, 0, 1, 2, 3/
      data (ipa42( 694,j),j=1,7)/  1, 1, 4, 1, 0, 2, 3/
      data (ipa42( 695,j),j=1,7)/  1, 1, 4, 1, 1, 1, 2/
      data (ipa42( 696,j),j=1,7)/  1, 1, 4, 2, 0, 1, 3/
      data (ipa42( 697,j),j=1,7)/  1, 1, 4, 3, 0, 0, 2/
      data (ipa42( 698,j),j=1,7)/  1, 1, 5, 0, 0, 2, 3/
      data (ipa42( 699,j),j=1,7)/  1, 1, 5, 0, 1, 1, 2/
      data (ipa42( 700,j),j=1,7)/  1, 1, 5, 1, 0, 1, 3/
      data (ipa42( 701,j),j=1,7)/  1, 1, 5, 2, 0, 0, 2/
      data (ipa42( 702,j),j=1,7)/  1, 1, 6, 0, 0, 1, 3/
      data (ipa42( 703,j),j=1,7)/  1, 1, 6, 1, 0, 0, 2/
      data (ipa42( 704,j),j=1,7)/  1, 1, 7, 0, 0, 0, 2/
      data (ipa42( 705,j),j=1,7)/  2, 2, 1, 4, 0, 0, 2/
      data (ipa42( 706,j),j=1,7)/  1, 2, 2, 0, 1, 3, 3/
      data (ipa42( 707,j),j=1,7)/  1, 2, 2, 0, 2, 2, 3/
      data (ipa42( 708,j),j=1,7)/  1, 2, 2, 0, 3, 1, 3/
      data (ipa42( 709,j),j=1,7)/  1, 2, 2, 0, 4, 0, 3/
      data (ipa42( 710,j),j=1,7)/  2, 2, 1, 2, 1, 1, 2/
      data (ipa42( 711,j),j=1,7)/  1, 2, 2, 1, 2, 1, 3/
      data (ipa42( 712,j),j=1,7)/  1, 2, 2, 1, 3, 0, 3/
      data (ipa42( 713,j),j=1,7)/  2, 2, 1, 0, 2, 2, 2/
      data (ipa42( 714,j),j=1,7)/  1, 2, 3, 0, 0, 3, 3/
      data (ipa42( 715,j),j=1,7)/  1, 2, 3, 0, 1, 2, 3/
      data (ipa42( 716,j),j=1,7)/  1, 2, 3, 0, 2, 1, 3/
      data (ipa42( 717,j),j=1,7)/  1, 2, 3, 0, 3, 0, 3/
      data (ipa42( 718,j),j=1,7)/  1, 2, 3, 1, 0, 2, 3/
      data (ipa42( 719,j),j=1,7)/  1, 2, 3, 1, 1, 1, 3/
      data (ipa42( 720,j),j=1,7)/  1, 2, 3, 1, 2, 0, 3/
      data (ipa42( 721,j),j=1,7)/  1, 2, 3, 2, 0, 1, 3/
      data (ipa42( 722,j),j=1,7)/  1, 2, 3, 2, 1, 0, 3/
      data (ipa42( 723,j),j=1,7)/  1, 2, 3, 3, 0, 0, 3/
      data (ipa42( 724,j),j=1,7)/  1, 2, 4, 0, 0, 2, 3/
      data (ipa42( 725,j),j=1,7)/  1, 2, 4, 0, 1, 1, 3/
      data (ipa42( 726,j),j=1,7)/  1, 2, 4, 0, 2, 0, 3/
      data (ipa42( 727,j),j=1,7)/  1, 2, 4, 1, 0, 1, 3/
      data (ipa42( 728,j),j=1,7)/  1, 2, 4, 1, 1, 0, 3/
      data (ipa42( 729,j),j=1,7)/  1, 2, 4, 2, 0, 0, 3/
      data (ipa42( 730,j),j=1,7)/  1, 2, 5, 0, 0, 1, 3/
      data (ipa42( 731,j),j=1,7)/  1, 2, 5, 0, 1, 0, 3/
      data (ipa42( 732,j),j=1,7)/  1, 2, 5, 1, 0, 0, 3/
      data (ipa42( 733,j),j=1,7)/  1, 2, 6, 0, 0, 0, 3/
      data (ipa42( 734,j),j=1,7)/  3, 3, 1, 2, 0, 0, 2/
      data (ipa42( 735,j),j=1,7)/  1, 3, 3, 0, 1, 1, 3/
      data (ipa42( 736,j),j=1,7)/  1, 3, 3, 0, 2, 0, 3/
      data (ipa42( 737,j),j=1,7)/  3, 3, 1, 0, 1, 1, 2/
      data (ipa42( 738,j),j=1,7)/  1, 3, 4, 0, 0, 1, 3/
      data (ipa42( 739,j),j=1,7)/  1, 3, 4, 0, 1, 0, 3/
      data (ipa42( 740,j),j=1,7)/  1, 3, 4, 1, 0, 0, 3/
      data (ipa42( 741,j),j=1,7)/  1, 3, 5, 0, 0, 0, 3/
      data (ipa42( 742,j),j=1,7)/  4, 4, 1, 0, 0, 0, 2/
      data (ipa42( 743,j),j=1,7)/  2, 2, 2, 3, 0, 0, 2/
      data (ipa42( 744,j),j=1,7)/  2, 2, 2, 0, 1, 2, 3/
      data (ipa42( 745,j),j=1,7)/  2, 2, 2, 1, 1, 1, 1/
      data (ipa42( 746,j),j=1,7)/  2, 2, 3, 0, 0, 2, 3/
      data (ipa42( 747,j),j=1,7)/  2, 2, 3, 0, 1, 1, 2/
      data (ipa42( 748,j),j=1,7)/  2, 2, 3, 1, 0, 1, 3/
      data (ipa42( 749,j),j=1,7)/  2, 2, 3, 2, 0, 0, 2/
      data (ipa42( 750,j),j=1,7)/  2, 2, 4, 0, 0, 1, 3/
      data (ipa42( 751,j),j=1,7)/  2, 2, 4, 1, 0, 0, 2/
      data (ipa42( 752,j),j=1,7)/  2, 2, 5, 0, 0, 0, 2/
      data (ipa42( 753,j),j=1,7)/  3, 3, 2, 1, 0, 0, 2/
      data (ipa42( 754,j),j=1,7)/  2, 3, 3, 0, 1, 0, 3/
      data (ipa42( 755,j),j=1,7)/  2, 3, 4, 0, 0, 0, 3/
      data (ipa42( 756,j),j=1,7)/  3, 3, 3, 0, 0, 0, 1/
*********** 9 321*******************************
      data (nfas42(i),i=1,9)/0,0,4,17,49,115,234,435,756/
*********** 0 0*******************************
*********** 1 0*******************************
*********** 2 0*******************************
      data (ipa43(   1,j),j=1,7)/  0, 0, 1, 0, 1, 1, 3/
      data (ipa43(   2,j),j=1,7)/  0, 0, 1, 1, 0, 1, 3/
      data (ipa43(   3,j),j=1,7)/  0, 0, 1, 1, 1, 0, 3/
      data (ipa43(   4,j),j=1,7)/  1, 0, 0, 0, 1, 1, 3/
      data (ipa43(   5,j),j=1,7)/  1, 0, 0, 1, 1, 0, 2/
      data (ipa43(   6,j),j=1,7)/  1, 0, 1, 1, 0, 0, 3/
*********** 3 6*******************************
      data (ipa43(   7,j),j=1,7)/  0, 0, 1, 0, 1, 2, 3/
      data (ipa43(   8,j),j=1,7)/  0, 0, 1, 0, 2, 1, 3/
      data (ipa43(   9,j),j=1,7)/  0, 0, 1, 1, 0, 2, 3/
      data (ipa43(  10,j),j=1,7)/  0, 0, 1, 1, 1, 1, 3/
      data (ipa43(  11,j),j=1,7)/  0, 0, 1, 1, 2, 0, 3/
      data (ipa43(  12,j),j=1,7)/  0, 0, 1, 2, 0, 1, 3/
      data (ipa43(  13,j),j=1,7)/  0, 0, 1, 2, 1, 0, 3/
      data (ipa43(  14,j),j=1,7)/  0, 0, 2, 1, 1, 0, 3/
      data (ipa43(  15,j),j=1,7)/  0, 1, 1, 1, 1, 0, 1/
      data (ipa43(  16,j),j=1,7)/  1, 0, 0, 0, 1, 2, 3/
      data (ipa43(  17,j),j=1,7)/  1, 0, 0, 0, 2, 1, 3/
      data (ipa43(  18,j),j=1,7)/  1, 0, 0, 1, 1, 1, 2/
      data (ipa43(  19,j),j=1,7)/  1, 0, 0, 1, 2, 0, 3/
      data (ipa43(  20,j),j=1,7)/  1, 0, 1, 0, 1, 1, 3/
      data (ipa43(  21,j),j=1,7)/  1, 0, 1, 1, 0, 1, 3/
      data (ipa43(  22,j),j=1,7)/  1, 0, 1, 1, 1, 0, 3/
      data (ipa43(  23,j),j=1,7)/  1, 0, 1, 2, 0, 0, 3/
      data (ipa43(  24,j),j=1,7)/  2, 0, 0, 0, 1, 1, 3/
      data (ipa43(  25,j),j=1,7)/  2, 0, 0, 1, 1, 0, 2/
      data (ipa43(  26,j),j=1,7)/  2, 0, 1, 1, 0, 0, 3/
*********** 4 20*******************************
      data (ipa43(  27,j),j=1,7)/  0, 0, 1, 0, 1, 3, 3/
      data (ipa43(  28,j),j=1,7)/  0, 0, 1, 0, 2, 2, 3/
      data (ipa43(  29,j),j=1,7)/  0, 0, 1, 0, 3, 1, 3/
      data (ipa43(  30,j),j=1,7)/  0, 0, 1, 1, 0, 3, 3/
      data (ipa43(  31,j),j=1,7)/  0, 0, 1, 1, 1, 2, 3/
      data (ipa43(  32,j),j=1,7)/  0, 0, 1, 1, 2, 1, 3/
      data (ipa43(  33,j),j=1,7)/  0, 0, 1, 1, 3, 0, 3/
      data (ipa43(  34,j),j=1,7)/  0, 0, 1, 2, 0, 2, 3/
      data (ipa43(  35,j),j=1,7)/  0, 0, 1, 2, 1, 1, 3/
      data (ipa43(  36,j),j=1,7)/  0, 0, 1, 2, 2, 0, 3/
      data (ipa43(  37,j),j=1,7)/  0, 0, 1, 3, 0, 1, 3/
      data (ipa43(  38,j),j=1,7)/  0, 0, 1, 3, 1, 0, 3/
      data (ipa43(  39,j),j=1,7)/  0, 0, 2, 0, 2, 1, 3/
      data (ipa43(  40,j),j=1,7)/  0, 0, 2, 1, 1, 1, 3/
      data (ipa43(  41,j),j=1,7)/  0, 0, 2, 1, 2, 0, 3/
      data (ipa43(  42,j),j=1,7)/  0, 0, 2, 2, 0, 1, 3/
      data (ipa43(  43,j),j=1,7)/  0, 0, 2, 2, 1, 0, 3/
      data (ipa43(  44,j),j=1,7)/  0, 0, 3, 1, 1, 0, 3/
      data (ipa43(  45,j),j=1,7)/  0, 1, 1, 1, 1, 1, 1/
      data (ipa43(  46,j),j=1,7)/  0, 1, 1, 1, 2, 0, 3/
      data (ipa43(  47,j),j=1,7)/  1, 0, 0, 0, 1, 3, 3/
      data (ipa43(  48,j),j=1,7)/  1, 0, 0, 0, 2, 2, 3/
      data (ipa43(  49,j),j=1,7)/  1, 0, 0, 0, 3, 1, 3/
      data (ipa43(  50,j),j=1,7)/  1, 0, 0, 1, 1, 2, 2/
      data (ipa43(  51,j),j=1,7)/  1, 0, 0, 1, 2, 1, 3/
      data (ipa43(  52,j),j=1,7)/  1, 0, 0, 1, 3, 0, 3/
      data (ipa43(  53,j),j=1,7)/  1, 0, 0, 2, 2, 0, 2/
      data (ipa43(  54,j),j=1,7)/  1, 0, 1, 0, 1, 2, 3/
      data (ipa43(  55,j),j=1,7)/  1, 0, 1, 0, 2, 1, 3/
      data (ipa43(  56,j),j=1,7)/  1, 0, 1, 1, 0, 2, 3/
      data (ipa43(  57,j),j=1,7)/  1, 0, 1, 1, 1, 1, 3/
      data (ipa43(  58,j),j=1,7)/  1, 0, 1, 1, 2, 0, 3/
      data (ipa43(  59,j),j=1,7)/  1, 0, 1, 2, 0, 1, 3/
      data (ipa43(  60,j),j=1,7)/  1, 0, 1, 2, 1, 0, 3/
      data (ipa43(  61,j),j=1,7)/  1, 0, 1, 3, 0, 0, 3/
      data (ipa43(  62,j),j=1,7)/  1, 0, 2, 1, 1, 0, 3/
      data (ipa43(  63,j),j=1,7)/  1, 0, 2, 2, 0, 0, 3/
      data (ipa43(  64,j),j=1,7)/  1, 1, 1, 1, 1, 0, 1/
      data (ipa43(  65,j),j=1,7)/  2, 0, 0, 0, 1, 2, 3/
      data (ipa43(  66,j),j=1,7)/  2, 0, 0, 0, 2, 1, 3/
      data (ipa43(  67,j),j=1,7)/  2, 0, 0, 1, 1, 1, 2/
      data (ipa43(  68,j),j=1,7)/  2, 0, 0, 1, 2, 0, 3/
      data (ipa43(  69,j),j=1,7)/  2, 0, 1, 0, 1, 1, 3/
      data (ipa43(  70,j),j=1,7)/  2, 0, 1, 1, 0, 1, 3/
      data (ipa43(  71,j),j=1,7)/  2, 0, 1, 1, 1, 0, 3/
      data (ipa43(  72,j),j=1,7)/  2, 0, 1, 2, 0, 0, 3/
      data (ipa43(  73,j),j=1,7)/  3, 0, 0, 0, 1, 1, 3/
      data (ipa43(  74,j),j=1,7)/  3, 0, 0, 1, 1, 0, 2/
      data (ipa43(  75,j),j=1,7)/  3, 0, 1, 1, 0, 0, 3/
*********** 5 49*******************************
      data (ipa43(  76,j),j=1,7)/  0, 0, 1, 0, 1, 4, 3/
      data (ipa43(  77,j),j=1,7)/  0, 0, 1, 0, 2, 3, 3/
      data (ipa43(  78,j),j=1,7)/  0, 0, 1, 0, 3, 2, 3/
      data (ipa43(  79,j),j=1,7)/  0, 0, 1, 0, 4, 1, 3/
      data (ipa43(  80,j),j=1,7)/  0, 0, 1, 1, 0, 4, 3/
      data (ipa43(  81,j),j=1,7)/  0, 0, 1, 1, 1, 3, 3/
      data (ipa43(  82,j),j=1,7)/  0, 0, 1, 1, 2, 2, 3/
      data (ipa43(  83,j),j=1,7)/  0, 0, 1, 1, 3, 1, 3/
      data (ipa43(  84,j),j=1,7)/  0, 0, 1, 1, 4, 0, 3/
      data (ipa43(  85,j),j=1,7)/  0, 0, 1, 2, 0, 3, 3/
      data (ipa43(  86,j),j=1,7)/  0, 0, 1, 2, 1, 2, 3/
      data (ipa43(  87,j),j=1,7)/  0, 0, 1, 2, 2, 1, 3/
      data (ipa43(  88,j),j=1,7)/  0, 0, 1, 2, 3, 0, 3/
      data (ipa43(  89,j),j=1,7)/  0, 0, 1, 3, 0, 2, 3/
      data (ipa43(  90,j),j=1,7)/  0, 0, 1, 3, 1, 1, 3/
      data (ipa43(  91,j),j=1,7)/  0, 0, 1, 3, 2, 0, 3/
      data (ipa43(  92,j),j=1,7)/  0, 0, 1, 4, 0, 1, 3/
      data (ipa43(  93,j),j=1,7)/  0, 0, 1, 4, 1, 0, 3/
      data (ipa43(  94,j),j=1,7)/  0, 0, 2, 0, 2, 2, 3/
      data (ipa43(  95,j),j=1,7)/  0, 0, 2, 0, 3, 1, 3/
      data (ipa43(  96,j),j=1,7)/  0, 0, 2, 1, 1, 2, 3/
      data (ipa43(  97,j),j=1,7)/  0, 0, 2, 1, 2, 1, 3/
      data (ipa43(  98,j),j=1,7)/  0, 0, 2, 1, 3, 0, 3/
      data (ipa43(  99,j),j=1,7)/  0, 0, 2, 2, 0, 2, 3/
      data (ipa43( 100,j),j=1,7)/  0, 0, 2, 2, 1, 1, 3/
      data (ipa43( 101,j),j=1,7)/  0, 0, 2, 2, 2, 0, 3/
      data (ipa43( 102,j),j=1,7)/  0, 0, 2, 3, 0, 1, 3/
      data (ipa43( 103,j),j=1,7)/  0, 0, 2, 3, 1, 0, 3/
      data (ipa43( 104,j),j=1,7)/  0, 0, 3, 1, 1, 1, 3/
      data (ipa43( 105,j),j=1,7)/  0, 0, 3, 1, 2, 0, 3/
      data (ipa43( 106,j),j=1,7)/  0, 0, 3, 2, 1, 0, 3/
      data (ipa43( 107,j),j=1,7)/  0, 0, 4, 1, 1, 0, 3/
      data (ipa43( 108,j),j=1,7)/  0, 1, 1, 1, 1, 2, 1/
      data (ipa43( 109,j),j=1,7)/  0, 1, 1, 1, 2, 1, 3/
      data (ipa43( 110,j),j=1,7)/  0, 1, 1, 1, 3, 0, 3/
      data (ipa43( 111,j),j=1,7)/  0, 1, 1, 2, 2, 0, 2/
      data (ipa43( 112,j),j=1,7)/  0, 1, 2, 1, 2, 0, 3/
      data (ipa43( 113,j),j=1,7)/  0, 1, 2, 2, 1, 0, 3/
      data (ipa43( 114,j),j=1,7)/  1, 0, 0, 0, 1, 4, 3/
      data (ipa43( 115,j),j=1,7)/  1, 0, 0, 0, 2, 3, 3/
      data (ipa43( 116,j),j=1,7)/  1, 0, 0, 0, 3, 2, 3/
      data (ipa43( 117,j),j=1,7)/  1, 0, 0, 0, 4, 1, 3/
      data (ipa43( 118,j),j=1,7)/  1, 0, 0, 1, 1, 3, 2/
      data (ipa43( 119,j),j=1,7)/  1, 0, 0, 1, 2, 2, 3/
      data (ipa43( 120,j),j=1,7)/  1, 0, 0, 1, 3, 1, 3/
      data (ipa43( 121,j),j=1,7)/  1, 0, 0, 1, 4, 0, 3/
      data (ipa43( 122,j),j=1,7)/  1, 0, 0, 2, 2, 1, 2/
      data (ipa43( 123,j),j=1,7)/  1, 0, 0, 2, 3, 0, 3/
      data (ipa43( 124,j),j=1,7)/  1, 0, 1, 0, 1, 3, 3/
      data (ipa43( 125,j),j=1,7)/  1, 0, 1, 0, 2, 2, 3/
      data (ipa43( 126,j),j=1,7)/  1, 0, 1, 0, 3, 1, 3/
      data (ipa43( 127,j),j=1,7)/  1, 0, 1, 1, 0, 3, 3/
      data (ipa43( 128,j),j=1,7)/  1, 0, 1, 1, 1, 2, 3/
      data (ipa43( 129,j),j=1,7)/  1, 0, 1, 1, 2, 1, 3/
      data (ipa43( 130,j),j=1,7)/  1, 0, 1, 1, 3, 0, 3/
      data (ipa43( 131,j),j=1,7)/  1, 0, 1, 2, 0, 2, 3/
      data (ipa43( 132,j),j=1,7)/  1, 0, 1, 2, 1, 1, 3/
      data (ipa43( 133,j),j=1,7)/  1, 0, 1, 2, 2, 0, 3/
      data (ipa43( 134,j),j=1,7)/  1, 0, 1, 3, 0, 1, 3/
      data (ipa43( 135,j),j=1,7)/  1, 0, 1, 3, 1, 0, 3/
      data (ipa43( 136,j),j=1,7)/  1, 0, 1, 4, 0, 0, 3/
      data (ipa43( 137,j),j=1,7)/  1, 0, 2, 0, 2, 1, 3/
      data (ipa43( 138,j),j=1,7)/  1, 0, 2, 1, 1, 1, 3/
      data (ipa43( 139,j),j=1,7)/  1, 0, 2, 1, 2, 0, 3/
      data (ipa43( 140,j),j=1,7)/  1, 0, 2, 2, 0, 1, 3/
      data (ipa43( 141,j),j=1,7)/  1, 0, 2, 2, 1, 0, 3/
      data (ipa43( 142,j),j=1,7)/  1, 0, 2, 3, 0, 0, 3/
      data (ipa43( 143,j),j=1,7)/  1, 0, 3, 1, 1, 0, 3/
      data (ipa43( 144,j),j=1,7)/  1, 1, 1, 1, 1, 1, 1/
      data (ipa43( 145,j),j=1,7)/  1, 1, 1, 1, 2, 0, 3/
      data (ipa43( 146,j),j=1,7)/  2, 0, 0, 0, 1, 3, 3/
      data (ipa43( 147,j),j=1,7)/  2, 0, 0, 0, 2, 2, 3/
      data (ipa43( 148,j),j=1,7)/  2, 0, 0, 0, 3, 1, 3/
      data (ipa43( 149,j),j=1,7)/  2, 0, 0, 1, 1, 2, 2/
      data (ipa43( 150,j),j=1,7)/  2, 0, 0, 1, 2, 1, 3/
      data (ipa43( 151,j),j=1,7)/  2, 0, 0, 1, 3, 0, 3/
      data (ipa43( 152,j),j=1,7)/  2, 0, 0, 2, 2, 0, 2/
      data (ipa43( 153,j),j=1,7)/  2, 0, 1, 0, 1, 2, 3/
      data (ipa43( 154,j),j=1,7)/  2, 0, 1, 0, 2, 1, 3/
      data (ipa43( 155,j),j=1,7)/  2, 0, 1, 1, 0, 2, 3/
      data (ipa43( 156,j),j=1,7)/  2, 0, 1, 1, 1, 1, 3/
      data (ipa43( 157,j),j=1,7)/  2, 0, 1, 1, 2, 0, 3/
      data (ipa43( 158,j),j=1,7)/  2, 0, 1, 2, 0, 1, 3/
      data (ipa43( 159,j),j=1,7)/  2, 0, 1, 2, 1, 0, 3/
      data (ipa43( 160,j),j=1,7)/  2, 0, 1, 3, 0, 0, 3/
      data (ipa43( 161,j),j=1,7)/  2, 0, 2, 1, 1, 0, 3/
      data (ipa43( 162,j),j=1,7)/  2, 0, 2, 2, 0, 0, 3/
      data (ipa43( 163,j),j=1,7)/  2, 1, 1, 1, 1, 0, 1/
      data (ipa43( 164,j),j=1,7)/  3, 0, 0, 0, 1, 2, 3/
      data (ipa43( 165,j),j=1,7)/  3, 0, 0, 0, 2, 1, 3/
      data (ipa43( 166,j),j=1,7)/  3, 0, 0, 1, 1, 1, 2/
      data (ipa43( 167,j),j=1,7)/  3, 0, 0, 1, 2, 0, 3/
      data (ipa43( 168,j),j=1,7)/  3, 0, 1, 0, 1, 1, 3/
      data (ipa43( 169,j),j=1,7)/  3, 0, 1, 1, 0, 1, 3/
      data (ipa43( 170,j),j=1,7)/  3, 0, 1, 1, 1, 0, 3/
      data (ipa43( 171,j),j=1,7)/  3, 0, 1, 2, 0, 0, 3/
      data (ipa43( 172,j),j=1,7)/  4, 0, 0, 0, 1, 1, 3/
      data (ipa43( 173,j),j=1,7)/  4, 0, 0, 1, 1, 0, 2/
      data (ipa43( 174,j),j=1,7)/  4, 0, 1, 1, 0, 0, 3/
*********** 6 99*******************************
      data (ipa43( 175,j),j=1,7)/  0, 0, 1, 0, 1, 5, 3/
      data (ipa43( 176,j),j=1,7)/  0, 0, 1, 0, 2, 4, 3/
      data (ipa43( 177,j),j=1,7)/  0, 0, 1, 0, 3, 3, 3/
      data (ipa43( 178,j),j=1,7)/  0, 0, 1, 0, 4, 2, 3/
      data (ipa43( 179,j),j=1,7)/  0, 0, 1, 0, 5, 1, 3/
      data (ipa43( 180,j),j=1,7)/  0, 0, 1, 1, 0, 5, 3/
      data (ipa43( 181,j),j=1,7)/  0, 0, 1, 1, 1, 4, 3/
      data (ipa43( 182,j),j=1,7)/  0, 0, 1, 1, 2, 3, 3/
      data (ipa43( 183,j),j=1,7)/  0, 0, 1, 1, 3, 2, 3/
      data (ipa43( 184,j),j=1,7)/  0, 0, 1, 1, 4, 1, 3/
      data (ipa43( 185,j),j=1,7)/  0, 0, 1, 1, 5, 0, 3/
      data (ipa43( 186,j),j=1,7)/  0, 0, 1, 2, 0, 4, 3/
      data (ipa43( 187,j),j=1,7)/  0, 0, 1, 2, 1, 3, 3/
      data (ipa43( 188,j),j=1,7)/  0, 0, 1, 2, 2, 2, 3/
      data (ipa43( 189,j),j=1,7)/  0, 0, 1, 2, 3, 1, 3/
      data (ipa43( 190,j),j=1,7)/  0, 0, 1, 2, 4, 0, 3/
      data (ipa43( 191,j),j=1,7)/  0, 0, 1, 3, 0, 3, 3/
      data (ipa43( 192,j),j=1,7)/  0, 0, 1, 3, 1, 2, 3/
      data (ipa43( 193,j),j=1,7)/  0, 0, 1, 3, 2, 1, 3/
      data (ipa43( 194,j),j=1,7)/  0, 0, 1, 3, 3, 0, 3/
      data (ipa43( 195,j),j=1,7)/  0, 0, 1, 4, 0, 2, 3/
      data (ipa43( 196,j),j=1,7)/  0, 0, 1, 4, 1, 1, 3/
      data (ipa43( 197,j),j=1,7)/  0, 0, 1, 4, 2, 0, 3/
      data (ipa43( 198,j),j=1,7)/  0, 0, 1, 5, 0, 1, 3/
      data (ipa43( 199,j),j=1,7)/  0, 0, 1, 5, 1, 0, 3/
      data (ipa43( 200,j),j=1,7)/  0, 0, 2, 0, 2, 3, 3/
      data (ipa43( 201,j),j=1,7)/  0, 0, 2, 0, 3, 2, 3/
      data (ipa43( 202,j),j=1,7)/  0, 0, 2, 0, 4, 1, 3/
      data (ipa43( 203,j),j=1,7)/  0, 0, 2, 1, 1, 3, 3/
      data (ipa43( 204,j),j=1,7)/  0, 0, 2, 1, 2, 2, 3/
      data (ipa43( 205,j),j=1,7)/  0, 0, 2, 1, 3, 1, 3/
      data (ipa43( 206,j),j=1,7)/  0, 0, 2, 1, 4, 0, 3/
      data (ipa43( 207,j),j=1,7)/  0, 0, 2, 2, 0, 3, 3/
      data (ipa43( 208,j),j=1,7)/  0, 0, 2, 2, 1, 2, 3/
      data (ipa43( 209,j),j=1,7)/  0, 0, 2, 2, 2, 1, 3/
      data (ipa43( 210,j),j=1,7)/  0, 0, 2, 2, 3, 0, 3/
      data (ipa43( 211,j),j=1,7)/  0, 0, 2, 3, 0, 2, 3/
      data (ipa43( 212,j),j=1,7)/  0, 0, 2, 3, 1, 1, 3/
      data (ipa43( 213,j),j=1,7)/  0, 0, 2, 3, 2, 0, 3/
      data (ipa43( 214,j),j=1,7)/  0, 0, 2, 4, 0, 1, 3/
      data (ipa43( 215,j),j=1,7)/  0, 0, 2, 4, 1, 0, 3/
      data (ipa43( 216,j),j=1,7)/  0, 0, 3, 0, 3, 1, 3/
      data (ipa43( 217,j),j=1,7)/  0, 0, 3, 1, 1, 2, 3/
      data (ipa43( 218,j),j=1,7)/  0, 0, 3, 1, 2, 1, 3/
      data (ipa43( 219,j),j=1,7)/  0, 0, 3, 1, 3, 0, 3/
      data (ipa43( 220,j),j=1,7)/  0, 0, 3, 2, 1, 1, 3/
      data (ipa43( 221,j),j=1,7)/  0, 0, 3, 2, 2, 0, 3/
      data (ipa43( 222,j),j=1,7)/  0, 0, 3, 3, 0, 1, 3/
      data (ipa43( 223,j),j=1,7)/  0, 0, 3, 3, 1, 0, 3/
      data (ipa43( 224,j),j=1,7)/  0, 0, 4, 1, 1, 1, 3/
      data (ipa43( 225,j),j=1,7)/  0, 0, 4, 1, 2, 0, 3/
      data (ipa43( 226,j),j=1,7)/  0, 0, 4, 2, 1, 0, 3/
      data (ipa43( 227,j),j=1,7)/  0, 0, 5, 1, 1, 0, 3/
      data (ipa43( 228,j),j=1,7)/  0, 1, 1, 1, 1, 3, 1/
      data (ipa43( 229,j),j=1,7)/  0, 1, 1, 1, 2, 2, 3/
      data (ipa43( 230,j),j=1,7)/  0, 1, 1, 1, 3, 1, 3/
      data (ipa43( 231,j),j=1,7)/  0, 1, 1, 1, 4, 0, 3/
      data (ipa43( 232,j),j=1,7)/  0, 1, 1, 2, 2, 1, 2/
      data (ipa43( 233,j),j=1,7)/  0, 1, 1, 2, 3, 0, 3/
      data (ipa43( 234,j),j=1,7)/  0, 1, 2, 1, 2, 1, 3/
      data (ipa43( 235,j),j=1,7)/  0, 1, 2, 1, 3, 0, 3/
      data (ipa43( 236,j),j=1,7)/  0, 1, 2, 2, 1, 1, 3/
      data (ipa43( 237,j),j=1,7)/  0, 1, 2, 2, 2, 0, 3/
      data (ipa43( 238,j),j=1,7)/  0, 1, 2, 3, 1, 0, 3/
      data (ipa43( 239,j),j=1,7)/  1, 0, 0, 0, 1, 5, 3/
      data (ipa43( 240,j),j=1,7)/  1, 0, 0, 0, 2, 4, 3/
      data (ipa43( 241,j),j=1,7)/  1, 0, 0, 0, 3, 3, 3/
      data (ipa43( 242,j),j=1,7)/  1, 0, 0, 0, 4, 2, 3/
      data (ipa43( 243,j),j=1,7)/  1, 0, 0, 0, 5, 1, 3/
      data (ipa43( 244,j),j=1,7)/  1, 0, 0, 1, 1, 4, 2/
      data (ipa43( 245,j),j=1,7)/  1, 0, 0, 1, 2, 3, 3/
      data (ipa43( 246,j),j=1,7)/  1, 0, 0, 1, 3, 2, 3/
      data (ipa43( 247,j),j=1,7)/  1, 0, 0, 1, 4, 1, 3/
      data (ipa43( 248,j),j=1,7)/  1, 0, 0, 1, 5, 0, 3/
      data (ipa43( 249,j),j=1,7)/  1, 0, 0, 2, 2, 2, 2/
      data (ipa43( 250,j),j=1,7)/  1, 0, 0, 2, 3, 1, 3/
      data (ipa43( 251,j),j=1,7)/  1, 0, 0, 2, 4, 0, 3/
      data (ipa43( 252,j),j=1,7)/  1, 0, 0, 3, 3, 0, 2/
      data (ipa43( 253,j),j=1,7)/  1, 0, 1, 0, 1, 4, 3/
      data (ipa43( 254,j),j=1,7)/  1, 0, 1, 0, 2, 3, 3/
      data (ipa43( 255,j),j=1,7)/  1, 0, 1, 0, 3, 2, 3/
      data (ipa43( 256,j),j=1,7)/  1, 0, 1, 0, 4, 1, 3/
      data (ipa43( 257,j),j=1,7)/  1, 0, 1, 1, 0, 4, 3/
      data (ipa43( 258,j),j=1,7)/  1, 0, 1, 1, 1, 3, 3/
      data (ipa43( 259,j),j=1,7)/  1, 0, 1, 1, 2, 2, 3/
      data (ipa43( 260,j),j=1,7)/  1, 0, 1, 1, 3, 1, 3/
      data (ipa43( 261,j),j=1,7)/  1, 0, 1, 1, 4, 0, 3/
      data (ipa43( 262,j),j=1,7)/  1, 0, 1, 2, 0, 3, 3/
      data (ipa43( 263,j),j=1,7)/  1, 0, 1, 2, 1, 2, 3/
      data (ipa43( 264,j),j=1,7)/  1, 0, 1, 2, 2, 1, 3/
      data (ipa43( 265,j),j=1,7)/  1, 0, 1, 2, 3, 0, 3/
      data (ipa43( 266,j),j=1,7)/  1, 0, 1, 3, 0, 2, 3/
      data (ipa43( 267,j),j=1,7)/  1, 0, 1, 3, 1, 1, 3/
      data (ipa43( 268,j),j=1,7)/  1, 0, 1, 3, 2, 0, 3/
      data (ipa43( 269,j),j=1,7)/  1, 0, 1, 4, 0, 1, 3/
      data (ipa43( 270,j),j=1,7)/  1, 0, 1, 4, 1, 0, 3/
      data (ipa43( 271,j),j=1,7)/  1, 0, 1, 5, 0, 0, 3/
      data (ipa43( 272,j),j=1,7)/  1, 0, 2, 0, 2, 2, 3/
      data (ipa43( 273,j),j=1,7)/  1, 0, 2, 0, 3, 1, 3/
      data (ipa43( 274,j),j=1,7)/  1, 0, 2, 1, 1, 2, 3/
      data (ipa43( 275,j),j=1,7)/  1, 0, 2, 1, 2, 1, 3/
      data (ipa43( 276,j),j=1,7)/  1, 0, 2, 1, 3, 0, 3/
      data (ipa43( 277,j),j=1,7)/  1, 0, 2, 2, 0, 2, 3/
      data (ipa43( 278,j),j=1,7)/  1, 0, 2, 2, 1, 1, 3/
      data (ipa43( 279,j),j=1,7)/  1, 0, 2, 2, 2, 0, 3/
      data (ipa43( 280,j),j=1,7)/  1, 0, 2, 3, 0, 1, 3/
      data (ipa43( 281,j),j=1,7)/  1, 0, 2, 3, 1, 0, 3/
      data (ipa43( 282,j),j=1,7)/  1, 0, 2, 4, 0, 0, 3/
      data (ipa43( 283,j),j=1,7)/  1, 0, 3, 1, 1, 1, 3/
      data (ipa43( 284,j),j=1,7)/  1, 0, 3, 1, 2, 0, 3/
      data (ipa43( 285,j),j=1,7)/  1, 0, 3, 2, 1, 0, 3/
      data (ipa43( 286,j),j=1,7)/  1, 0, 3, 3, 0, 0, 3/
      data (ipa43( 287,j),j=1,7)/  1, 0, 4, 1, 1, 0, 3/
      data (ipa43( 288,j),j=1,7)/  1, 1, 1, 1, 1, 2, 1/
      data (ipa43( 289,j),j=1,7)/  1, 1, 1, 1, 2, 1, 3/
      data (ipa43( 290,j),j=1,7)/  1, 1, 1, 1, 3, 0, 3/
      data (ipa43( 291,j),j=1,7)/  1, 1, 1, 2, 2, 0, 2/
      data (ipa43( 292,j),j=1,7)/  1, 1, 2, 1, 2, 0, 3/
      data (ipa43( 293,j),j=1,7)/  1, 1, 2, 2, 1, 0, 3/
      data (ipa43( 294,j),j=1,7)/  2, 0, 0, 0, 1, 4, 3/
      data (ipa43( 295,j),j=1,7)/  2, 0, 0, 0, 2, 3, 3/
      data (ipa43( 296,j),j=1,7)/  2, 0, 0, 0, 3, 2, 3/
      data (ipa43( 297,j),j=1,7)/  2, 0, 0, 0, 4, 1, 3/
      data (ipa43( 298,j),j=1,7)/  2, 0, 0, 1, 1, 3, 2/
      data (ipa43( 299,j),j=1,7)/  2, 0, 0, 1, 2, 2, 3/
      data (ipa43( 300,j),j=1,7)/  2, 0, 0, 1, 3, 1, 3/
      data (ipa43( 301,j),j=1,7)/  2, 0, 0, 1, 4, 0, 3/
      data (ipa43( 302,j),j=1,7)/  2, 0, 0, 2, 2, 1, 2/
      data (ipa43( 303,j),j=1,7)/  2, 0, 0, 2, 3, 0, 3/
      data (ipa43( 304,j),j=1,7)/  2, 0, 1, 0, 1, 3, 3/
      data (ipa43( 305,j),j=1,7)/  2, 0, 1, 0, 2, 2, 3/
      data (ipa43( 306,j),j=1,7)/  2, 0, 1, 0, 3, 1, 3/
      data (ipa43( 307,j),j=1,7)/  2, 0, 1, 1, 0, 3, 3/
      data (ipa43( 308,j),j=1,7)/  2, 0, 1, 1, 1, 2, 3/
      data (ipa43( 309,j),j=1,7)/  2, 0, 1, 1, 2, 1, 3/
      data (ipa43( 310,j),j=1,7)/  2, 0, 1, 1, 3, 0, 3/
      data (ipa43( 311,j),j=1,7)/  2, 0, 1, 2, 0, 2, 3/
      data (ipa43( 312,j),j=1,7)/  2, 0, 1, 2, 1, 1, 3/
      data (ipa43( 313,j),j=1,7)/  2, 0, 1, 2, 2, 0, 3/
      data (ipa43( 314,j),j=1,7)/  2, 0, 1, 3, 0, 1, 3/
      data (ipa43( 315,j),j=1,7)/  2, 0, 1, 3, 1, 0, 3/
      data (ipa43( 316,j),j=1,7)/  2, 0, 1, 4, 0, 0, 3/
      data (ipa43( 317,j),j=1,7)/  2, 0, 2, 0, 2, 1, 3/
      data (ipa43( 318,j),j=1,7)/  2, 0, 2, 1, 1, 1, 3/
      data (ipa43( 319,j),j=1,7)/  2, 0, 2, 1, 2, 0, 3/
      data (ipa43( 320,j),j=1,7)/  2, 0, 2, 2, 0, 1, 3/
      data (ipa43( 321,j),j=1,7)/  2, 0, 2, 2, 1, 0, 3/
      data (ipa43( 322,j),j=1,7)/  2, 0, 2, 3, 0, 0, 3/
      data (ipa43( 323,j),j=1,7)/  2, 0, 3, 1, 1, 0, 3/
      data (ipa43( 324,j),j=1,7)/  2, 1, 1, 1, 1, 1, 1/
      data (ipa43( 325,j),j=1,7)/  2, 1, 1, 1, 2, 0, 3/
      data (ipa43( 326,j),j=1,7)/  3, 0, 0, 0, 1, 3, 3/
      data (ipa43( 327,j),j=1,7)/  3, 0, 0, 0, 2, 2, 3/
      data (ipa43( 328,j),j=1,7)/  3, 0, 0, 0, 3, 1, 3/
      data (ipa43( 329,j),j=1,7)/  3, 0, 0, 1, 1, 2, 2/
      data (ipa43( 330,j),j=1,7)/  3, 0, 0, 1, 2, 1, 3/
      data (ipa43( 331,j),j=1,7)/  3, 0, 0, 1, 3, 0, 3/
      data (ipa43( 332,j),j=1,7)/  3, 0, 0, 2, 2, 0, 2/
      data (ipa43( 333,j),j=1,7)/  3, 0, 1, 0, 1, 2, 3/
      data (ipa43( 334,j),j=1,7)/  3, 0, 1, 0, 2, 1, 3/
      data (ipa43( 335,j),j=1,7)/  3, 0, 1, 1, 0, 2, 3/
      data (ipa43( 336,j),j=1,7)/  3, 0, 1, 1, 1, 1, 3/
      data (ipa43( 337,j),j=1,7)/  3, 0, 1, 1, 2, 0, 3/
      data (ipa43( 338,j),j=1,7)/  3, 0, 1, 2, 0, 1, 3/
      data (ipa43( 339,j),j=1,7)/  3, 0, 1, 2, 1, 0, 3/
      data (ipa43( 340,j),j=1,7)/  3, 0, 1, 3, 0, 0, 3/
      data (ipa43( 341,j),j=1,7)/  3, 0, 2, 1, 1, 0, 3/
      data (ipa43( 342,j),j=1,7)/  3, 0, 2, 2, 0, 0, 3/
      data (ipa43( 343,j),j=1,7)/  3, 1, 1, 1, 1, 0, 1/
      data (ipa43( 344,j),j=1,7)/  4, 0, 0, 0, 1, 2, 3/
      data (ipa43( 345,j),j=1,7)/  4, 0, 0, 0, 2, 1, 3/
      data (ipa43( 346,j),j=1,7)/  4, 0, 0, 1, 1, 1, 2/
      data (ipa43( 347,j),j=1,7)/  4, 0, 0, 1, 2, 0, 3/
      data (ipa43( 348,j),j=1,7)/  4, 0, 1, 0, 1, 1, 3/
      data (ipa43( 349,j),j=1,7)/  4, 0, 1, 1, 0, 1, 3/
      data (ipa43( 350,j),j=1,7)/  4, 0, 1, 1, 1, 0, 3/
      data (ipa43( 351,j),j=1,7)/  4, 0, 1, 2, 0, 0, 3/
      data (ipa43( 352,j),j=1,7)/  5, 0, 0, 0, 1, 1, 3/
      data (ipa43( 353,j),j=1,7)/  5, 0, 0, 1, 1, 0, 2/
      data (ipa43( 354,j),j=1,7)/  5, 0, 1, 1, 0, 0, 3/
*********** 7 180*******************************
      data (ipa43( 355,j),j=1,7)/  0, 0, 1, 0, 1, 6, 3/
      data (ipa43( 356,j),j=1,7)/  0, 0, 1, 0, 2, 5, 3/
      data (ipa43( 357,j),j=1,7)/  0, 0, 1, 0, 3, 4, 3/
      data (ipa43( 358,j),j=1,7)/  0, 0, 1, 0, 4, 3, 3/
      data (ipa43( 359,j),j=1,7)/  0, 0, 1, 0, 5, 2, 3/
      data (ipa43( 360,j),j=1,7)/  0, 0, 1, 0, 6, 1, 3/
      data (ipa43( 361,j),j=1,7)/  0, 0, 1, 1, 0, 6, 3/
      data (ipa43( 362,j),j=1,7)/  0, 0, 1, 1, 1, 5, 3/
      data (ipa43( 363,j),j=1,7)/  0, 0, 1, 1, 2, 4, 3/
      data (ipa43( 364,j),j=1,7)/  0, 0, 1, 1, 3, 3, 3/
      data (ipa43( 365,j),j=1,7)/  0, 0, 1, 1, 4, 2, 3/
      data (ipa43( 366,j),j=1,7)/  0, 0, 1, 1, 5, 1, 3/
      data (ipa43( 367,j),j=1,7)/  0, 0, 1, 1, 6, 0, 3/
      data (ipa43( 368,j),j=1,7)/  0, 0, 1, 2, 0, 5, 3/
      data (ipa43( 369,j),j=1,7)/  0, 0, 1, 2, 1, 4, 3/
      data (ipa43( 370,j),j=1,7)/  0, 0, 1, 2, 2, 3, 3/
      data (ipa43( 371,j),j=1,7)/  0, 0, 1, 2, 3, 2, 3/
      data (ipa43( 372,j),j=1,7)/  0, 0, 1, 2, 4, 1, 3/
      data (ipa43( 373,j),j=1,7)/  0, 0, 1, 2, 5, 0, 3/
      data (ipa43( 374,j),j=1,7)/  0, 0, 1, 3, 0, 4, 3/
      data (ipa43( 375,j),j=1,7)/  0, 0, 1, 3, 1, 3, 3/
      data (ipa43( 376,j),j=1,7)/  0, 0, 1, 3, 2, 2, 3/
      data (ipa43( 377,j),j=1,7)/  0, 0, 1, 3, 3, 1, 3/
      data (ipa43( 378,j),j=1,7)/  0, 0, 1, 3, 4, 0, 3/
      data (ipa43( 379,j),j=1,7)/  0, 0, 1, 4, 0, 3, 3/
      data (ipa43( 380,j),j=1,7)/  0, 0, 1, 4, 1, 2, 3/
      data (ipa43( 381,j),j=1,7)/  0, 0, 1, 4, 2, 1, 3/
      data (ipa43( 382,j),j=1,7)/  0, 0, 1, 4, 3, 0, 3/
      data (ipa43( 383,j),j=1,7)/  0, 0, 1, 5, 0, 2, 3/
      data (ipa43( 384,j),j=1,7)/  0, 0, 1, 5, 1, 1, 3/
      data (ipa43( 385,j),j=1,7)/  0, 0, 1, 5, 2, 0, 3/
      data (ipa43( 386,j),j=1,7)/  0, 0, 1, 6, 0, 1, 3/
      data (ipa43( 387,j),j=1,7)/  0, 0, 1, 6, 1, 0, 3/
      data (ipa43( 388,j),j=1,7)/  0, 0, 2, 0, 2, 4, 3/
      data (ipa43( 389,j),j=1,7)/  0, 0, 2, 0, 3, 3, 3/
      data (ipa43( 390,j),j=1,7)/  0, 0, 2, 0, 4, 2, 3/
      data (ipa43( 391,j),j=1,7)/  0, 0, 2, 0, 5, 1, 3/
      data (ipa43( 392,j),j=1,7)/  0, 0, 2, 1, 1, 4, 3/
      data (ipa43( 393,j),j=1,7)/  0, 0, 2, 1, 2, 3, 3/
      data (ipa43( 394,j),j=1,7)/  0, 0, 2, 1, 3, 2, 3/
      data (ipa43( 395,j),j=1,7)/  0, 0, 2, 1, 4, 1, 3/
      data (ipa43( 396,j),j=1,7)/  0, 0, 2, 1, 5, 0, 3/
      data (ipa43( 397,j),j=1,7)/  0, 0, 2, 2, 0, 4, 3/
      data (ipa43( 398,j),j=1,7)/  0, 0, 2, 2, 1, 3, 3/
      data (ipa43( 399,j),j=1,7)/  0, 0, 2, 2, 2, 2, 3/
      data (ipa43( 400,j),j=1,7)/  0, 0, 2, 2, 3, 1, 3/
      data (ipa43( 401,j),j=1,7)/  0, 0, 2, 2, 4, 0, 3/
      data (ipa43( 402,j),j=1,7)/  0, 0, 2, 3, 0, 3, 3/
      data (ipa43( 403,j),j=1,7)/  0, 0, 2, 3, 1, 2, 3/
      data (ipa43( 404,j),j=1,7)/  0, 0, 2, 3, 2, 1, 3/
      data (ipa43( 405,j),j=1,7)/  0, 0, 2, 3, 3, 0, 3/
      data (ipa43( 406,j),j=1,7)/  0, 0, 2, 4, 0, 2, 3/
      data (ipa43( 407,j),j=1,7)/  0, 0, 2, 4, 1, 1, 3/
      data (ipa43( 408,j),j=1,7)/  0, 0, 2, 4, 2, 0, 3/
      data (ipa43( 409,j),j=1,7)/  0, 0, 2, 5, 0, 1, 3/
      data (ipa43( 410,j),j=1,7)/  0, 0, 2, 5, 1, 0, 3/
      data (ipa43( 411,j),j=1,7)/  0, 0, 3, 0, 3, 2, 3/
      data (ipa43( 412,j),j=1,7)/  0, 0, 3, 0, 4, 1, 3/
      data (ipa43( 413,j),j=1,7)/  0, 0, 3, 1, 1, 3, 3/
      data (ipa43( 414,j),j=1,7)/  0, 0, 3, 1, 2, 2, 3/
      data (ipa43( 415,j),j=1,7)/  0, 0, 3, 1, 3, 1, 3/
      data (ipa43( 416,j),j=1,7)/  0, 0, 3, 1, 4, 0, 3/
      data (ipa43( 417,j),j=1,7)/  0, 0, 3, 2, 1, 2, 3/
      data (ipa43( 418,j),j=1,7)/  0, 0, 3, 2, 2, 1, 3/
      data (ipa43( 419,j),j=1,7)/  0, 0, 3, 2, 3, 0, 3/
      data (ipa43( 420,j),j=1,7)/  0, 0, 3, 3, 0, 2, 3/
      data (ipa43( 421,j),j=1,7)/  0, 0, 3, 3, 1, 1, 3/
      data (ipa43( 422,j),j=1,7)/  0, 0, 3, 3, 2, 0, 3/
      data (ipa43( 423,j),j=1,7)/  0, 0, 3, 4, 0, 1, 3/
      data (ipa43( 424,j),j=1,7)/  0, 0, 3, 4, 1, 0, 3/
      data (ipa43( 425,j),j=1,7)/  0, 0, 4, 1, 1, 2, 3/
      data (ipa43( 426,j),j=1,7)/  0, 0, 4, 1, 2, 1, 3/
      data (ipa43( 427,j),j=1,7)/  0, 0, 4, 1, 3, 0, 3/
      data (ipa43( 428,j),j=1,7)/  0, 0, 4, 2, 1, 1, 3/
      data (ipa43( 429,j),j=1,7)/  0, 0, 4, 2, 2, 0, 3/
      data (ipa43( 430,j),j=1,7)/  0, 0, 4, 3, 1, 0, 3/
      data (ipa43( 431,j),j=1,7)/  0, 0, 5, 1, 1, 1, 3/
      data (ipa43( 432,j),j=1,7)/  0, 0, 5, 1, 2, 0, 3/
      data (ipa43( 433,j),j=1,7)/  0, 0, 5, 2, 1, 0, 3/
      data (ipa43( 434,j),j=1,7)/  0, 0, 6, 1, 1, 0, 3/
      data (ipa43( 435,j),j=1,7)/  0, 1, 1, 1, 1, 4, 1/
      data (ipa43( 436,j),j=1,7)/  0, 1, 1, 1, 2, 3, 3/
      data (ipa43( 437,j),j=1,7)/  0, 1, 1, 1, 3, 2, 3/
      data (ipa43( 438,j),j=1,7)/  0, 1, 1, 1, 4, 1, 3/
      data (ipa43( 439,j),j=1,7)/  0, 1, 1, 1, 5, 0, 3/
      data (ipa43( 440,j),j=1,7)/  0, 1, 1, 2, 2, 2, 2/
      data (ipa43( 441,j),j=1,7)/  0, 1, 1, 2, 3, 1, 3/
      data (ipa43( 442,j),j=1,7)/  0, 1, 1, 2, 4, 0, 3/
      data (ipa43( 443,j),j=1,7)/  0, 1, 1, 3, 3, 0, 2/
      data (ipa43( 444,j),j=1,7)/  0, 1, 2, 1, 2, 2, 3/
      data (ipa43( 445,j),j=1,7)/  0, 1, 2, 1, 3, 1, 3/
      data (ipa43( 446,j),j=1,7)/  0, 1, 2, 1, 4, 0, 3/
      data (ipa43( 447,j),j=1,7)/  0, 1, 2, 2, 1, 2, 3/
      data (ipa43( 448,j),j=1,7)/  0, 1, 2, 2, 2, 1, 3/
      data (ipa43( 449,j),j=1,7)/  0, 1, 2, 2, 3, 0, 3/
      data (ipa43( 450,j),j=1,7)/  0, 1, 2, 3, 1, 1, 3/
      data (ipa43( 451,j),j=1,7)/  0, 1, 2, 3, 2, 0, 3/
      data (ipa43( 452,j),j=1,7)/  0, 1, 2, 4, 1, 0, 3/
      data (ipa43( 453,j),j=1,7)/  0, 1, 3, 1, 3, 0, 3/
      data (ipa43( 454,j),j=1,7)/  0, 1, 3, 2, 2, 0, 3/
      data (ipa43( 455,j),j=1,7)/  0, 1, 3, 3, 1, 0, 3/
      data (ipa43( 456,j),j=1,7)/  0, 2, 2, 2, 2, 0, 1/
      data (ipa43( 457,j),j=1,7)/  1, 0, 0, 0, 1, 6, 3/
      data (ipa43( 458,j),j=1,7)/  1, 0, 0, 0, 2, 5, 3/
      data (ipa43( 459,j),j=1,7)/  1, 0, 0, 0, 3, 4, 3/
      data (ipa43( 460,j),j=1,7)/  1, 0, 0, 0, 4, 3, 3/
      data (ipa43( 461,j),j=1,7)/  1, 0, 0, 0, 5, 2, 3/
      data (ipa43( 462,j),j=1,7)/  1, 0, 0, 0, 6, 1, 3/
      data (ipa43( 463,j),j=1,7)/  1, 0, 0, 1, 1, 5, 2/
      data (ipa43( 464,j),j=1,7)/  1, 0, 0, 1, 2, 4, 3/
      data (ipa43( 465,j),j=1,7)/  1, 0, 0, 1, 3, 3, 3/
      data (ipa43( 466,j),j=1,7)/  1, 0, 0, 1, 4, 2, 3/
      data (ipa43( 467,j),j=1,7)/  1, 0, 0, 1, 5, 1, 3/
      data (ipa43( 468,j),j=1,7)/  1, 0, 0, 1, 6, 0, 3/
      data (ipa43( 469,j),j=1,7)/  1, 0, 0, 2, 2, 3, 2/
      data (ipa43( 470,j),j=1,7)/  1, 0, 0, 2, 3, 2, 3/
      data (ipa43( 471,j),j=1,7)/  1, 0, 0, 2, 4, 1, 3/
      data (ipa43( 472,j),j=1,7)/  1, 0, 0, 2, 5, 0, 3/
      data (ipa43( 473,j),j=1,7)/  1, 0, 0, 3, 3, 1, 2/
      data (ipa43( 474,j),j=1,7)/  1, 0, 0, 3, 4, 0, 3/
      data (ipa43( 475,j),j=1,7)/  1, 0, 1, 0, 1, 5, 3/
      data (ipa43( 476,j),j=1,7)/  1, 0, 1, 0, 2, 4, 3/
      data (ipa43( 477,j),j=1,7)/  1, 0, 1, 0, 3, 3, 3/
      data (ipa43( 478,j),j=1,7)/  1, 0, 1, 0, 4, 2, 3/
      data (ipa43( 479,j),j=1,7)/  1, 0, 1, 0, 5, 1, 3/
      data (ipa43( 480,j),j=1,7)/  1, 0, 1, 1, 0, 5, 3/
      data (ipa43( 481,j),j=1,7)/  1, 0, 1, 1, 1, 4, 3/
      data (ipa43( 482,j),j=1,7)/  1, 0, 1, 1, 2, 3, 3/
      data (ipa43( 483,j),j=1,7)/  1, 0, 1, 1, 3, 2, 3/
      data (ipa43( 484,j),j=1,7)/  1, 0, 1, 1, 4, 1, 3/
      data (ipa43( 485,j),j=1,7)/  1, 0, 1, 1, 5, 0, 3/
      data (ipa43( 486,j),j=1,7)/  1, 0, 1, 2, 0, 4, 3/
      data (ipa43( 487,j),j=1,7)/  1, 0, 1, 2, 1, 3, 3/
      data (ipa43( 488,j),j=1,7)/  1, 0, 1, 2, 2, 2, 3/
      data (ipa43( 489,j),j=1,7)/  1, 0, 1, 2, 3, 1, 3/
      data (ipa43( 490,j),j=1,7)/  1, 0, 1, 2, 4, 0, 3/
      data (ipa43( 491,j),j=1,7)/  1, 0, 1, 3, 0, 3, 3/
      data (ipa43( 492,j),j=1,7)/  1, 0, 1, 3, 1, 2, 3/
      data (ipa43( 493,j),j=1,7)/  1, 0, 1, 3, 2, 1, 3/
      data (ipa43( 494,j),j=1,7)/  1, 0, 1, 3, 3, 0, 3/
      data (ipa43( 495,j),j=1,7)/  1, 0, 1, 4, 0, 2, 3/
      data (ipa43( 496,j),j=1,7)/  1, 0, 1, 4, 1, 1, 3/
      data (ipa43( 497,j),j=1,7)/  1, 0, 1, 4, 2, 0, 3/
      data (ipa43( 498,j),j=1,7)/  1, 0, 1, 5, 0, 1, 3/
      data (ipa43( 499,j),j=1,7)/  1, 0, 1, 5, 1, 0, 3/
      data (ipa43( 500,j),j=1,7)/  1, 0, 1, 6, 0, 0, 3/
      data (ipa43( 501,j),j=1,7)/  1, 0, 2, 0, 2, 3, 3/
      data (ipa43( 502,j),j=1,7)/  1, 0, 2, 0, 3, 2, 3/
      data (ipa43( 503,j),j=1,7)/  1, 0, 2, 0, 4, 1, 3/
      data (ipa43( 504,j),j=1,7)/  1, 0, 2, 1, 1, 3, 3/
      data (ipa43( 505,j),j=1,7)/  1, 0, 2, 1, 2, 2, 3/
      data (ipa43( 506,j),j=1,7)/  1, 0, 2, 1, 3, 1, 3/
      data (ipa43( 507,j),j=1,7)/  1, 0, 2, 1, 4, 0, 3/
      data (ipa43( 508,j),j=1,7)/  1, 0, 2, 2, 0, 3, 3/
      data (ipa43( 509,j),j=1,7)/  1, 0, 2, 2, 1, 2, 3/
      data (ipa43( 510,j),j=1,7)/  1, 0, 2, 2, 2, 1, 3/
      data (ipa43( 511,j),j=1,7)/  1, 0, 2, 2, 3, 0, 3/
      data (ipa43( 512,j),j=1,7)/  1, 0, 2, 3, 0, 2, 3/
      data (ipa43( 513,j),j=1,7)/  1, 0, 2, 3, 1, 1, 3/
      data (ipa43( 514,j),j=1,7)/  1, 0, 2, 3, 2, 0, 3/
      data (ipa43( 515,j),j=1,7)/  1, 0, 2, 4, 0, 1, 3/
      data (ipa43( 516,j),j=1,7)/  1, 0, 2, 4, 1, 0, 3/
      data (ipa43( 517,j),j=1,7)/  1, 0, 2, 5, 0, 0, 3/
      data (ipa43( 518,j),j=1,7)/  1, 0, 3, 0, 3, 1, 3/
      data (ipa43( 519,j),j=1,7)/  1, 0, 3, 1, 1, 2, 3/
      data (ipa43( 520,j),j=1,7)/  1, 0, 3, 1, 2, 1, 3/
      data (ipa43( 521,j),j=1,7)/  1, 0, 3, 1, 3, 0, 3/
      data (ipa43( 522,j),j=1,7)/  1, 0, 3, 2, 1, 1, 3/
      data (ipa43( 523,j),j=1,7)/  1, 0, 3, 2, 2, 0, 3/
      data (ipa43( 524,j),j=1,7)/  1, 0, 3, 3, 0, 1, 3/
      data (ipa43( 525,j),j=1,7)/  1, 0, 3, 3, 1, 0, 3/
      data (ipa43( 526,j),j=1,7)/  1, 0, 3, 4, 0, 0, 3/
      data (ipa43( 527,j),j=1,7)/  1, 0, 4, 1, 1, 1, 3/
      data (ipa43( 528,j),j=1,7)/  1, 0, 4, 1, 2, 0, 3/
      data (ipa43( 529,j),j=1,7)/  1, 0, 4, 2, 1, 0, 3/
      data (ipa43( 530,j),j=1,7)/  1, 0, 5, 1, 1, 0, 3/
      data (ipa43( 531,j),j=1,7)/  1, 1, 1, 1, 1, 3, 1/
      data (ipa43( 532,j),j=1,7)/  1, 1, 1, 1, 2, 2, 3/
      data (ipa43( 533,j),j=1,7)/  1, 1, 1, 1, 3, 1, 3/
      data (ipa43( 534,j),j=1,7)/  1, 1, 1, 1, 4, 0, 3/
      data (ipa43( 535,j),j=1,7)/  1, 1, 1, 2, 2, 1, 2/
      data (ipa43( 536,j),j=1,7)/  1, 1, 1, 2, 3, 0, 3/
      data (ipa43( 537,j),j=1,7)/  1, 1, 2, 1, 2, 1, 3/
      data (ipa43( 538,j),j=1,7)/  1, 1, 2, 1, 3, 0, 3/
      data (ipa43( 539,j),j=1,7)/  1, 1, 2, 2, 1, 1, 3/
      data (ipa43( 540,j),j=1,7)/  1, 1, 2, 2, 2, 0, 3/
      data (ipa43( 541,j),j=1,7)/  1, 1, 2, 3, 1, 0, 3/
      data (ipa43( 542,j),j=1,7)/  2, 0, 0, 0, 1, 5, 3/
      data (ipa43( 543,j),j=1,7)/  2, 0, 0, 0, 2, 4, 3/
      data (ipa43( 544,j),j=1,7)/  2, 0, 0, 0, 3, 3, 3/
      data (ipa43( 545,j),j=1,7)/  2, 0, 0, 0, 4, 2, 3/
      data (ipa43( 546,j),j=1,7)/  2, 0, 0, 0, 5, 1, 3/
      data (ipa43( 547,j),j=1,7)/  2, 0, 0, 1, 1, 4, 2/
      data (ipa43( 548,j),j=1,7)/  2, 0, 0, 1, 2, 3, 3/
      data (ipa43( 549,j),j=1,7)/  2, 0, 0, 1, 3, 2, 3/
      data (ipa43( 550,j),j=1,7)/  2, 0, 0, 1, 4, 1, 3/
      data (ipa43( 551,j),j=1,7)/  2, 0, 0, 1, 5, 0, 3/
      data (ipa43( 552,j),j=1,7)/  2, 0, 0, 2, 2, 2, 2/
      data (ipa43( 553,j),j=1,7)/  2, 0, 0, 2, 3, 1, 3/
      data (ipa43( 554,j),j=1,7)/  2, 0, 0, 2, 4, 0, 3/
      data (ipa43( 555,j),j=1,7)/  2, 0, 0, 3, 3, 0, 2/
      data (ipa43( 556,j),j=1,7)/  2, 0, 1, 0, 1, 4, 3/
      data (ipa43( 557,j),j=1,7)/  2, 0, 1, 0, 2, 3, 3/
      data (ipa43( 558,j),j=1,7)/  2, 0, 1, 0, 3, 2, 3/
      data (ipa43( 559,j),j=1,7)/  2, 0, 1, 0, 4, 1, 3/
      data (ipa43( 560,j),j=1,7)/  2, 0, 1, 1, 0, 4, 3/
      data (ipa43( 561,j),j=1,7)/  2, 0, 1, 1, 1, 3, 3/
      data (ipa43( 562,j),j=1,7)/  2, 0, 1, 1, 2, 2, 3/
      data (ipa43( 563,j),j=1,7)/  2, 0, 1, 1, 3, 1, 3/
      data (ipa43( 564,j),j=1,7)/  2, 0, 1, 1, 4, 0, 3/
      data (ipa43( 565,j),j=1,7)/  2, 0, 1, 2, 0, 3, 3/
      data (ipa43( 566,j),j=1,7)/  2, 0, 1, 2, 1, 2, 3/
      data (ipa43( 567,j),j=1,7)/  2, 0, 1, 2, 2, 1, 3/
      data (ipa43( 568,j),j=1,7)/  2, 0, 1, 2, 3, 0, 3/
      data (ipa43( 569,j),j=1,7)/  2, 0, 1, 3, 0, 2, 3/
      data (ipa43( 570,j),j=1,7)/  2, 0, 1, 3, 1, 1, 3/
      data (ipa43( 571,j),j=1,7)/  2, 0, 1, 3, 2, 0, 3/
      data (ipa43( 572,j),j=1,7)/  2, 0, 1, 4, 0, 1, 3/
      data (ipa43( 573,j),j=1,7)/  2, 0, 1, 4, 1, 0, 3/
      data (ipa43( 574,j),j=1,7)/  2, 0, 1, 5, 0, 0, 3/
      data (ipa43( 575,j),j=1,7)/  2, 0, 2, 0, 2, 2, 3/
      data (ipa43( 576,j),j=1,7)/  2, 0, 2, 0, 3, 1, 3/
      data (ipa43( 577,j),j=1,7)/  2, 0, 2, 1, 1, 2, 3/
      data (ipa43( 578,j),j=1,7)/  2, 0, 2, 1, 2, 1, 3/
      data (ipa43( 579,j),j=1,7)/  2, 0, 2, 1, 3, 0, 3/
      data (ipa43( 580,j),j=1,7)/  2, 0, 2, 2, 0, 2, 3/
      data (ipa43( 581,j),j=1,7)/  2, 0, 2, 2, 1, 1, 3/
      data (ipa43( 582,j),j=1,7)/  2, 0, 2, 2, 2, 0, 3/
      data (ipa43( 583,j),j=1,7)/  2, 0, 2, 3, 0, 1, 3/
      data (ipa43( 584,j),j=1,7)/  2, 0, 2, 3, 1, 0, 3/
      data (ipa43( 585,j),j=1,7)/  2, 0, 2, 4, 0, 0, 3/
      data (ipa43( 586,j),j=1,7)/  2, 0, 3, 1, 1, 1, 3/
      data (ipa43( 587,j),j=1,7)/  2, 0, 3, 1, 2, 0, 3/
      data (ipa43( 588,j),j=1,7)/  2, 0, 3, 2, 1, 0, 3/
      data (ipa43( 589,j),j=1,7)/  2, 0, 3, 3, 0, 0, 3/
      data (ipa43( 590,j),j=1,7)/  2, 0, 4, 1, 1, 0, 3/
      data (ipa43( 591,j),j=1,7)/  2, 1, 1, 1, 1, 2, 1/
      data (ipa43( 592,j),j=1,7)/  2, 1, 1, 1, 2, 1, 3/
      data (ipa43( 593,j),j=1,7)/  2, 1, 1, 1, 3, 0, 3/
      data (ipa43( 594,j),j=1,7)/  2, 1, 1, 2, 2, 0, 2/
      data (ipa43( 595,j),j=1,7)/  2, 1, 2, 1, 2, 0, 3/
      data (ipa43( 596,j),j=1,7)/  2, 1, 2, 2, 1, 0, 3/
      data (ipa43( 597,j),j=1,7)/  3, 0, 0, 0, 1, 4, 3/
      data (ipa43( 598,j),j=1,7)/  3, 0, 0, 0, 2, 3, 3/
      data (ipa43( 599,j),j=1,7)/  3, 0, 0, 0, 3, 2, 3/
      data (ipa43( 600,j),j=1,7)/  3, 0, 0, 0, 4, 1, 3/
      data (ipa43( 601,j),j=1,7)/  3, 0, 0, 1, 1, 3, 2/
      data (ipa43( 602,j),j=1,7)/  3, 0, 0, 1, 2, 2, 3/
      data (ipa43( 603,j),j=1,7)/  3, 0, 0, 1, 3, 1, 3/
      data (ipa43( 604,j),j=1,7)/  3, 0, 0, 1, 4, 0, 3/
      data (ipa43( 605,j),j=1,7)/  3, 0, 0, 2, 2, 1, 2/
      data (ipa43( 606,j),j=1,7)/  3, 0, 0, 2, 3, 0, 3/
      data (ipa43( 607,j),j=1,7)/  3, 0, 1, 0, 1, 3, 3/
      data (ipa43( 608,j),j=1,7)/  3, 0, 1, 0, 2, 2, 3/
      data (ipa43( 609,j),j=1,7)/  3, 0, 1, 0, 3, 1, 3/
      data (ipa43( 610,j),j=1,7)/  3, 0, 1, 1, 0, 3, 3/
      data (ipa43( 611,j),j=1,7)/  3, 0, 1, 1, 1, 2, 3/
      data (ipa43( 612,j),j=1,7)/  3, 0, 1, 1, 2, 1, 3/
      data (ipa43( 613,j),j=1,7)/  3, 0, 1, 1, 3, 0, 3/
      data (ipa43( 614,j),j=1,7)/  3, 0, 1, 2, 0, 2, 3/
      data (ipa43( 615,j),j=1,7)/  3, 0, 1, 2, 1, 1, 3/
      data (ipa43( 616,j),j=1,7)/  3, 0, 1, 2, 2, 0, 3/
      data (ipa43( 617,j),j=1,7)/  3, 0, 1, 3, 0, 1, 3/
      data (ipa43( 618,j),j=1,7)/  3, 0, 1, 3, 1, 0, 3/
      data (ipa43( 619,j),j=1,7)/  3, 0, 1, 4, 0, 0, 3/
      data (ipa43( 620,j),j=1,7)/  3, 0, 2, 0, 2, 1, 3/
      data (ipa43( 621,j),j=1,7)/  3, 0, 2, 1, 1, 1, 3/
      data (ipa43( 622,j),j=1,7)/  3, 0, 2, 1, 2, 0, 3/
      data (ipa43( 623,j),j=1,7)/  3, 0, 2, 2, 0, 1, 3/
      data (ipa43( 624,j),j=1,7)/  3, 0, 2, 2, 1, 0, 3/
      data (ipa43( 625,j),j=1,7)/  3, 0, 2, 3, 0, 0, 3/
      data (ipa43( 626,j),j=1,7)/  3, 0, 3, 1, 1, 0, 3/
      data (ipa43( 627,j),j=1,7)/  3, 1, 1, 1, 1, 1, 1/
      data (ipa43( 628,j),j=1,7)/  3, 1, 1, 1, 2, 0, 3/
      data (ipa43( 629,j),j=1,7)/  4, 0, 0, 0, 1, 3, 3/
      data (ipa43( 630,j),j=1,7)/  4, 0, 0, 0, 2, 2, 3/
      data (ipa43( 631,j),j=1,7)/  4, 0, 0, 0, 3, 1, 3/
      data (ipa43( 632,j),j=1,7)/  4, 0, 0, 1, 1, 2, 2/
      data (ipa43( 633,j),j=1,7)/  4, 0, 0, 1, 2, 1, 3/
      data (ipa43( 634,j),j=1,7)/  4, 0, 0, 1, 3, 0, 3/
      data (ipa43( 635,j),j=1,7)/  4, 0, 0, 2, 2, 0, 2/
      data (ipa43( 636,j),j=1,7)/  4, 0, 1, 0, 1, 2, 3/
      data (ipa43( 637,j),j=1,7)/  4, 0, 1, 0, 2, 1, 3/
      data (ipa43( 638,j),j=1,7)/  4, 0, 1, 1, 0, 2, 3/
      data (ipa43( 639,j),j=1,7)/  4, 0, 1, 1, 1, 1, 3/
      data (ipa43( 640,j),j=1,7)/  4, 0, 1, 1, 2, 0, 3/
      data (ipa43( 641,j),j=1,7)/  4, 0, 1, 2, 0, 1, 3/
      data (ipa43( 642,j),j=1,7)/  4, 0, 1, 2, 1, 0, 3/
      data (ipa43( 643,j),j=1,7)/  4, 0, 1, 3, 0, 0, 3/
      data (ipa43( 644,j),j=1,7)/  4, 0, 2, 1, 1, 0, 3/
      data (ipa43( 645,j),j=1,7)/  4, 0, 2, 2, 0, 0, 3/
      data (ipa43( 646,j),j=1,7)/  4, 1, 1, 1, 1, 0, 1/
      data (ipa43( 647,j),j=1,7)/  5, 0, 0, 0, 1, 2, 3/
      data (ipa43( 648,j),j=1,7)/  5, 0, 0, 0, 2, 1, 3/
      data (ipa43( 649,j),j=1,7)/  5, 0, 0, 1, 1, 1, 2/
      data (ipa43( 650,j),j=1,7)/  5, 0, 0, 1, 2, 0, 3/
      data (ipa43( 651,j),j=1,7)/  5, 0, 1, 0, 1, 1, 3/
      data (ipa43( 652,j),j=1,7)/  5, 0, 1, 1, 0, 1, 3/
      data (ipa43( 653,j),j=1,7)/  5, 0, 1, 1, 1, 0, 3/
      data (ipa43( 654,j),j=1,7)/  5, 0, 1, 2, 0, 0, 3/
      data (ipa43( 655,j),j=1,7)/  6, 0, 0, 0, 1, 1, 3/
      data (ipa43( 656,j),j=1,7)/  6, 0, 0, 1, 1, 0, 2/
      data (ipa43( 657,j),j=1,7)/  6, 0, 1, 1, 0, 0, 3/
*********** 8 303*******************************
      data (nfas43(i),i=1,8)/0,0,6,26,75,174,354,657/
*********** 0 0*******************************
*********** 1 0*******************************
*********** 2 0*******************************
      data (ipa44(   1,j),j=1,7)/  0, 0, 1, 0, 1, 1, 2/
      data (ipa44(   2,j),j=1,7)/  0, 0, 1, 1, 0, 1, 2/
      data (ipa44(   3,j),j=1,7)/  0, 0, 1, 1, 1, 0, 2/
      data (ipa44(   4,j),j=1,7)/  0, 1, 1, 0, 1, 0, 2/
      data (ipa44(   5,j),j=1,7)/  1, 0, 0, 0, 1, 1, 2/
      data (ipa44(   6,j),j=1,7)/  1, 0, 0, 1, 1, 0, 1/
      data (ipa44(   7,j),j=1,7)/  1, 0, 1, 0, 0, 1, 2/
      data (ipa44(   8,j),j=1,7)/  1, 0, 1, 1, 0, 0, 2/
      data (ipa44(   9,j),j=1,7)/  1, 1, 1, 0, 0, 0, 1/
*********** 3 9*******************************
      data (ipa44(  10,j),j=1,7)/  0, 0, 1, 0, 1, 2, 2/
      data (ipa44(  11,j),j=1,7)/  0, 0, 1, 0, 2, 1, 2/
      data (ipa44(  12,j),j=1,7)/  0, 0, 1, 1, 0, 2, 2/
      data (ipa44(  13,j),j=1,7)/  0, 0, 1, 1, 1, 1, 2/
      data (ipa44(  14,j),j=1,7)/  0, 0, 1, 1, 2, 0, 2/
      data (ipa44(  15,j),j=1,7)/  0, 0, 1, 2, 0, 1, 2/
      data (ipa44(  16,j),j=1,7)/  0, 0, 1, 2, 1, 0, 2/
      data (ipa44(  17,j),j=1,7)/  0, 0, 2, 0, 1, 1, 2/
      data (ipa44(  18,j),j=1,7)/  0, 0, 2, 1, 0, 1, 2/
      data (ipa44(  19,j),j=1,7)/  0, 0, 2, 1, 1, 0, 2/
      data (ipa44(  20,j),j=1,7)/  0, 1, 1, 0, 1, 1, 2/
      data (ipa44(  21,j),j=1,7)/  0, 1, 1, 0, 2, 0, 2/
      data (ipa44(  22,j),j=1,7)/  0, 1, 1, 1, 1, 0, 1/
      data (ipa44(  23,j),j=1,7)/  0, 1, 2, 0, 1, 0, 2/
      data (ipa44(  24,j),j=1,7)/  0, 1, 2, 1, 0, 0, 2/
      data (ipa44(  25,j),j=1,7)/  1, 0, 0, 0, 1, 2, 2/
      data (ipa44(  26,j),j=1,7)/  1, 0, 0, 0, 2, 1, 2/
      data (ipa44(  27,j),j=1,7)/  1, 0, 0, 1, 1, 1, 1/
      data (ipa44(  28,j),j=1,7)/  1, 0, 0, 1, 2, 0, 2/
      data (ipa44(  29,j),j=1,7)/  1, 0, 1, 0, 0, 2, 2/
      data (ipa44(  30,j),j=1,7)/  1, 0, 1, 0, 1, 1, 2/
      data (ipa44(  31,j),j=1,7)/  1, 0, 1, 1, 0, 1, 2/
      data (ipa44(  32,j),j=1,7)/  1, 0, 1, 1, 1, 0, 2/
      data (ipa44(  33,j),j=1,7)/  1, 0, 1, 2, 0, 0, 2/
      data (ipa44(  34,j),j=1,7)/  1, 0, 2, 0, 0, 1, 2/
      data (ipa44(  35,j),j=1,7)/  1, 0, 2, 1, 0, 0, 2/
      data (ipa44(  36,j),j=1,7)/  1, 1, 1, 0, 0, 1, 1/
      data (ipa44(  37,j),j=1,7)/  1, 1, 1, 0, 1, 0, 2/
      data (ipa44(  38,j),j=1,7)/  1, 1, 2, 0, 0, 0, 2/
      data (ipa44(  39,j),j=1,7)/  2, 0, 0, 0, 1, 1, 2/
      data (ipa44(  40,j),j=1,7)/  2, 0, 0, 1, 1, 0, 1/
      data (ipa44(  41,j),j=1,7)/  2, 0, 1, 0, 0, 1, 2/
      data (ipa44(  42,j),j=1,7)/  2, 0, 1, 1, 0, 0, 2/
      data (ipa44(  43,j),j=1,7)/  2, 1, 1, 0, 0, 0, 1/
*********** 4 34*******************************
      data (ipa44(  44,j),j=1,7)/  0, 0, 1, 0, 1, 3, 2/
      data (ipa44(  45,j),j=1,7)/  0, 0, 1, 0, 2, 2, 2/
      data (ipa44(  46,j),j=1,7)/  0, 0, 1, 0, 3, 1, 2/
      data (ipa44(  47,j),j=1,7)/  0, 0, 1, 1, 0, 3, 2/
      data (ipa44(  48,j),j=1,7)/  0, 0, 1, 1, 1, 2, 2/
      data (ipa44(  49,j),j=1,7)/  0, 0, 1, 1, 2, 1, 2/
      data (ipa44(  50,j),j=1,7)/  0, 0, 1, 1, 3, 0, 2/
      data (ipa44(  51,j),j=1,7)/  0, 0, 1, 2, 0, 2, 2/
      data (ipa44(  52,j),j=1,7)/  0, 0, 1, 2, 1, 1, 2/
      data (ipa44(  53,j),j=1,7)/  0, 0, 1, 2, 2, 0, 2/
      data (ipa44(  54,j),j=1,7)/  0, 0, 1, 3, 0, 1, 2/
      data (ipa44(  55,j),j=1,7)/  0, 0, 1, 3, 1, 0, 2/
      data (ipa44(  56,j),j=1,7)/  0, 0, 2, 0, 1, 2, 2/
      data (ipa44(  57,j),j=1,7)/  0, 0, 2, 0, 2, 1, 2/
      data (ipa44(  58,j),j=1,7)/  0, 0, 2, 1, 0, 2, 2/
      data (ipa44(  59,j),j=1,7)/  0, 0, 2, 1, 1, 1, 2/
      data (ipa44(  60,j),j=1,7)/  0, 0, 2, 1, 2, 0, 2/
      data (ipa44(  61,j),j=1,7)/  0, 0, 2, 2, 0, 1, 2/
      data (ipa44(  62,j),j=1,7)/  0, 0, 2, 2, 1, 0, 2/
      data (ipa44(  63,j),j=1,7)/  0, 0, 3, 0, 1, 1, 2/
      data (ipa44(  64,j),j=1,7)/  0, 0, 3, 1, 0, 1, 2/
      data (ipa44(  65,j),j=1,7)/  0, 0, 3, 1, 1, 0, 2/
      data (ipa44(  66,j),j=1,7)/  0, 1, 1, 0, 1, 2, 2/
      data (ipa44(  67,j),j=1,7)/  0, 1, 1, 0, 2, 1, 2/
      data (ipa44(  68,j),j=1,7)/  0, 1, 1, 0, 3, 0, 2/
      data (ipa44(  69,j),j=1,7)/  0, 1, 1, 1, 1, 1, 1/
      data (ipa44(  70,j),j=1,7)/  0, 1, 1, 1, 2, 0, 2/
      data (ipa44(  71,j),j=1,7)/  0, 1, 2, 0, 1, 1, 2/
      data (ipa44(  72,j),j=1,7)/  0, 1, 2, 0, 2, 0, 2/
      data (ipa44(  73,j),j=1,7)/  0, 1, 2, 1, 0, 1, 2/
      data (ipa44(  74,j),j=1,7)/  0, 1, 2, 1, 1, 0, 2/
      data (ipa44(  75,j),j=1,7)/  0, 1, 2, 2, 0, 0, 2/
      data (ipa44(  76,j),j=1,7)/  0, 1, 3, 0, 1, 0, 2/
      data (ipa44(  77,j),j=1,7)/  0, 1, 3, 1, 0, 0, 2/
      data (ipa44(  78,j),j=1,7)/  0, 2, 2, 0, 1, 0, 2/
      data (ipa44(  79,j),j=1,7)/  1, 0, 0, 0, 1, 3, 2/
      data (ipa44(  80,j),j=1,7)/  1, 0, 0, 0, 2, 2, 2/
      data (ipa44(  81,j),j=1,7)/  1, 0, 0, 0, 3, 1, 2/
      data (ipa44(  82,j),j=1,7)/  1, 0, 0, 1, 1, 2, 1/
      data (ipa44(  83,j),j=1,7)/  1, 0, 0, 1, 2, 1, 2/
      data (ipa44(  84,j),j=1,7)/  1, 0, 0, 1, 3, 0, 2/
      data (ipa44(  85,j),j=1,7)/  1, 0, 0, 2, 2, 0, 1/
      data (ipa44(  86,j),j=1,7)/  1, 0, 1, 0, 0, 3, 2/
      data (ipa44(  87,j),j=1,7)/  1, 0, 1, 0, 1, 2, 2/
      data (ipa44(  88,j),j=1,7)/  1, 0, 1, 0, 2, 1, 2/
      data (ipa44(  89,j),j=1,7)/  1, 0, 1, 1, 0, 2, 2/
      data (ipa44(  90,j),j=1,7)/  1, 0, 1, 1, 1, 1, 2/
      data (ipa44(  91,j),j=1,7)/  1, 0, 1, 1, 2, 0, 2/
      data (ipa44(  92,j),j=1,7)/  1, 0, 1, 2, 0, 1, 2/
      data (ipa44(  93,j),j=1,7)/  1, 0, 1, 2, 1, 0, 2/
      data (ipa44(  94,j),j=1,7)/  1, 0, 1, 3, 0, 0, 2/
      data (ipa44(  95,j),j=1,7)/  1, 0, 2, 0, 0, 2, 2/
      data (ipa44(  96,j),j=1,7)/  1, 0, 2, 0, 1, 1, 2/
      data (ipa44(  97,j),j=1,7)/  1, 0, 2, 1, 0, 1, 2/
      data (ipa44(  98,j),j=1,7)/  1, 0, 2, 1, 1, 0, 2/
      data (ipa44(  99,j),j=1,7)/  1, 0, 2, 2, 0, 0, 2/
      data (ipa44( 100,j),j=1,7)/  1, 0, 3, 0, 0, 1, 2/
      data (ipa44( 101,j),j=1,7)/  1, 0, 3, 1, 0, 0, 2/
      data (ipa44( 102,j),j=1,7)/  1, 1, 1, 0, 0, 2, 1/
      data (ipa44( 103,j),j=1,7)/  1, 1, 1, 0, 1, 1, 2/
      data (ipa44( 104,j),j=1,7)/  1, 1, 1, 0, 2, 0, 2/
      data (ipa44( 105,j),j=1,7)/  1, 1, 1, 1, 1, 0, 1/
      data (ipa44( 106,j),j=1,7)/  1, 1, 2, 0, 0, 1, 2/
      data (ipa44( 107,j),j=1,7)/  1, 1, 2, 0, 1, 0, 2/
      data (ipa44( 108,j),j=1,7)/  1, 1, 2, 1, 0, 0, 2/
      data (ipa44( 109,j),j=1,7)/  1, 1, 3, 0, 0, 0, 2/
      data (ipa44( 110,j),j=1,7)/  1, 2, 2, 0, 0, 0, 1/
      data (ipa44( 111,j),j=1,7)/  2, 0, 0, 0, 1, 2, 2/
      data (ipa44( 112,j),j=1,7)/  2, 0, 0, 0, 2, 1, 2/
      data (ipa44( 113,j),j=1,7)/  2, 0, 0, 1, 1, 1, 1/
      data (ipa44( 114,j),j=1,7)/  2, 0, 0, 1, 2, 0, 2/
      data (ipa44( 115,j),j=1,7)/  2, 0, 1, 0, 0, 2, 2/
      data (ipa44( 116,j),j=1,7)/  2, 0, 1, 0, 1, 1, 2/
      data (ipa44( 117,j),j=1,7)/  2, 0, 1, 1, 0, 1, 2/
      data (ipa44( 118,j),j=1,7)/  2, 0, 1, 1, 1, 0, 2/
      data (ipa44( 119,j),j=1,7)/  2, 0, 1, 2, 0, 0, 2/
      data (ipa44( 120,j),j=1,7)/  2, 0, 2, 0, 0, 1, 2/
      data (ipa44( 121,j),j=1,7)/  2, 0, 2, 1, 0, 0, 2/
      data (ipa44( 122,j),j=1,7)/  2, 1, 1, 0, 0, 1, 1/
      data (ipa44( 123,j),j=1,7)/  2, 1, 1, 0, 1, 0, 2/
      data (ipa44( 124,j),j=1,7)/  2, 1, 2, 0, 0, 0, 2/
      data (ipa44( 125,j),j=1,7)/  3, 0, 0, 0, 1, 1, 2/
      data (ipa44( 126,j),j=1,7)/  3, 0, 0, 1, 1, 0, 1/
      data (ipa44( 127,j),j=1,7)/  3, 0, 1, 0, 0, 1, 2/
      data (ipa44( 128,j),j=1,7)/  3, 0, 1, 1, 0, 0, 2/
      data (ipa44( 129,j),j=1,7)/  3, 1, 1, 0, 0, 0, 1/
*********** 5 86*******************************
      data (ipa44( 130,j),j=1,7)/  0, 0, 1, 0, 1, 4, 2/
      data (ipa44( 131,j),j=1,7)/  0, 0, 1, 0, 2, 3, 2/
      data (ipa44( 132,j),j=1,7)/  0, 0, 1, 0, 3, 2, 2/
      data (ipa44( 133,j),j=1,7)/  0, 0, 1, 0, 4, 1, 2/
      data (ipa44( 134,j),j=1,7)/  0, 0, 1, 1, 0, 4, 2/
      data (ipa44( 135,j),j=1,7)/  0, 0, 1, 1, 1, 3, 2/
      data (ipa44( 136,j),j=1,7)/  0, 0, 1, 1, 2, 2, 2/
      data (ipa44( 137,j),j=1,7)/  0, 0, 1, 1, 3, 1, 2/
      data (ipa44( 138,j),j=1,7)/  0, 0, 1, 1, 4, 0, 2/
      data (ipa44( 139,j),j=1,7)/  0, 0, 1, 2, 0, 3, 2/
      data (ipa44( 140,j),j=1,7)/  0, 0, 1, 2, 1, 2, 2/
      data (ipa44( 141,j),j=1,7)/  0, 0, 1, 2, 2, 1, 2/
      data (ipa44( 142,j),j=1,7)/  0, 0, 1, 2, 3, 0, 2/
      data (ipa44( 143,j),j=1,7)/  0, 0, 1, 3, 0, 2, 2/
      data (ipa44( 144,j),j=1,7)/  0, 0, 1, 3, 1, 1, 2/
      data (ipa44( 145,j),j=1,7)/  0, 0, 1, 3, 2, 0, 2/
      data (ipa44( 146,j),j=1,7)/  0, 0, 1, 4, 0, 1, 2/
      data (ipa44( 147,j),j=1,7)/  0, 0, 1, 4, 1, 0, 2/
      data (ipa44( 148,j),j=1,7)/  0, 0, 2, 0, 1, 3, 2/
      data (ipa44( 149,j),j=1,7)/  0, 0, 2, 0, 2, 2, 2/
      data (ipa44( 150,j),j=1,7)/  0, 0, 2, 0, 3, 1, 2/
      data (ipa44( 151,j),j=1,7)/  0, 0, 2, 1, 0, 3, 2/
      data (ipa44( 152,j),j=1,7)/  0, 0, 2, 1, 1, 2, 2/
      data (ipa44( 153,j),j=1,7)/  0, 0, 2, 1, 2, 1, 2/
      data (ipa44( 154,j),j=1,7)/  0, 0, 2, 1, 3, 0, 2/
      data (ipa44( 155,j),j=1,7)/  0, 0, 2, 2, 0, 2, 2/
      data (ipa44( 156,j),j=1,7)/  0, 0, 2, 2, 1, 1, 2/
      data (ipa44( 157,j),j=1,7)/  0, 0, 2, 2, 2, 0, 2/
      data (ipa44( 158,j),j=1,7)/  0, 0, 2, 3, 0, 1, 2/
      data (ipa44( 159,j),j=1,7)/  0, 0, 2, 3, 1, 0, 2/
      data (ipa44( 160,j),j=1,7)/  0, 0, 3, 0, 1, 2, 2/
      data (ipa44( 161,j),j=1,7)/  0, 0, 3, 0, 2, 1, 2/
      data (ipa44( 162,j),j=1,7)/  0, 0, 3, 1, 0, 2, 2/
      data (ipa44( 163,j),j=1,7)/  0, 0, 3, 1, 1, 1, 2/
      data (ipa44( 164,j),j=1,7)/  0, 0, 3, 1, 2, 0, 2/
      data (ipa44( 165,j),j=1,7)/  0, 0, 3, 2, 0, 1, 2/
      data (ipa44( 166,j),j=1,7)/  0, 0, 3, 2, 1, 0, 2/
      data (ipa44( 167,j),j=1,7)/  0, 0, 4, 0, 1, 1, 2/
      data (ipa44( 168,j),j=1,7)/  0, 0, 4, 1, 0, 1, 2/
      data (ipa44( 169,j),j=1,7)/  0, 0, 4, 1, 1, 0, 2/
      data (ipa44( 170,j),j=1,7)/  0, 1, 1, 0, 1, 3, 2/
      data (ipa44( 171,j),j=1,7)/  0, 1, 1, 0, 2, 2, 2/
      data (ipa44( 172,j),j=1,7)/  0, 1, 1, 0, 3, 1, 2/
      data (ipa44( 173,j),j=1,7)/  0, 1, 1, 0, 4, 0, 2/
      data (ipa44( 174,j),j=1,7)/  0, 1, 1, 1, 1, 2, 1/
      data (ipa44( 175,j),j=1,7)/  0, 1, 1, 1, 2, 1, 2/
      data (ipa44( 176,j),j=1,7)/  0, 1, 1, 1, 3, 0, 2/
      data (ipa44( 177,j),j=1,7)/  0, 1, 1, 2, 2, 0, 1/
      data (ipa44( 178,j),j=1,7)/  0, 1, 2, 0, 1, 2, 2/
      data (ipa44( 179,j),j=1,7)/  0, 1, 2, 0, 2, 1, 2/
      data (ipa44( 180,j),j=1,7)/  0, 1, 2, 0, 3, 0, 2/
      data (ipa44( 181,j),j=1,7)/  0, 1, 2, 1, 0, 2, 2/
      data (ipa44( 182,j),j=1,7)/  0, 1, 2, 1, 1, 1, 2/
      data (ipa44( 183,j),j=1,7)/  0, 1, 2, 1, 2, 0, 2/
      data (ipa44( 184,j),j=1,7)/  0, 1, 2, 2, 0, 1, 2/
      data (ipa44( 185,j),j=1,7)/  0, 1, 2, 2, 1, 0, 2/
      data (ipa44( 186,j),j=1,7)/  0, 1, 2, 3, 0, 0, 2/
      data (ipa44( 187,j),j=1,7)/  0, 1, 3, 0, 1, 1, 2/
      data (ipa44( 188,j),j=1,7)/  0, 1, 3, 0, 2, 0, 2/
      data (ipa44( 189,j),j=1,7)/  0, 1, 3, 1, 0, 1, 2/
      data (ipa44( 190,j),j=1,7)/  0, 1, 3, 1, 1, 0, 2/
      data (ipa44( 191,j),j=1,7)/  0, 1, 3, 2, 0, 0, 2/
      data (ipa44( 192,j),j=1,7)/  0, 1, 4, 0, 1, 0, 2/
      data (ipa44( 193,j),j=1,7)/  0, 1, 4, 1, 0, 0, 2/
      data (ipa44( 194,j),j=1,7)/  0, 2, 2, 0, 1, 1, 2/
      data (ipa44( 195,j),j=1,7)/  0, 2, 2, 0, 2, 0, 2/
      data (ipa44( 196,j),j=1,7)/  0, 2, 2, 1, 1, 0, 1/
      data (ipa44( 197,j),j=1,7)/  0, 2, 3, 0, 1, 0, 2/
      data (ipa44( 198,j),j=1,7)/  0, 2, 3, 1, 0, 0, 2/
      data (ipa44( 199,j),j=1,7)/  1, 0, 0, 0, 1, 4, 2/
      data (ipa44( 200,j),j=1,7)/  1, 0, 0, 0, 2, 3, 2/
      data (ipa44( 201,j),j=1,7)/  1, 0, 0, 0, 3, 2, 2/
      data (ipa44( 202,j),j=1,7)/  1, 0, 0, 0, 4, 1, 2/
      data (ipa44( 203,j),j=1,7)/  1, 0, 0, 1, 1, 3, 1/
      data (ipa44( 204,j),j=1,7)/  1, 0, 0, 1, 2, 2, 2/
      data (ipa44( 205,j),j=1,7)/  1, 0, 0, 1, 3, 1, 2/
      data (ipa44( 206,j),j=1,7)/  1, 0, 0, 1, 4, 0, 2/
      data (ipa44( 207,j),j=1,7)/  1, 0, 0, 2, 2, 1, 1/
      data (ipa44( 208,j),j=1,7)/  1, 0, 0, 2, 3, 0, 2/
      data (ipa44( 209,j),j=1,7)/  1, 0, 1, 0, 0, 4, 2/
      data (ipa44( 210,j),j=1,7)/  1, 0, 1, 0, 1, 3, 2/
      data (ipa44( 211,j),j=1,7)/  1, 0, 1, 0, 2, 2, 2/
      data (ipa44( 212,j),j=1,7)/  1, 0, 1, 0, 3, 1, 2/
      data (ipa44( 213,j),j=1,7)/  1, 0, 1, 1, 0, 3, 2/
      data (ipa44( 214,j),j=1,7)/  1, 0, 1, 1, 1, 2, 2/
      data (ipa44( 215,j),j=1,7)/  1, 0, 1, 1, 2, 1, 2/
      data (ipa44( 216,j),j=1,7)/  1, 0, 1, 1, 3, 0, 2/
      data (ipa44( 217,j),j=1,7)/  1, 0, 1, 2, 0, 2, 2/
      data (ipa44( 218,j),j=1,7)/  1, 0, 1, 2, 1, 1, 2/
      data (ipa44( 219,j),j=1,7)/  1, 0, 1, 2, 2, 0, 2/
      data (ipa44( 220,j),j=1,7)/  1, 0, 1, 3, 0, 1, 2/
      data (ipa44( 221,j),j=1,7)/  1, 0, 1, 3, 1, 0, 2/
      data (ipa44( 222,j),j=1,7)/  1, 0, 1, 4, 0, 0, 2/
      data (ipa44( 223,j),j=1,7)/  1, 0, 2, 0, 0, 3, 2/
      data (ipa44( 224,j),j=1,7)/  1, 0, 2, 0, 1, 2, 2/
      data (ipa44( 225,j),j=1,7)/  1, 0, 2, 0, 2, 1, 2/
      data (ipa44( 226,j),j=1,7)/  1, 0, 2, 1, 0, 2, 2/
      data (ipa44( 227,j),j=1,7)/  1, 0, 2, 1, 1, 1, 2/
      data (ipa44( 228,j),j=1,7)/  1, 0, 2, 1, 2, 0, 2/
      data (ipa44( 229,j),j=1,7)/  1, 0, 2, 2, 0, 1, 2/
      data (ipa44( 230,j),j=1,7)/  1, 0, 2, 2, 1, 0, 2/
      data (ipa44( 231,j),j=1,7)/  1, 0, 2, 3, 0, 0, 2/
      data (ipa44( 232,j),j=1,7)/  1, 0, 3, 0, 0, 2, 2/
      data (ipa44( 233,j),j=1,7)/  1, 0, 3, 0, 1, 1, 2/
      data (ipa44( 234,j),j=1,7)/  1, 0, 3, 1, 0, 1, 2/
      data (ipa44( 235,j),j=1,7)/  1, 0, 3, 1, 1, 0, 2/
      data (ipa44( 236,j),j=1,7)/  1, 0, 3, 2, 0, 0, 2/
      data (ipa44( 237,j),j=1,7)/  1, 0, 4, 0, 0, 1, 2/
      data (ipa44( 238,j),j=1,7)/  1, 0, 4, 1, 0, 0, 2/
      data (ipa44( 239,j),j=1,7)/  1, 1, 1, 0, 0, 3, 1/
      data (ipa44( 240,j),j=1,7)/  1, 1, 1, 0, 1, 2, 2/
      data (ipa44( 241,j),j=1,7)/  1, 1, 1, 0, 2, 1, 2/
      data (ipa44( 242,j),j=1,7)/  1, 1, 1, 0, 3, 0, 2/
      data (ipa44( 243,j),j=1,7)/  1, 1, 1, 1, 1, 1, 1/
      data (ipa44( 244,j),j=1,7)/  1, 1, 1, 1, 2, 0, 2/
      data (ipa44( 245,j),j=1,7)/  1, 1, 2, 0, 0, 2, 2/
      data (ipa44( 246,j),j=1,7)/  1, 1, 2, 0, 1, 1, 2/
      data (ipa44( 247,j),j=1,7)/  1, 1, 2, 0, 2, 0, 2/
      data (ipa44( 248,j),j=1,7)/  1, 1, 2, 1, 0, 1, 2/
      data (ipa44( 249,j),j=1,7)/  1, 1, 2, 1, 1, 0, 2/
      data (ipa44( 250,j),j=1,7)/  1, 1, 2, 2, 0, 0, 2/
      data (ipa44( 251,j),j=1,7)/  1, 1, 3, 0, 0, 1, 2/
      data (ipa44( 252,j),j=1,7)/  1, 1, 3, 0, 1, 0, 2/
      data (ipa44( 253,j),j=1,7)/  1, 1, 3, 1, 0, 0, 2/
      data (ipa44( 254,j),j=1,7)/  1, 1, 4, 0, 0, 0, 2/
      data (ipa44( 255,j),j=1,7)/  1, 2, 2, 0, 0, 1, 1/
      data (ipa44( 256,j),j=1,7)/  1, 2, 2, 0, 1, 0, 2/
      data (ipa44( 257,j),j=1,7)/  1, 2, 3, 0, 0, 0, 2/
      data (ipa44( 258,j),j=1,7)/  2, 0, 0, 0, 1, 3, 2/
      data (ipa44( 259,j),j=1,7)/  2, 0, 0, 0, 2, 2, 2/
      data (ipa44( 260,j),j=1,7)/  2, 0, 0, 0, 3, 1, 2/
      data (ipa44( 261,j),j=1,7)/  2, 0, 0, 1, 1, 2, 1/
      data (ipa44( 262,j),j=1,7)/  2, 0, 0, 1, 2, 1, 2/
      data (ipa44( 263,j),j=1,7)/  2, 0, 0, 1, 3, 0, 2/
      data (ipa44( 264,j),j=1,7)/  2, 0, 0, 2, 2, 0, 1/
      data (ipa44( 265,j),j=1,7)/  2, 0, 1, 0, 0, 3, 2/
      data (ipa44( 266,j),j=1,7)/  2, 0, 1, 0, 1, 2, 2/
      data (ipa44( 267,j),j=1,7)/  2, 0, 1, 0, 2, 1, 2/
      data (ipa44( 268,j),j=1,7)/  2, 0, 1, 1, 0, 2, 2/
      data (ipa44( 269,j),j=1,7)/  2, 0, 1, 1, 1, 1, 2/
      data (ipa44( 270,j),j=1,7)/  2, 0, 1, 1, 2, 0, 2/
      data (ipa44( 271,j),j=1,7)/  2, 0, 1, 2, 0, 1, 2/
      data (ipa44( 272,j),j=1,7)/  2, 0, 1, 2, 1, 0, 2/
      data (ipa44( 273,j),j=1,7)/  2, 0, 1, 3, 0, 0, 2/
      data (ipa44( 274,j),j=1,7)/  2, 0, 2, 0, 0, 2, 2/
      data (ipa44( 275,j),j=1,7)/  2, 0, 2, 0, 1, 1, 2/
      data (ipa44( 276,j),j=1,7)/  2, 0, 2, 1, 0, 1, 2/
      data (ipa44( 277,j),j=1,7)/  2, 0, 2, 1, 1, 0, 2/
      data (ipa44( 278,j),j=1,7)/  2, 0, 2, 2, 0, 0, 2/
      data (ipa44( 279,j),j=1,7)/  2, 0, 3, 0, 0, 1, 2/
      data (ipa44( 280,j),j=1,7)/  2, 0, 3, 1, 0, 0, 2/
      data (ipa44( 281,j),j=1,7)/  2, 1, 1, 0, 0, 2, 1/
      data (ipa44( 282,j),j=1,7)/  2, 1, 1, 0, 1, 1, 2/
      data (ipa44( 283,j),j=1,7)/  2, 1, 1, 0, 2, 0, 2/
      data (ipa44( 284,j),j=1,7)/  2, 1, 1, 1, 1, 0, 1/
      data (ipa44( 285,j),j=1,7)/  2, 1, 2, 0, 0, 1, 2/
      data (ipa44( 286,j),j=1,7)/  2, 1, 2, 0, 1, 0, 2/
      data (ipa44( 287,j),j=1,7)/  2, 1, 2, 1, 0, 0, 2/
      data (ipa44( 288,j),j=1,7)/  2, 1, 3, 0, 0, 0, 2/
      data (ipa44( 289,j),j=1,7)/  2, 2, 2, 0, 0, 0, 1/
      data (ipa44( 290,j),j=1,7)/  3, 0, 0, 0, 1, 2, 2/
      data (ipa44( 291,j),j=1,7)/  3, 0, 0, 0, 2, 1, 2/
      data (ipa44( 292,j),j=1,7)/  3, 0, 0, 1, 1, 1, 1/
      data (ipa44( 293,j),j=1,7)/  3, 0, 0, 1, 2, 0, 2/
      data (ipa44( 294,j),j=1,7)/  3, 0, 1, 0, 0, 2, 2/
      data (ipa44( 295,j),j=1,7)/  3, 0, 1, 0, 1, 1, 2/
      data (ipa44( 296,j),j=1,7)/  3, 0, 1, 1, 0, 1, 2/
      data (ipa44( 297,j),j=1,7)/  3, 0, 1, 1, 1, 0, 2/
      data (ipa44( 298,j),j=1,7)/  3, 0, 1, 2, 0, 0, 2/
      data (ipa44( 299,j),j=1,7)/  3, 0, 2, 0, 0, 1, 2/
      data (ipa44( 300,j),j=1,7)/  3, 0, 2, 1, 0, 0, 2/
      data (ipa44( 301,j),j=1,7)/  3, 1, 1, 0, 0, 1, 1/
      data (ipa44( 302,j),j=1,7)/  3, 1, 1, 0, 1, 0, 2/
      data (ipa44( 303,j),j=1,7)/  3, 1, 2, 0, 0, 0, 2/
      data (ipa44( 304,j),j=1,7)/  4, 0, 0, 0, 1, 1, 2/
      data (ipa44( 305,j),j=1,7)/  4, 0, 0, 1, 1, 0, 1/
      data (ipa44( 306,j),j=1,7)/  4, 0, 1, 0, 0, 1, 2/
      data (ipa44( 307,j),j=1,7)/  4, 0, 1, 1, 0, 0, 2/
      data (ipa44( 308,j),j=1,7)/  4, 1, 1, 0, 0, 0, 1/
*********** 6 179*******************************
      data (ipa44( 309,j),j=1,7)/  0, 0, 1, 0, 1, 5, 2/
      data (ipa44( 310,j),j=1,7)/  0, 0, 1, 0, 2, 4, 2/
      data (ipa44( 311,j),j=1,7)/  0, 0, 1, 0, 3, 3, 2/
      data (ipa44( 312,j),j=1,7)/  0, 0, 1, 0, 4, 2, 2/
      data (ipa44( 313,j),j=1,7)/  0, 0, 1, 0, 5, 1, 2/
      data (ipa44( 314,j),j=1,7)/  0, 0, 1, 1, 0, 5, 2/
      data (ipa44( 315,j),j=1,7)/  0, 0, 1, 1, 1, 4, 2/
      data (ipa44( 316,j),j=1,7)/  0, 0, 1, 1, 2, 3, 2/
      data (ipa44( 317,j),j=1,7)/  0, 0, 1, 1, 3, 2, 2/
      data (ipa44( 318,j),j=1,7)/  0, 0, 1, 1, 4, 1, 2/
      data (ipa44( 319,j),j=1,7)/  0, 0, 1, 1, 5, 0, 2/
      data (ipa44( 320,j),j=1,7)/  0, 0, 1, 2, 0, 4, 2/
      data (ipa44( 321,j),j=1,7)/  0, 0, 1, 2, 1, 3, 2/
      data (ipa44( 322,j),j=1,7)/  0, 0, 1, 2, 2, 2, 2/
      data (ipa44( 323,j),j=1,7)/  0, 0, 1, 2, 3, 1, 2/
      data (ipa44( 324,j),j=1,7)/  0, 0, 1, 2, 4, 0, 2/
      data (ipa44( 325,j),j=1,7)/  0, 0, 1, 3, 0, 3, 2/
      data (ipa44( 326,j),j=1,7)/  0, 0, 1, 3, 1, 2, 2/
      data (ipa44( 327,j),j=1,7)/  0, 0, 1, 3, 2, 1, 2/
      data (ipa44( 328,j),j=1,7)/  0, 0, 1, 3, 3, 0, 2/
      data (ipa44( 329,j),j=1,7)/  0, 0, 1, 4, 0, 2, 2/
      data (ipa44( 330,j),j=1,7)/  0, 0, 1, 4, 1, 1, 2/
      data (ipa44( 331,j),j=1,7)/  0, 0, 1, 4, 2, 0, 2/
      data (ipa44( 332,j),j=1,7)/  0, 0, 1, 5, 0, 1, 2/
      data (ipa44( 333,j),j=1,7)/  0, 0, 1, 5, 1, 0, 2/
      data (ipa44( 334,j),j=1,7)/  0, 0, 2, 0, 1, 4, 2/
      data (ipa44( 335,j),j=1,7)/  0, 0, 2, 0, 2, 3, 2/
      data (ipa44( 336,j),j=1,7)/  0, 0, 2, 0, 3, 2, 2/
      data (ipa44( 337,j),j=1,7)/  0, 0, 2, 0, 4, 1, 2/
      data (ipa44( 338,j),j=1,7)/  0, 0, 2, 1, 0, 4, 2/
      data (ipa44( 339,j),j=1,7)/  0, 0, 2, 1, 1, 3, 2/
      data (ipa44( 340,j),j=1,7)/  0, 0, 2, 1, 2, 2, 2/
      data (ipa44( 341,j),j=1,7)/  0, 0, 2, 1, 3, 1, 2/
      data (ipa44( 342,j),j=1,7)/  0, 0, 2, 1, 4, 0, 2/
      data (ipa44( 343,j),j=1,7)/  0, 0, 2, 2, 0, 3, 2/
      data (ipa44( 344,j),j=1,7)/  0, 0, 2, 2, 1, 2, 2/
      data (ipa44( 345,j),j=1,7)/  0, 0, 2, 2, 2, 1, 2/
      data (ipa44( 346,j),j=1,7)/  0, 0, 2, 2, 3, 0, 2/
      data (ipa44( 347,j),j=1,7)/  0, 0, 2, 3, 0, 2, 2/
      data (ipa44( 348,j),j=1,7)/  0, 0, 2, 3, 1, 1, 2/
      data (ipa44( 349,j),j=1,7)/  0, 0, 2, 3, 2, 0, 2/
      data (ipa44( 350,j),j=1,7)/  0, 0, 2, 4, 0, 1, 2/
      data (ipa44( 351,j),j=1,7)/  0, 0, 2, 4, 1, 0, 2/
      data (ipa44( 352,j),j=1,7)/  0, 0, 3, 0, 1, 3, 2/
      data (ipa44( 353,j),j=1,7)/  0, 0, 3, 0, 2, 2, 2/
      data (ipa44( 354,j),j=1,7)/  0, 0, 3, 0, 3, 1, 2/
      data (ipa44( 355,j),j=1,7)/  0, 0, 3, 1, 0, 3, 2/
      data (ipa44( 356,j),j=1,7)/  0, 0, 3, 1, 1, 2, 2/
      data (ipa44( 357,j),j=1,7)/  0, 0, 3, 1, 2, 1, 2/
      data (ipa44( 358,j),j=1,7)/  0, 0, 3, 1, 3, 0, 2/
      data (ipa44( 359,j),j=1,7)/  0, 0, 3, 2, 0, 2, 2/
      data (ipa44( 360,j),j=1,7)/  0, 0, 3, 2, 1, 1, 2/
      data (ipa44( 361,j),j=1,7)/  0, 0, 3, 2, 2, 0, 2/
      data (ipa44( 362,j),j=1,7)/  0, 0, 3, 3, 0, 1, 2/
      data (ipa44( 363,j),j=1,7)/  0, 0, 3, 3, 1, 0, 2/
      data (ipa44( 364,j),j=1,7)/  0, 0, 4, 0, 1, 2, 2/
      data (ipa44( 365,j),j=1,7)/  0, 0, 4, 0, 2, 1, 2/
      data (ipa44( 366,j),j=1,7)/  0, 0, 4, 1, 0, 2, 2/
      data (ipa44( 367,j),j=1,7)/  0, 0, 4, 1, 1, 1, 2/
      data (ipa44( 368,j),j=1,7)/  0, 0, 4, 1, 2, 0, 2/
      data (ipa44( 369,j),j=1,7)/  0, 0, 4, 2, 0, 1, 2/
      data (ipa44( 370,j),j=1,7)/  0, 0, 4, 2, 1, 0, 2/
      data (ipa44( 371,j),j=1,7)/  0, 0, 5, 0, 1, 1, 2/
      data (ipa44( 372,j),j=1,7)/  0, 0, 5, 1, 0, 1, 2/
      data (ipa44( 373,j),j=1,7)/  0, 0, 5, 1, 1, 0, 2/
      data (ipa44( 374,j),j=1,7)/  0, 1, 1, 0, 1, 4, 2/
      data (ipa44( 375,j),j=1,7)/  0, 1, 1, 0, 2, 3, 2/
      data (ipa44( 376,j),j=1,7)/  0, 1, 1, 0, 3, 2, 2/
      data (ipa44( 377,j),j=1,7)/  0, 1, 1, 0, 4, 1, 2/
      data (ipa44( 378,j),j=1,7)/  0, 1, 1, 0, 5, 0, 2/
      data (ipa44( 379,j),j=1,7)/  0, 1, 1, 1, 1, 3, 1/
      data (ipa44( 380,j),j=1,7)/  0, 1, 1, 1, 2, 2, 2/
      data (ipa44( 381,j),j=1,7)/  0, 1, 1, 1, 3, 1, 2/
      data (ipa44( 382,j),j=1,7)/  0, 1, 1, 1, 4, 0, 2/
      data (ipa44( 383,j),j=1,7)/  0, 1, 1, 2, 2, 1, 1/
      data (ipa44( 384,j),j=1,7)/  0, 1, 1, 2, 3, 0, 2/
      data (ipa44( 385,j),j=1,7)/  0, 1, 2, 0, 1, 3, 2/
      data (ipa44( 386,j),j=1,7)/  0, 1, 2, 0, 2, 2, 2/
      data (ipa44( 387,j),j=1,7)/  0, 1, 2, 0, 3, 1, 2/
      data (ipa44( 388,j),j=1,7)/  0, 1, 2, 0, 4, 0, 2/
      data (ipa44( 389,j),j=1,7)/  0, 1, 2, 1, 0, 3, 2/
      data (ipa44( 390,j),j=1,7)/  0, 1, 2, 1, 1, 2, 2/
      data (ipa44( 391,j),j=1,7)/  0, 1, 2, 1, 2, 1, 2/
      data (ipa44( 392,j),j=1,7)/  0, 1, 2, 1, 3, 0, 2/
      data (ipa44( 393,j),j=1,7)/  0, 1, 2, 2, 0, 2, 2/
      data (ipa44( 394,j),j=1,7)/  0, 1, 2, 2, 1, 1, 2/
      data (ipa44( 395,j),j=1,7)/  0, 1, 2, 2, 2, 0, 2/
      data (ipa44( 396,j),j=1,7)/  0, 1, 2, 3, 0, 1, 2/
      data (ipa44( 397,j),j=1,7)/  0, 1, 2, 3, 1, 0, 2/
      data (ipa44( 398,j),j=1,7)/  0, 1, 2, 4, 0, 0, 2/
      data (ipa44( 399,j),j=1,7)/  0, 1, 3, 0, 1, 2, 2/
      data (ipa44( 400,j),j=1,7)/  0, 1, 3, 0, 2, 1, 2/
      data (ipa44( 401,j),j=1,7)/  0, 1, 3, 0, 3, 0, 2/
      data (ipa44( 402,j),j=1,7)/  0, 1, 3, 1, 0, 2, 2/
      data (ipa44( 403,j),j=1,7)/  0, 1, 3, 1, 1, 1, 2/
      data (ipa44( 404,j),j=1,7)/  0, 1, 3, 1, 2, 0, 2/
      data (ipa44( 405,j),j=1,7)/  0, 1, 3, 2, 0, 1, 2/
      data (ipa44( 406,j),j=1,7)/  0, 1, 3, 2, 1, 0, 2/
      data (ipa44( 407,j),j=1,7)/  0, 1, 3, 3, 0, 0, 2/
      data (ipa44( 408,j),j=1,7)/  0, 1, 4, 0, 1, 1, 2/
      data (ipa44( 409,j),j=1,7)/  0, 1, 4, 0, 2, 0, 2/
      data (ipa44( 410,j),j=1,7)/  0, 1, 4, 1, 0, 1, 2/
      data (ipa44( 411,j),j=1,7)/  0, 1, 4, 1, 1, 0, 2/
      data (ipa44( 412,j),j=1,7)/  0, 1, 4, 2, 0, 0, 2/
      data (ipa44( 413,j),j=1,7)/  0, 1, 5, 0, 1, 0, 2/
      data (ipa44( 414,j),j=1,7)/  0, 1, 5, 1, 0, 0, 2/
      data (ipa44( 415,j),j=1,7)/  0, 2, 2, 0, 1, 2, 2/
      data (ipa44( 416,j),j=1,7)/  0, 2, 2, 0, 2, 1, 2/
      data (ipa44( 417,j),j=1,7)/  0, 2, 2, 0, 3, 0, 2/
      data (ipa44( 418,j),j=1,7)/  0, 2, 2, 1, 1, 1, 1/
      data (ipa44( 419,j),j=1,7)/  0, 2, 2, 1, 2, 0, 2/
      data (ipa44( 420,j),j=1,7)/  0, 2, 3, 0, 1, 1, 2/
      data (ipa44( 421,j),j=1,7)/  0, 2, 3, 0, 2, 0, 2/
      data (ipa44( 422,j),j=1,7)/  0, 2, 3, 1, 0, 1, 2/
      data (ipa44( 423,j),j=1,7)/  0, 2, 3, 1, 1, 0, 2/
      data (ipa44( 424,j),j=1,7)/  0, 2, 3, 2, 0, 0, 2/
      data (ipa44( 425,j),j=1,7)/  0, 2, 4, 0, 1, 0, 2/
      data (ipa44( 426,j),j=1,7)/  0, 2, 4, 1, 0, 0, 2/
      data (ipa44( 427,j),j=1,7)/  0, 3, 3, 0, 1, 0, 2/
      data (ipa44( 428,j),j=1,7)/  1, 0, 0, 0, 1, 5, 2/
      data (ipa44( 429,j),j=1,7)/  1, 0, 0, 0, 2, 4, 2/
      data (ipa44( 430,j),j=1,7)/  1, 0, 0, 0, 3, 3, 2/
      data (ipa44( 431,j),j=1,7)/  1, 0, 0, 0, 4, 2, 2/
      data (ipa44( 432,j),j=1,7)/  1, 0, 0, 0, 5, 1, 2/
      data (ipa44( 433,j),j=1,7)/  1, 0, 0, 1, 1, 4, 1/
      data (ipa44( 434,j),j=1,7)/  1, 0, 0, 1, 2, 3, 2/
      data (ipa44( 435,j),j=1,7)/  1, 0, 0, 1, 3, 2, 2/
      data (ipa44( 436,j),j=1,7)/  1, 0, 0, 1, 4, 1, 2/
      data (ipa44( 437,j),j=1,7)/  1, 0, 0, 1, 5, 0, 2/
      data (ipa44( 438,j),j=1,7)/  1, 0, 0, 2, 2, 2, 1/
      data (ipa44( 439,j),j=1,7)/  1, 0, 0, 2, 3, 1, 2/
      data (ipa44( 440,j),j=1,7)/  1, 0, 0, 2, 4, 0, 2/
      data (ipa44( 441,j),j=1,7)/  1, 0, 0, 3, 3, 0, 1/
      data (ipa44( 442,j),j=1,7)/  1, 0, 1, 0, 0, 5, 2/
      data (ipa44( 443,j),j=1,7)/  1, 0, 1, 0, 1, 4, 2/
      data (ipa44( 444,j),j=1,7)/  1, 0, 1, 0, 2, 3, 2/
      data (ipa44( 445,j),j=1,7)/  1, 0, 1, 0, 3, 2, 2/
      data (ipa44( 446,j),j=1,7)/  1, 0, 1, 0, 4, 1, 2/
      data (ipa44( 447,j),j=1,7)/  1, 0, 1, 1, 0, 4, 2/
      data (ipa44( 448,j),j=1,7)/  1, 0, 1, 1, 1, 3, 2/
      data (ipa44( 449,j),j=1,7)/  1, 0, 1, 1, 2, 2, 2/
      data (ipa44( 450,j),j=1,7)/  1, 0, 1, 1, 3, 1, 2/
      data (ipa44( 451,j),j=1,7)/  1, 0, 1, 1, 4, 0, 2/
      data (ipa44( 452,j),j=1,7)/  1, 0, 1, 2, 0, 3, 2/
      data (ipa44( 453,j),j=1,7)/  1, 0, 1, 2, 1, 2, 2/
      data (ipa44( 454,j),j=1,7)/  1, 0, 1, 2, 2, 1, 2/
      data (ipa44( 455,j),j=1,7)/  1, 0, 1, 2, 3, 0, 2/
      data (ipa44( 456,j),j=1,7)/  1, 0, 1, 3, 0, 2, 2/
      data (ipa44( 457,j),j=1,7)/  1, 0, 1, 3, 1, 1, 2/
      data (ipa44( 458,j),j=1,7)/  1, 0, 1, 3, 2, 0, 2/
      data (ipa44( 459,j),j=1,7)/  1, 0, 1, 4, 0, 1, 2/
      data (ipa44( 460,j),j=1,7)/  1, 0, 1, 4, 1, 0, 2/
      data (ipa44( 461,j),j=1,7)/  1, 0, 1, 5, 0, 0, 2/
      data (ipa44( 462,j),j=1,7)/  1, 0, 2, 0, 0, 4, 2/
      data (ipa44( 463,j),j=1,7)/  1, 0, 2, 0, 1, 3, 2/
      data (ipa44( 464,j),j=1,7)/  1, 0, 2, 0, 2, 2, 2/
      data (ipa44( 465,j),j=1,7)/  1, 0, 2, 0, 3, 1, 2/
      data (ipa44( 466,j),j=1,7)/  1, 0, 2, 1, 0, 3, 2/
      data (ipa44( 467,j),j=1,7)/  1, 0, 2, 1, 1, 2, 2/
      data (ipa44( 468,j),j=1,7)/  1, 0, 2, 1, 2, 1, 2/
      data (ipa44( 469,j),j=1,7)/  1, 0, 2, 1, 3, 0, 2/
      data (ipa44( 470,j),j=1,7)/  1, 0, 2, 2, 0, 2, 2/
      data (ipa44( 471,j),j=1,7)/  1, 0, 2, 2, 1, 1, 2/
      data (ipa44( 472,j),j=1,7)/  1, 0, 2, 2, 2, 0, 2/
      data (ipa44( 473,j),j=1,7)/  1, 0, 2, 3, 0, 1, 2/
      data (ipa44( 474,j),j=1,7)/  1, 0, 2, 3, 1, 0, 2/
      data (ipa44( 475,j),j=1,7)/  1, 0, 2, 4, 0, 0, 2/
      data (ipa44( 476,j),j=1,7)/  1, 0, 3, 0, 0, 3, 2/
      data (ipa44( 477,j),j=1,7)/  1, 0, 3, 0, 1, 2, 2/
      data (ipa44( 478,j),j=1,7)/  1, 0, 3, 0, 2, 1, 2/
      data (ipa44( 479,j),j=1,7)/  1, 0, 3, 1, 0, 2, 2/
      data (ipa44( 480,j),j=1,7)/  1, 0, 3, 1, 1, 1, 2/
      data (ipa44( 481,j),j=1,7)/  1, 0, 3, 1, 2, 0, 2/
      data (ipa44( 482,j),j=1,7)/  1, 0, 3, 2, 0, 1, 2/
      data (ipa44( 483,j),j=1,7)/  1, 0, 3, 2, 1, 0, 2/
      data (ipa44( 484,j),j=1,7)/  1, 0, 3, 3, 0, 0, 2/
      data (ipa44( 485,j),j=1,7)/  1, 0, 4, 0, 0, 2, 2/
      data (ipa44( 486,j),j=1,7)/  1, 0, 4, 0, 1, 1, 2/
      data (ipa44( 487,j),j=1,7)/  1, 0, 4, 1, 0, 1, 2/
      data (ipa44( 488,j),j=1,7)/  1, 0, 4, 1, 1, 0, 2/
      data (ipa44( 489,j),j=1,7)/  1, 0, 4, 2, 0, 0, 2/
      data (ipa44( 490,j),j=1,7)/  1, 0, 5, 0, 0, 1, 2/
      data (ipa44( 491,j),j=1,7)/  1, 0, 5, 1, 0, 0, 2/
      data (ipa44( 492,j),j=1,7)/  1, 1, 1, 0, 0, 4, 1/
      data (ipa44( 493,j),j=1,7)/  1, 1, 1, 0, 1, 3, 2/
      data (ipa44( 494,j),j=1,7)/  1, 1, 1, 0, 2, 2, 2/
      data (ipa44( 495,j),j=1,7)/  1, 1, 1, 0, 3, 1, 2/
      data (ipa44( 496,j),j=1,7)/  1, 1, 1, 0, 4, 0, 2/
      data (ipa44( 497,j),j=1,7)/  1, 1, 1, 1, 1, 2, 1/
      data (ipa44( 498,j),j=1,7)/  1, 1, 1, 1, 2, 1, 2/
      data (ipa44( 499,j),j=1,7)/  1, 1, 1, 1, 3, 0, 2/
      data (ipa44( 500,j),j=1,7)/  1, 1, 1, 2, 2, 0, 1/
      data (ipa44( 501,j),j=1,7)/  1, 1, 2, 0, 0, 3, 2/
      data (ipa44( 502,j),j=1,7)/  1, 1, 2, 0, 1, 2, 2/
      data (ipa44( 503,j),j=1,7)/  1, 1, 2, 0, 2, 1, 2/
      data (ipa44( 504,j),j=1,7)/  1, 1, 2, 0, 3, 0, 2/
      data (ipa44( 505,j),j=1,7)/  1, 1, 2, 1, 0, 2, 2/
      data (ipa44( 506,j),j=1,7)/  1, 1, 2, 1, 1, 1, 2/
      data (ipa44( 507,j),j=1,7)/  1, 1, 2, 1, 2, 0, 2/
      data (ipa44( 508,j),j=1,7)/  1, 1, 2, 2, 0, 1, 2/
      data (ipa44( 509,j),j=1,7)/  1, 1, 2, 2, 1, 0, 2/
      data (ipa44( 510,j),j=1,7)/  1, 1, 2, 3, 0, 0, 2/
      data (ipa44( 511,j),j=1,7)/  1, 1, 3, 0, 0, 2, 2/
      data (ipa44( 512,j),j=1,7)/  1, 1, 3, 0, 1, 1, 2/
      data (ipa44( 513,j),j=1,7)/  1, 1, 3, 0, 2, 0, 2/
      data (ipa44( 514,j),j=1,7)/  1, 1, 3, 1, 0, 1, 2/
      data (ipa44( 515,j),j=1,7)/  1, 1, 3, 1, 1, 0, 2/
      data (ipa44( 516,j),j=1,7)/  1, 1, 3, 2, 0, 0, 2/
      data (ipa44( 517,j),j=1,7)/  1, 1, 4, 0, 0, 1, 2/
      data (ipa44( 518,j),j=1,7)/  1, 1, 4, 0, 1, 0, 2/
      data (ipa44( 519,j),j=1,7)/  1, 1, 4, 1, 0, 0, 2/
      data (ipa44( 520,j),j=1,7)/  1, 1, 5, 0, 0, 0, 2/
      data (ipa44( 521,j),j=1,7)/  1, 2, 2, 0, 0, 2, 1/
      data (ipa44( 522,j),j=1,7)/  1, 2, 2, 0, 1, 1, 2/
      data (ipa44( 523,j),j=1,7)/  1, 2, 2, 0, 2, 0, 2/
      data (ipa44( 524,j),j=1,7)/  1, 2, 2, 1, 1, 0, 1/
      data (ipa44( 525,j),j=1,7)/  1, 2, 3, 0, 0, 1, 2/
      data (ipa44( 526,j),j=1,7)/  1, 2, 3, 0, 1, 0, 2/
      data (ipa44( 527,j),j=1,7)/  1, 2, 3, 1, 0, 0, 2/
      data (ipa44( 528,j),j=1,7)/  1, 2, 4, 0, 0, 0, 2/
      data (ipa44( 529,j),j=1,7)/  1, 3, 3, 0, 0, 0, 1/
      data (ipa44( 530,j),j=1,7)/  2, 0, 0, 0, 1, 4, 2/
      data (ipa44( 531,j),j=1,7)/  2, 0, 0, 0, 2, 3, 2/
      data (ipa44( 532,j),j=1,7)/  2, 0, 0, 0, 3, 2, 2/
      data (ipa44( 533,j),j=1,7)/  2, 0, 0, 0, 4, 1, 2/
      data (ipa44( 534,j),j=1,7)/  2, 0, 0, 1, 1, 3, 1/
      data (ipa44( 535,j),j=1,7)/  2, 0, 0, 1, 2, 2, 2/
      data (ipa44( 536,j),j=1,7)/  2, 0, 0, 1, 3, 1, 2/
      data (ipa44( 537,j),j=1,7)/  2, 0, 0, 1, 4, 0, 2/
      data (ipa44( 538,j),j=1,7)/  2, 0, 0, 2, 2, 1, 1/
      data (ipa44( 539,j),j=1,7)/  2, 0, 0, 2, 3, 0, 2/
      data (ipa44( 540,j),j=1,7)/  2, 0, 1, 0, 0, 4, 2/
      data (ipa44( 541,j),j=1,7)/  2, 0, 1, 0, 1, 3, 2/
      data (ipa44( 542,j),j=1,7)/  2, 0, 1, 0, 2, 2, 2/
      data (ipa44( 543,j),j=1,7)/  2, 0, 1, 0, 3, 1, 2/
      data (ipa44( 544,j),j=1,7)/  2, 0, 1, 1, 0, 3, 2/
      data (ipa44( 545,j),j=1,7)/  2, 0, 1, 1, 1, 2, 2/
      data (ipa44( 546,j),j=1,7)/  2, 0, 1, 1, 2, 1, 2/
      data (ipa44( 547,j),j=1,7)/  2, 0, 1, 1, 3, 0, 2/
      data (ipa44( 548,j),j=1,7)/  2, 0, 1, 2, 0, 2, 2/
      data (ipa44( 549,j),j=1,7)/  2, 0, 1, 2, 1, 1, 2/
      data (ipa44( 550,j),j=1,7)/  2, 0, 1, 2, 2, 0, 2/
      data (ipa44( 551,j),j=1,7)/  2, 0, 1, 3, 0, 1, 2/
      data (ipa44( 552,j),j=1,7)/  2, 0, 1, 3, 1, 0, 2/
      data (ipa44( 553,j),j=1,7)/  2, 0, 1, 4, 0, 0, 2/
      data (ipa44( 554,j),j=1,7)/  2, 0, 2, 0, 0, 3, 2/
      data (ipa44( 555,j),j=1,7)/  2, 0, 2, 0, 1, 2, 2/
      data (ipa44( 556,j),j=1,7)/  2, 0, 2, 0, 2, 1, 2/
      data (ipa44( 557,j),j=1,7)/  2, 0, 2, 1, 0, 2, 2/
      data (ipa44( 558,j),j=1,7)/  2, 0, 2, 1, 1, 1, 2/
      data (ipa44( 559,j),j=1,7)/  2, 0, 2, 1, 2, 0, 2/
      data (ipa44( 560,j),j=1,7)/  2, 0, 2, 2, 0, 1, 2/
      data (ipa44( 561,j),j=1,7)/  2, 0, 2, 2, 1, 0, 2/
      data (ipa44( 562,j),j=1,7)/  2, 0, 2, 3, 0, 0, 2/
      data (ipa44( 563,j),j=1,7)/  2, 0, 3, 0, 0, 2, 2/
      data (ipa44( 564,j),j=1,7)/  2, 0, 3, 0, 1, 1, 2/
      data (ipa44( 565,j),j=1,7)/  2, 0, 3, 1, 0, 1, 2/
      data (ipa44( 566,j),j=1,7)/  2, 0, 3, 1, 1, 0, 2/
      data (ipa44( 567,j),j=1,7)/  2, 0, 3, 2, 0, 0, 2/
      data (ipa44( 568,j),j=1,7)/  2, 0, 4, 0, 0, 1, 2/
      data (ipa44( 569,j),j=1,7)/  2, 0, 4, 1, 0, 0, 2/
      data (ipa44( 570,j),j=1,7)/  2, 1, 1, 0, 0, 3, 1/
      data (ipa44( 571,j),j=1,7)/  2, 1, 1, 0, 1, 2, 2/
      data (ipa44( 572,j),j=1,7)/  2, 1, 1, 0, 2, 1, 2/
      data (ipa44( 573,j),j=1,7)/  2, 1, 1, 0, 3, 0, 2/
      data (ipa44( 574,j),j=1,7)/  2, 1, 1, 1, 1, 1, 1/
      data (ipa44( 575,j),j=1,7)/  2, 1, 1, 1, 2, 0, 2/
      data (ipa44( 576,j),j=1,7)/  2, 1, 2, 0, 0, 2, 2/
      data (ipa44( 577,j),j=1,7)/  2, 1, 2, 0, 1, 1, 2/
      data (ipa44( 578,j),j=1,7)/  2, 1, 2, 0, 2, 0, 2/
      data (ipa44( 579,j),j=1,7)/  2, 1, 2, 1, 0, 1, 2/
      data (ipa44( 580,j),j=1,7)/  2, 1, 2, 1, 1, 0, 2/
      data (ipa44( 581,j),j=1,7)/  2, 1, 2, 2, 0, 0, 2/
      data (ipa44( 582,j),j=1,7)/  2, 1, 3, 0, 0, 1, 2/
      data (ipa44( 583,j),j=1,7)/  2, 1, 3, 0, 1, 0, 2/
      data (ipa44( 584,j),j=1,7)/  2, 1, 3, 1, 0, 0, 2/
      data (ipa44( 585,j),j=1,7)/  2, 1, 4, 0, 0, 0, 2/
      data (ipa44( 586,j),j=1,7)/  2, 2, 2, 0, 0, 1, 1/
      data (ipa44( 587,j),j=1,7)/  2, 2, 2, 0, 1, 0, 2/
      data (ipa44( 588,j),j=1,7)/  2, 2, 3, 0, 0, 0, 2/
      data (ipa44( 589,j),j=1,7)/  3, 0, 0, 0, 1, 3, 2/
      data (ipa44( 590,j),j=1,7)/  3, 0, 0, 0, 2, 2, 2/
      data (ipa44( 591,j),j=1,7)/  3, 0, 0, 0, 3, 1, 2/
      data (ipa44( 592,j),j=1,7)/  3, 0, 0, 1, 1, 2, 1/
      data (ipa44( 593,j),j=1,7)/  3, 0, 0, 1, 2, 1, 2/
      data (ipa44( 594,j),j=1,7)/  3, 0, 0, 1, 3, 0, 2/
      data (ipa44( 595,j),j=1,7)/  3, 0, 0, 2, 2, 0, 1/
      data (ipa44( 596,j),j=1,7)/  3, 0, 1, 0, 0, 3, 2/
      data (ipa44( 597,j),j=1,7)/  3, 0, 1, 0, 1, 2, 2/
      data (ipa44( 598,j),j=1,7)/  3, 0, 1, 0, 2, 1, 2/
      data (ipa44( 599,j),j=1,7)/  3, 0, 1, 1, 0, 2, 2/
      data (ipa44( 600,j),j=1,7)/  3, 0, 1, 1, 1, 1, 2/
      data (ipa44( 601,j),j=1,7)/  3, 0, 1, 1, 2, 0, 2/
      data (ipa44( 602,j),j=1,7)/  3, 0, 1, 2, 0, 1, 2/
      data (ipa44( 603,j),j=1,7)/  3, 0, 1, 2, 1, 0, 2/
      data (ipa44( 604,j),j=1,7)/  3, 0, 1, 3, 0, 0, 2/
      data (ipa44( 605,j),j=1,7)/  3, 0, 2, 0, 0, 2, 2/
      data (ipa44( 606,j),j=1,7)/  3, 0, 2, 0, 1, 1, 2/
      data (ipa44( 607,j),j=1,7)/  3, 0, 2, 1, 0, 1, 2/
      data (ipa44( 608,j),j=1,7)/  3, 0, 2, 1, 1, 0, 2/
      data (ipa44( 609,j),j=1,7)/  3, 0, 2, 2, 0, 0, 2/
      data (ipa44( 610,j),j=1,7)/  3, 0, 3, 0, 0, 1, 2/
      data (ipa44( 611,j),j=1,7)/  3, 0, 3, 1, 0, 0, 2/
      data (ipa44( 612,j),j=1,7)/  3, 1, 1, 0, 0, 2, 1/
      data (ipa44( 613,j),j=1,7)/  3, 1, 1, 0, 1, 1, 2/
      data (ipa44( 614,j),j=1,7)/  3, 1, 1, 0, 2, 0, 2/
      data (ipa44( 615,j),j=1,7)/  3, 1, 1, 1, 1, 0, 1/
      data (ipa44( 616,j),j=1,7)/  3, 1, 2, 0, 0, 1, 2/
      data (ipa44( 617,j),j=1,7)/  3, 1, 2, 0, 1, 0, 2/
      data (ipa44( 618,j),j=1,7)/  3, 1, 2, 1, 0, 0, 2/
      data (ipa44( 619,j),j=1,7)/  3, 1, 3, 0, 0, 0, 2/
      data (ipa44( 620,j),j=1,7)/  3, 2, 2, 0, 0, 0, 1/
      data (ipa44( 621,j),j=1,7)/  4, 0, 0, 0, 1, 2, 2/
      data (ipa44( 622,j),j=1,7)/  4, 0, 0, 0, 2, 1, 2/
      data (ipa44( 623,j),j=1,7)/  4, 0, 0, 1, 1, 1, 1/
      data (ipa44( 624,j),j=1,7)/  4, 0, 0, 1, 2, 0, 2/
      data (ipa44( 625,j),j=1,7)/  4, 0, 1, 0, 0, 2, 2/
      data (ipa44( 626,j),j=1,7)/  4, 0, 1, 0, 1, 1, 2/
      data (ipa44( 627,j),j=1,7)/  4, 0, 1, 1, 0, 1, 2/
      data (ipa44( 628,j),j=1,7)/  4, 0, 1, 1, 1, 0, 2/
      data (ipa44( 629,j),j=1,7)/  4, 0, 1, 2, 0, 0, 2/
      data (ipa44( 630,j),j=1,7)/  4, 0, 2, 0, 0, 1, 2/
      data (ipa44( 631,j),j=1,7)/  4, 0, 2, 1, 0, 0, 2/
      data (ipa44( 632,j),j=1,7)/  4, 1, 1, 0, 0, 1, 1/
      data (ipa44( 633,j),j=1,7)/  4, 1, 1, 0, 1, 0, 2/
      data (ipa44( 634,j),j=1,7)/  4, 1, 2, 0, 0, 0, 2/
      data (ipa44( 635,j),j=1,7)/  5, 0, 0, 0, 1, 1, 2/
      data (ipa44( 636,j),j=1,7)/  5, 0, 0, 1, 1, 0, 1/
      data (ipa44( 637,j),j=1,7)/  5, 0, 1, 0, 0, 1, 2/
      data (ipa44( 638,j),j=1,7)/  5, 0, 1, 1, 0, 0, 2/
      data (ipa44( 639,j),j=1,7)/  5, 1, 1, 0, 0, 0, 1/
*********** 7 331*******************************
      data (nfas44(i),i=1,7)/0,0,9,43,129,308,639/
*********** 0 0*******************************
*********** 1 0*******************************
*********** 2 0*******************************
      data (ipa45(   1,j),j=1,6)/  0, 0, 1, 0, 1, 1/
      data (ipa45(   2,j),j=1,6)/  0, 0, 1, 1, 0, 1/
      data (ipa45(   3,j),j=1,6)/  0, 0, 1, 1, 1, 0/
      data (ipa45(   4,j),j=1,6)/  0, 1, 0, 0, 1, 1/
      data (ipa45(   5,j),j=1,6)/  0, 1, 0, 1, 0, 1/
      data (ipa45(   6,j),j=1,6)/  0, 1, 0, 1, 1, 0/
      data (ipa45(   7,j),j=1,6)/  0, 1, 1, 0, 1, 0/
      data (ipa45(   8,j),j=1,6)/  0, 1, 1, 1, 0, 0/
      data (ipa45(   9,j),j=1,6)/  1, 0, 0, 0, 1, 1/
      data (ipa45(  10,j),j=1,6)/  1, 0, 0, 1, 0, 1/
      data (ipa45(  11,j),j=1,6)/  1, 0, 0, 1, 1, 0/
      data (ipa45(  12,j),j=1,6)/  1, 0, 1, 0, 0, 1/
      data (ipa45(  13,j),j=1,6)/  1, 0, 1, 1, 0, 0/
      data (ipa45(  14,j),j=1,6)/  1, 1, 0, 0, 0, 1/
      data (ipa45(  15,j),j=1,6)/  1, 1, 0, 0, 1, 0/
      data (ipa45(  16,j),j=1,6)/  1, 1, 1, 0, 0, 0/
*********** 3 16*******************************
      data (ipa45(  17,j),j=1,6)/  0, 0, 1, 0, 1, 2/
      data (ipa45(  18,j),j=1,6)/  0, 0, 1, 0, 2, 1/
      data (ipa45(  19,j),j=1,6)/  0, 0, 1, 1, 0, 2/
      data (ipa45(  20,j),j=1,6)/  0, 0, 1, 1, 1, 1/
      data (ipa45(  21,j),j=1,6)/  0, 0, 1, 1, 2, 0/
      data (ipa45(  22,j),j=1,6)/  0, 0, 1, 2, 0, 1/
      data (ipa45(  23,j),j=1,6)/  0, 0, 1, 2, 1, 0/
      data (ipa45(  24,j),j=1,6)/  0, 0, 2, 0, 1, 1/
      data (ipa45(  25,j),j=1,6)/  0, 0, 2, 1, 0, 1/
      data (ipa45(  26,j),j=1,6)/  0, 0, 2, 1, 1, 0/
      data (ipa45(  27,j),j=1,6)/  0, 1, 0, 0, 1, 2/
      data (ipa45(  28,j),j=1,6)/  0, 1, 0, 0, 2, 1/
      data (ipa45(  29,j),j=1,6)/  0, 1, 0, 1, 0, 2/
      data (ipa45(  30,j),j=1,6)/  0, 1, 0, 1, 1, 1/
      data (ipa45(  31,j),j=1,6)/  0, 1, 0, 1, 2, 0/
      data (ipa45(  32,j),j=1,6)/  0, 1, 0, 2, 0, 1/
      data (ipa45(  33,j),j=1,6)/  0, 1, 0, 2, 1, 0/
      data (ipa45(  34,j),j=1,6)/  0, 1, 1, 0, 1, 1/
      data (ipa45(  35,j),j=1,6)/  0, 1, 1, 0, 2, 0/
      data (ipa45(  36,j),j=1,6)/  0, 1, 1, 1, 0, 1/
      data (ipa45(  37,j),j=1,6)/  0, 1, 1, 1, 1, 0/
      data (ipa45(  38,j),j=1,6)/  0, 1, 1, 2, 0, 0/
      data (ipa45(  39,j),j=1,6)/  0, 1, 2, 0, 1, 0/
      data (ipa45(  40,j),j=1,6)/  0, 1, 2, 1, 0, 0/
      data (ipa45(  41,j),j=1,6)/  0, 2, 0, 0, 1, 1/
      data (ipa45(  42,j),j=1,6)/  0, 2, 0, 1, 0, 1/
      data (ipa45(  43,j),j=1,6)/  0, 2, 0, 1, 1, 0/
      data (ipa45(  44,j),j=1,6)/  0, 2, 1, 0, 1, 0/
      data (ipa45(  45,j),j=1,6)/  0, 2, 1, 1, 0, 0/
      data (ipa45(  46,j),j=1,6)/  1, 0, 0, 0, 1, 2/
      data (ipa45(  47,j),j=1,6)/  1, 0, 0, 0, 2, 1/
      data (ipa45(  48,j),j=1,6)/  1, 0, 0, 1, 0, 2/
      data (ipa45(  49,j),j=1,6)/  1, 0, 0, 1, 1, 1/
      data (ipa45(  50,j),j=1,6)/  1, 0, 0, 1, 2, 0/
      data (ipa45(  51,j),j=1,6)/  1, 0, 0, 2, 0, 1/
      data (ipa45(  52,j),j=1,6)/  1, 0, 0, 2, 1, 0/
      data (ipa45(  53,j),j=1,6)/  1, 0, 1, 0, 0, 2/
      data (ipa45(  54,j),j=1,6)/  1, 0, 1, 0, 1, 1/
      data (ipa45(  55,j),j=1,6)/  1, 0, 1, 1, 0, 1/
      data (ipa45(  56,j),j=1,6)/  1, 0, 1, 1, 1, 0/
      data (ipa45(  57,j),j=1,6)/  1, 0, 1, 2, 0, 0/
      data (ipa45(  58,j),j=1,6)/  1, 0, 2, 0, 0, 1/
      data (ipa45(  59,j),j=1,6)/  1, 0, 2, 1, 0, 0/
      data (ipa45(  60,j),j=1,6)/  1, 1, 0, 0, 0, 2/
      data (ipa45(  61,j),j=1,6)/  1, 1, 0, 0, 1, 1/
      data (ipa45(  62,j),j=1,6)/  1, 1, 0, 0, 2, 0/
      data (ipa45(  63,j),j=1,6)/  1, 1, 0, 1, 0, 1/
      data (ipa45(  64,j),j=1,6)/  1, 1, 0, 1, 1, 0/
      data (ipa45(  65,j),j=1,6)/  1, 1, 1, 0, 0, 1/
      data (ipa45(  66,j),j=1,6)/  1, 1, 1, 0, 1, 0/
      data (ipa45(  67,j),j=1,6)/  1, 1, 1, 1, 0, 0/
      data (ipa45(  68,j),j=1,6)/  1, 1, 2, 0, 0, 0/
      data (ipa45(  69,j),j=1,6)/  1, 2, 0, 0, 0, 1/
      data (ipa45(  70,j),j=1,6)/  1, 2, 0, 0, 1, 0/
      data (ipa45(  71,j),j=1,6)/  1, 2, 1, 0, 0, 0/
      data (ipa45(  72,j),j=1,6)/  2, 0, 0, 0, 1, 1/
      data (ipa45(  73,j),j=1,6)/  2, 0, 0, 1, 0, 1/
      data (ipa45(  74,j),j=1,6)/  2, 0, 0, 1, 1, 0/
      data (ipa45(  75,j),j=1,6)/  2, 0, 1, 0, 0, 1/
      data (ipa45(  76,j),j=1,6)/  2, 0, 1, 1, 0, 0/
      data (ipa45(  77,j),j=1,6)/  2, 1, 0, 0, 0, 1/
      data (ipa45(  78,j),j=1,6)/  2, 1, 0, 0, 1, 0/
      data (ipa45(  79,j),j=1,6)/  2, 1, 1, 0, 0, 0/
*********** 4 63*******************************
      data (ipa45(  80,j),j=1,6)/  0, 0, 1, 0, 1, 3/
      data (ipa45(  81,j),j=1,6)/  0, 0, 1, 0, 2, 2/
      data (ipa45(  82,j),j=1,6)/  0, 0, 1, 0, 3, 1/
      data (ipa45(  83,j),j=1,6)/  0, 0, 1, 1, 0, 3/
      data (ipa45(  84,j),j=1,6)/  0, 0, 1, 1, 1, 2/
      data (ipa45(  85,j),j=1,6)/  0, 0, 1, 1, 2, 1/
      data (ipa45(  86,j),j=1,6)/  0, 0, 1, 1, 3, 0/
      data (ipa45(  87,j),j=1,6)/  0, 0, 1, 2, 0, 2/
      data (ipa45(  88,j),j=1,6)/  0, 0, 1, 2, 1, 1/
      data (ipa45(  89,j),j=1,6)/  0, 0, 1, 2, 2, 0/
      data (ipa45(  90,j),j=1,6)/  0, 0, 1, 3, 0, 1/
      data (ipa45(  91,j),j=1,6)/  0, 0, 1, 3, 1, 0/
      data (ipa45(  92,j),j=1,6)/  0, 0, 2, 0, 1, 2/
      data (ipa45(  93,j),j=1,6)/  0, 0, 2, 0, 2, 1/
      data (ipa45(  94,j),j=1,6)/  0, 0, 2, 1, 0, 2/
      data (ipa45(  95,j),j=1,6)/  0, 0, 2, 1, 1, 1/
      data (ipa45(  96,j),j=1,6)/  0, 0, 2, 1, 2, 0/
      data (ipa45(  97,j),j=1,6)/  0, 0, 2, 2, 0, 1/
      data (ipa45(  98,j),j=1,6)/  0, 0, 2, 2, 1, 0/
      data (ipa45(  99,j),j=1,6)/  0, 0, 3, 0, 1, 1/
      data (ipa45( 100,j),j=1,6)/  0, 0, 3, 1, 0, 1/
      data (ipa45( 101,j),j=1,6)/  0, 0, 3, 1, 1, 0/
      data (ipa45( 102,j),j=1,6)/  0, 1, 0, 0, 1, 3/
      data (ipa45( 103,j),j=1,6)/  0, 1, 0, 0, 2, 2/
      data (ipa45( 104,j),j=1,6)/  0, 1, 0, 0, 3, 1/
      data (ipa45( 105,j),j=1,6)/  0, 1, 0, 1, 0, 3/
      data (ipa45( 106,j),j=1,6)/  0, 1, 0, 1, 1, 2/
      data (ipa45( 107,j),j=1,6)/  0, 1, 0, 1, 2, 1/
      data (ipa45( 108,j),j=1,6)/  0, 1, 0, 1, 3, 0/
      data (ipa45( 109,j),j=1,6)/  0, 1, 0, 2, 0, 2/
      data (ipa45( 110,j),j=1,6)/  0, 1, 0, 2, 1, 1/
      data (ipa45( 111,j),j=1,6)/  0, 1, 0, 2, 2, 0/
      data (ipa45( 112,j),j=1,6)/  0, 1, 0, 3, 0, 1/
      data (ipa45( 113,j),j=1,6)/  0, 1, 0, 3, 1, 0/
      data (ipa45( 114,j),j=1,6)/  0, 1, 1, 0, 1, 2/
      data (ipa45( 115,j),j=1,6)/  0, 1, 1, 0, 2, 1/
      data (ipa45( 116,j),j=1,6)/  0, 1, 1, 0, 3, 0/
      data (ipa45( 117,j),j=1,6)/  0, 1, 1, 1, 0, 2/
      data (ipa45( 118,j),j=1,6)/  0, 1, 1, 1, 1, 1/
      data (ipa45( 119,j),j=1,6)/  0, 1, 1, 1, 2, 0/
      data (ipa45( 120,j),j=1,6)/  0, 1, 1, 2, 0, 1/
      data (ipa45( 121,j),j=1,6)/  0, 1, 1, 2, 1, 0/
      data (ipa45( 122,j),j=1,6)/  0, 1, 1, 3, 0, 0/
      data (ipa45( 123,j),j=1,6)/  0, 1, 2, 0, 1, 1/
      data (ipa45( 124,j),j=1,6)/  0, 1, 2, 0, 2, 0/
      data (ipa45( 125,j),j=1,6)/  0, 1, 2, 1, 0, 1/
      data (ipa45( 126,j),j=1,6)/  0, 1, 2, 1, 1, 0/
      data (ipa45( 127,j),j=1,6)/  0, 1, 2, 2, 0, 0/
      data (ipa45( 128,j),j=1,6)/  0, 1, 3, 0, 1, 0/
      data (ipa45( 129,j),j=1,6)/  0, 1, 3, 1, 0, 0/
      data (ipa45( 130,j),j=1,6)/  0, 2, 0, 0, 1, 2/
      data (ipa45( 131,j),j=1,6)/  0, 2, 0, 0, 2, 1/
      data (ipa45( 132,j),j=1,6)/  0, 2, 0, 1, 0, 2/
      data (ipa45( 133,j),j=1,6)/  0, 2, 0, 1, 1, 1/
      data (ipa45( 134,j),j=1,6)/  0, 2, 0, 1, 2, 0/
      data (ipa45( 135,j),j=1,6)/  0, 2, 0, 2, 0, 1/
      data (ipa45( 136,j),j=1,6)/  0, 2, 0, 2, 1, 0/
      data (ipa45( 137,j),j=1,6)/  0, 2, 1, 0, 1, 1/
      data (ipa45( 138,j),j=1,6)/  0, 2, 1, 0, 2, 0/
      data (ipa45( 139,j),j=1,6)/  0, 2, 1, 1, 0, 1/
      data (ipa45( 140,j),j=1,6)/  0, 2, 1, 1, 1, 0/
      data (ipa45( 141,j),j=1,6)/  0, 2, 1, 2, 0, 0/
      data (ipa45( 142,j),j=1,6)/  0, 2, 2, 0, 1, 0/
      data (ipa45( 143,j),j=1,6)/  0, 2, 2, 1, 0, 0/
      data (ipa45( 144,j),j=1,6)/  0, 3, 0, 0, 1, 1/
      data (ipa45( 145,j),j=1,6)/  0, 3, 0, 1, 0, 1/
      data (ipa45( 146,j),j=1,6)/  0, 3, 0, 1, 1, 0/
      data (ipa45( 147,j),j=1,6)/  0, 3, 1, 0, 1, 0/
      data (ipa45( 148,j),j=1,6)/  0, 3, 1, 1, 0, 0/
      data (ipa45( 149,j),j=1,6)/  1, 0, 0, 0, 1, 3/
      data (ipa45( 150,j),j=1,6)/  1, 0, 0, 0, 2, 2/
      data (ipa45( 151,j),j=1,6)/  1, 0, 0, 0, 3, 1/
      data (ipa45( 152,j),j=1,6)/  1, 0, 0, 1, 0, 3/
      data (ipa45( 153,j),j=1,6)/  1, 0, 0, 1, 1, 2/
      data (ipa45( 154,j),j=1,6)/  1, 0, 0, 1, 2, 1/
      data (ipa45( 155,j),j=1,6)/  1, 0, 0, 1, 3, 0/
      data (ipa45( 156,j),j=1,6)/  1, 0, 0, 2, 0, 2/
      data (ipa45( 157,j),j=1,6)/  1, 0, 0, 2, 1, 1/
      data (ipa45( 158,j),j=1,6)/  1, 0, 0, 2, 2, 0/
      data (ipa45( 159,j),j=1,6)/  1, 0, 0, 3, 0, 1/
      data (ipa45( 160,j),j=1,6)/  1, 0, 0, 3, 1, 0/
      data (ipa45( 161,j),j=1,6)/  1, 0, 1, 0, 0, 3/
      data (ipa45( 162,j),j=1,6)/  1, 0, 1, 0, 1, 2/
      data (ipa45( 163,j),j=1,6)/  1, 0, 1, 0, 2, 1/
      data (ipa45( 164,j),j=1,6)/  1, 0, 1, 1, 0, 2/
      data (ipa45( 165,j),j=1,6)/  1, 0, 1, 1, 1, 1/
      data (ipa45( 166,j),j=1,6)/  1, 0, 1, 1, 2, 0/
      data (ipa45( 167,j),j=1,6)/  1, 0, 1, 2, 0, 1/
      data (ipa45( 168,j),j=1,6)/  1, 0, 1, 2, 1, 0/
      data (ipa45( 169,j),j=1,6)/  1, 0, 1, 3, 0, 0/
      data (ipa45( 170,j),j=1,6)/  1, 0, 2, 0, 0, 2/
      data (ipa45( 171,j),j=1,6)/  1, 0, 2, 0, 1, 1/
      data (ipa45( 172,j),j=1,6)/  1, 0, 2, 1, 0, 1/
      data (ipa45( 173,j),j=1,6)/  1, 0, 2, 1, 1, 0/
      data (ipa45( 174,j),j=1,6)/  1, 0, 2, 2, 0, 0/
      data (ipa45( 175,j),j=1,6)/  1, 0, 3, 0, 0, 1/
      data (ipa45( 176,j),j=1,6)/  1, 0, 3, 1, 0, 0/
      data (ipa45( 177,j),j=1,6)/  1, 1, 0, 0, 0, 3/
      data (ipa45( 178,j),j=1,6)/  1, 1, 0, 0, 1, 2/
      data (ipa45( 179,j),j=1,6)/  1, 1, 0, 0, 2, 1/
      data (ipa45( 180,j),j=1,6)/  1, 1, 0, 0, 3, 0/
      data (ipa45( 181,j),j=1,6)/  1, 1, 0, 1, 0, 2/
      data (ipa45( 182,j),j=1,6)/  1, 1, 0, 1, 1, 1/
      data (ipa45( 183,j),j=1,6)/  1, 1, 0, 1, 2, 0/
      data (ipa45( 184,j),j=1,6)/  1, 1, 0, 2, 0, 1/
      data (ipa45( 185,j),j=1,6)/  1, 1, 0, 2, 1, 0/
      data (ipa45( 186,j),j=1,6)/  1, 1, 1, 0, 0, 2/
      data (ipa45( 187,j),j=1,6)/  1, 1, 1, 0, 1, 1/
      data (ipa45( 188,j),j=1,6)/  1, 1, 1, 0, 2, 0/
      data (ipa45( 189,j),j=1,6)/  1, 1, 1, 1, 0, 1/
      data (ipa45( 190,j),j=1,6)/  1, 1, 1, 1, 1, 0/
      data (ipa45( 191,j),j=1,6)/  1, 1, 1, 2, 0, 0/
      data (ipa45( 192,j),j=1,6)/  1, 1, 2, 0, 0, 1/
      data (ipa45( 193,j),j=1,6)/  1, 1, 2, 0, 1, 0/
      data (ipa45( 194,j),j=1,6)/  1, 1, 2, 1, 0, 0/
      data (ipa45( 195,j),j=1,6)/  1, 1, 3, 0, 0, 0/
      data (ipa45( 196,j),j=1,6)/  1, 2, 0, 0, 0, 2/
      data (ipa45( 197,j),j=1,6)/  1, 2, 0, 0, 1, 1/
      data (ipa45( 198,j),j=1,6)/  1, 2, 0, 0, 2, 0/
      data (ipa45( 199,j),j=1,6)/  1, 2, 0, 1, 0, 1/
      data (ipa45( 200,j),j=1,6)/  1, 2, 0, 1, 1, 0/
      data (ipa45( 201,j),j=1,6)/  1, 2, 1, 0, 0, 1/
      data (ipa45( 202,j),j=1,6)/  1, 2, 1, 0, 1, 0/
      data (ipa45( 203,j),j=1,6)/  1, 2, 1, 1, 0, 0/
      data (ipa45( 204,j),j=1,6)/  1, 2, 2, 0, 0, 0/
      data (ipa45( 205,j),j=1,6)/  1, 3, 0, 0, 0, 1/
      data (ipa45( 206,j),j=1,6)/  1, 3, 0, 0, 1, 0/
      data (ipa45( 207,j),j=1,6)/  1, 3, 1, 0, 0, 0/
      data (ipa45( 208,j),j=1,6)/  2, 0, 0, 0, 1, 2/
      data (ipa45( 209,j),j=1,6)/  2, 0, 0, 0, 2, 1/
      data (ipa45( 210,j),j=1,6)/  2, 0, 0, 1, 0, 2/
      data (ipa45( 211,j),j=1,6)/  2, 0, 0, 1, 1, 1/
      data (ipa45( 212,j),j=1,6)/  2, 0, 0, 1, 2, 0/
      data (ipa45( 213,j),j=1,6)/  2, 0, 0, 2, 0, 1/
      data (ipa45( 214,j),j=1,6)/  2, 0, 0, 2, 1, 0/
      data (ipa45( 215,j),j=1,6)/  2, 0, 1, 0, 0, 2/
      data (ipa45( 216,j),j=1,6)/  2, 0, 1, 0, 1, 1/
      data (ipa45( 217,j),j=1,6)/  2, 0, 1, 1, 0, 1/
      data (ipa45( 218,j),j=1,6)/  2, 0, 1, 1, 1, 0/
      data (ipa45( 219,j),j=1,6)/  2, 0, 1, 2, 0, 0/
      data (ipa45( 220,j),j=1,6)/  2, 0, 2, 0, 0, 1/
      data (ipa45( 221,j),j=1,6)/  2, 0, 2, 1, 0, 0/
      data (ipa45( 222,j),j=1,6)/  2, 1, 0, 0, 0, 2/
      data (ipa45( 223,j),j=1,6)/  2, 1, 0, 0, 1, 1/
      data (ipa45( 224,j),j=1,6)/  2, 1, 0, 0, 2, 0/
      data (ipa45( 225,j),j=1,6)/  2, 1, 0, 1, 0, 1/
      data (ipa45( 226,j),j=1,6)/  2, 1, 0, 1, 1, 0/
      data (ipa45( 227,j),j=1,6)/  2, 1, 1, 0, 0, 1/
      data (ipa45( 228,j),j=1,6)/  2, 1, 1, 0, 1, 0/
      data (ipa45( 229,j),j=1,6)/  2, 1, 1, 1, 0, 0/
      data (ipa45( 230,j),j=1,6)/  2, 1, 2, 0, 0, 0/
      data (ipa45( 231,j),j=1,6)/  2, 2, 0, 0, 0, 1/
      data (ipa45( 232,j),j=1,6)/  2, 2, 0, 0, 1, 0/
      data (ipa45( 233,j),j=1,6)/  2, 2, 1, 0, 0, 0/
      data (ipa45( 234,j),j=1,6)/  3, 0, 0, 0, 1, 1/
      data (ipa45( 235,j),j=1,6)/  3, 0, 0, 1, 0, 1/
      data (ipa45( 236,j),j=1,6)/  3, 0, 0, 1, 1, 0/
      data (ipa45( 237,j),j=1,6)/  3, 0, 1, 0, 0, 1/
      data (ipa45( 238,j),j=1,6)/  3, 0, 1, 1, 0, 0/
      data (ipa45( 239,j),j=1,6)/  3, 1, 0, 0, 0, 1/
      data (ipa45( 240,j),j=1,6)/  3, 1, 0, 0, 1, 0/
      data (ipa45( 241,j),j=1,6)/  3, 1, 1, 0, 0, 0/
*********** 5 162*******************************
      data (ipa45( 242,j),j=1,6)/  0, 0, 1, 0, 1, 4/
      data (ipa45( 243,j),j=1,6)/  0, 0, 1, 0, 2, 3/
      data (ipa45( 244,j),j=1,6)/  0, 0, 1, 0, 3, 2/
      data (ipa45( 245,j),j=1,6)/  0, 0, 1, 0, 4, 1/
      data (ipa45( 246,j),j=1,6)/  0, 0, 1, 1, 0, 4/
      data (ipa45( 247,j),j=1,6)/  0, 0, 1, 1, 1, 3/
      data (ipa45( 248,j),j=1,6)/  0, 0, 1, 1, 2, 2/
      data (ipa45( 249,j),j=1,6)/  0, 0, 1, 1, 3, 1/
      data (ipa45( 250,j),j=1,6)/  0, 0, 1, 1, 4, 0/
      data (ipa45( 251,j),j=1,6)/  0, 0, 1, 2, 0, 3/
      data (ipa45( 252,j),j=1,6)/  0, 0, 1, 2, 1, 2/
      data (ipa45( 253,j),j=1,6)/  0, 0, 1, 2, 2, 1/
      data (ipa45( 254,j),j=1,6)/  0, 0, 1, 2, 3, 0/
      data (ipa45( 255,j),j=1,6)/  0, 0, 1, 3, 0, 2/
      data (ipa45( 256,j),j=1,6)/  0, 0, 1, 3, 1, 1/
      data (ipa45( 257,j),j=1,6)/  0, 0, 1, 3, 2, 0/
      data (ipa45( 258,j),j=1,6)/  0, 0, 1, 4, 0, 1/
      data (ipa45( 259,j),j=1,6)/  0, 0, 1, 4, 1, 0/
      data (ipa45( 260,j),j=1,6)/  0, 0, 2, 0, 1, 3/
      data (ipa45( 261,j),j=1,6)/  0, 0, 2, 0, 2, 2/
      data (ipa45( 262,j),j=1,6)/  0, 0, 2, 0, 3, 1/
      data (ipa45( 263,j),j=1,6)/  0, 0, 2, 1, 0, 3/
      data (ipa45( 264,j),j=1,6)/  0, 0, 2, 1, 1, 2/
      data (ipa45( 265,j),j=1,6)/  0, 0, 2, 1, 2, 1/
      data (ipa45( 266,j),j=1,6)/  0, 0, 2, 1, 3, 0/
      data (ipa45( 267,j),j=1,6)/  0, 0, 2, 2, 0, 2/
      data (ipa45( 268,j),j=1,6)/  0, 0, 2, 2, 1, 1/
      data (ipa45( 269,j),j=1,6)/  0, 0, 2, 2, 2, 0/
      data (ipa45( 270,j),j=1,6)/  0, 0, 2, 3, 0, 1/
      data (ipa45( 271,j),j=1,6)/  0, 0, 2, 3, 1, 0/
      data (ipa45( 272,j),j=1,6)/  0, 0, 3, 0, 1, 2/
      data (ipa45( 273,j),j=1,6)/  0, 0, 3, 0, 2, 1/
      data (ipa45( 274,j),j=1,6)/  0, 0, 3, 1, 0, 2/
      data (ipa45( 275,j),j=1,6)/  0, 0, 3, 1, 1, 1/
      data (ipa45( 276,j),j=1,6)/  0, 0, 3, 1, 2, 0/
      data (ipa45( 277,j),j=1,6)/  0, 0, 3, 2, 0, 1/
      data (ipa45( 278,j),j=1,6)/  0, 0, 3, 2, 1, 0/
      data (ipa45( 279,j),j=1,6)/  0, 0, 4, 0, 1, 1/
      data (ipa45( 280,j),j=1,6)/  0, 0, 4, 1, 0, 1/
      data (ipa45( 281,j),j=1,6)/  0, 0, 4, 1, 1, 0/
      data (ipa45( 282,j),j=1,6)/  0, 1, 0, 0, 1, 4/
      data (ipa45( 283,j),j=1,6)/  0, 1, 0, 0, 2, 3/
      data (ipa45( 284,j),j=1,6)/  0, 1, 0, 0, 3, 2/
      data (ipa45( 285,j),j=1,6)/  0, 1, 0, 0, 4, 1/
      data (ipa45( 286,j),j=1,6)/  0, 1, 0, 1, 0, 4/
      data (ipa45( 287,j),j=1,6)/  0, 1, 0, 1, 1, 3/
      data (ipa45( 288,j),j=1,6)/  0, 1, 0, 1, 2, 2/
      data (ipa45( 289,j),j=1,6)/  0, 1, 0, 1, 3, 1/
      data (ipa45( 290,j),j=1,6)/  0, 1, 0, 1, 4, 0/
      data (ipa45( 291,j),j=1,6)/  0, 1, 0, 2, 0, 3/
      data (ipa45( 292,j),j=1,6)/  0, 1, 0, 2, 1, 2/
      data (ipa45( 293,j),j=1,6)/  0, 1, 0, 2, 2, 1/
      data (ipa45( 294,j),j=1,6)/  0, 1, 0, 2, 3, 0/
      data (ipa45( 295,j),j=1,6)/  0, 1, 0, 3, 0, 2/
      data (ipa45( 296,j),j=1,6)/  0, 1, 0, 3, 1, 1/
      data (ipa45( 297,j),j=1,6)/  0, 1, 0, 3, 2, 0/
      data (ipa45( 298,j),j=1,6)/  0, 1, 0, 4, 0, 1/
      data (ipa45( 299,j),j=1,6)/  0, 1, 0, 4, 1, 0/
      data (ipa45( 300,j),j=1,6)/  0, 1, 1, 0, 1, 3/
      data (ipa45( 301,j),j=1,6)/  0, 1, 1, 0, 2, 2/
      data (ipa45( 302,j),j=1,6)/  0, 1, 1, 0, 3, 1/
      data (ipa45( 303,j),j=1,6)/  0, 1, 1, 0, 4, 0/
      data (ipa45( 304,j),j=1,6)/  0, 1, 1, 1, 0, 3/
      data (ipa45( 305,j),j=1,6)/  0, 1, 1, 1, 1, 2/
      data (ipa45( 306,j),j=1,6)/  0, 1, 1, 1, 2, 1/
      data (ipa45( 307,j),j=1,6)/  0, 1, 1, 1, 3, 0/
      data (ipa45( 308,j),j=1,6)/  0, 1, 1, 2, 0, 2/
      data (ipa45( 309,j),j=1,6)/  0, 1, 1, 2, 1, 1/
      data (ipa45( 310,j),j=1,6)/  0, 1, 1, 2, 2, 0/
      data (ipa45( 311,j),j=1,6)/  0, 1, 1, 3, 0, 1/
      data (ipa45( 312,j),j=1,6)/  0, 1, 1, 3, 1, 0/
      data (ipa45( 313,j),j=1,6)/  0, 1, 1, 4, 0, 0/
      data (ipa45( 314,j),j=1,6)/  0, 1, 2, 0, 1, 2/
      data (ipa45( 315,j),j=1,6)/  0, 1, 2, 0, 2, 1/
      data (ipa45( 316,j),j=1,6)/  0, 1, 2, 0, 3, 0/
      data (ipa45( 317,j),j=1,6)/  0, 1, 2, 1, 0, 2/
      data (ipa45( 318,j),j=1,6)/  0, 1, 2, 1, 1, 1/
      data (ipa45( 319,j),j=1,6)/  0, 1, 2, 1, 2, 0/
      data (ipa45( 320,j),j=1,6)/  0, 1, 2, 2, 0, 1/
      data (ipa45( 321,j),j=1,6)/  0, 1, 2, 2, 1, 0/
      data (ipa45( 322,j),j=1,6)/  0, 1, 2, 3, 0, 0/
      data (ipa45( 323,j),j=1,6)/  0, 1, 3, 0, 1, 1/
      data (ipa45( 324,j),j=1,6)/  0, 1, 3, 0, 2, 0/
      data (ipa45( 325,j),j=1,6)/  0, 1, 3, 1, 0, 1/
      data (ipa45( 326,j),j=1,6)/  0, 1, 3, 1, 1, 0/
      data (ipa45( 327,j),j=1,6)/  0, 1, 3, 2, 0, 0/
      data (ipa45( 328,j),j=1,6)/  0, 1, 4, 0, 1, 0/
      data (ipa45( 329,j),j=1,6)/  0, 1, 4, 1, 0, 0/
      data (ipa45( 330,j),j=1,6)/  0, 2, 0, 0, 1, 3/
      data (ipa45( 331,j),j=1,6)/  0, 2, 0, 0, 2, 2/
      data (ipa45( 332,j),j=1,6)/  0, 2, 0, 0, 3, 1/
      data (ipa45( 333,j),j=1,6)/  0, 2, 0, 1, 0, 3/
      data (ipa45( 334,j),j=1,6)/  0, 2, 0, 1, 1, 2/
      data (ipa45( 335,j),j=1,6)/  0, 2, 0, 1, 2, 1/
      data (ipa45( 336,j),j=1,6)/  0, 2, 0, 1, 3, 0/
      data (ipa45( 337,j),j=1,6)/  0, 2, 0, 2, 0, 2/
      data (ipa45( 338,j),j=1,6)/  0, 2, 0, 2, 1, 1/
      data (ipa45( 339,j),j=1,6)/  0, 2, 0, 2, 2, 0/
      data (ipa45( 340,j),j=1,6)/  0, 2, 0, 3, 0, 1/
      data (ipa45( 341,j),j=1,6)/  0, 2, 0, 3, 1, 0/
      data (ipa45( 342,j),j=1,6)/  0, 2, 1, 0, 1, 2/
      data (ipa45( 343,j),j=1,6)/  0, 2, 1, 0, 2, 1/
      data (ipa45( 344,j),j=1,6)/  0, 2, 1, 0, 3, 0/
      data (ipa45( 345,j),j=1,6)/  0, 2, 1, 1, 0, 2/
      data (ipa45( 346,j),j=1,6)/  0, 2, 1, 1, 1, 1/
      data (ipa45( 347,j),j=1,6)/  0, 2, 1, 1, 2, 0/
      data (ipa45( 348,j),j=1,6)/  0, 2, 1, 2, 0, 1/
      data (ipa45( 349,j),j=1,6)/  0, 2, 1, 2, 1, 0/
      data (ipa45( 350,j),j=1,6)/  0, 2, 1, 3, 0, 0/
      data (ipa45( 351,j),j=1,6)/  0, 2, 2, 0, 1, 1/
      data (ipa45( 352,j),j=1,6)/  0, 2, 2, 0, 2, 0/
      data (ipa45( 353,j),j=1,6)/  0, 2, 2, 1, 0, 1/
      data (ipa45( 354,j),j=1,6)/  0, 2, 2, 1, 1, 0/
      data (ipa45( 355,j),j=1,6)/  0, 2, 2, 2, 0, 0/
      data (ipa45( 356,j),j=1,6)/  0, 2, 3, 0, 1, 0/
      data (ipa45( 357,j),j=1,6)/  0, 2, 3, 1, 0, 0/
      data (ipa45( 358,j),j=1,6)/  0, 3, 0, 0, 1, 2/
      data (ipa45( 359,j),j=1,6)/  0, 3, 0, 0, 2, 1/
      data (ipa45( 360,j),j=1,6)/  0, 3, 0, 1, 0, 2/
      data (ipa45( 361,j),j=1,6)/  0, 3, 0, 1, 1, 1/
      data (ipa45( 362,j),j=1,6)/  0, 3, 0, 1, 2, 0/
      data (ipa45( 363,j),j=1,6)/  0, 3, 0, 2, 0, 1/
      data (ipa45( 364,j),j=1,6)/  0, 3, 0, 2, 1, 0/
      data (ipa45( 365,j),j=1,6)/  0, 3, 1, 0, 1, 1/
      data (ipa45( 366,j),j=1,6)/  0, 3, 1, 0, 2, 0/
      data (ipa45( 367,j),j=1,6)/  0, 3, 1, 1, 0, 1/
      data (ipa45( 368,j),j=1,6)/  0, 3, 1, 1, 1, 0/
      data (ipa45( 369,j),j=1,6)/  0, 3, 1, 2, 0, 0/
      data (ipa45( 370,j),j=1,6)/  0, 3, 2, 0, 1, 0/
      data (ipa45( 371,j),j=1,6)/  0, 3, 2, 1, 0, 0/
      data (ipa45( 372,j),j=1,6)/  0, 4, 0, 0, 1, 1/
      data (ipa45( 373,j),j=1,6)/  0, 4, 0, 1, 0, 1/
      data (ipa45( 374,j),j=1,6)/  0, 4, 0, 1, 1, 0/
      data (ipa45( 375,j),j=1,6)/  0, 4, 1, 0, 1, 0/
      data (ipa45( 376,j),j=1,6)/  0, 4, 1, 1, 0, 0/
      data (ipa45( 377,j),j=1,6)/  1, 0, 0, 0, 1, 4/
      data (ipa45( 378,j),j=1,6)/  1, 0, 0, 0, 2, 3/
      data (ipa45( 379,j),j=1,6)/  1, 0, 0, 0, 3, 2/
      data (ipa45( 380,j),j=1,6)/  1, 0, 0, 0, 4, 1/
      data (ipa45( 381,j),j=1,6)/  1, 0, 0, 1, 0, 4/
      data (ipa45( 382,j),j=1,6)/  1, 0, 0, 1, 1, 3/
      data (ipa45( 383,j),j=1,6)/  1, 0, 0, 1, 2, 2/
      data (ipa45( 384,j),j=1,6)/  1, 0, 0, 1, 3, 1/
      data (ipa45( 385,j),j=1,6)/  1, 0, 0, 1, 4, 0/
      data (ipa45( 386,j),j=1,6)/  1, 0, 0, 2, 0, 3/
      data (ipa45( 387,j),j=1,6)/  1, 0, 0, 2, 1, 2/
      data (ipa45( 388,j),j=1,6)/  1, 0, 0, 2, 2, 1/
      data (ipa45( 389,j),j=1,6)/  1, 0, 0, 2, 3, 0/
      data (ipa45( 390,j),j=1,6)/  1, 0, 0, 3, 0, 2/
      data (ipa45( 391,j),j=1,6)/  1, 0, 0, 3, 1, 1/
      data (ipa45( 392,j),j=1,6)/  1, 0, 0, 3, 2, 0/
      data (ipa45( 393,j),j=1,6)/  1, 0, 0, 4, 0, 1/
      data (ipa45( 394,j),j=1,6)/  1, 0, 0, 4, 1, 0/
      data (ipa45( 395,j),j=1,6)/  1, 0, 1, 0, 0, 4/
      data (ipa45( 396,j),j=1,6)/  1, 0, 1, 0, 1, 3/
      data (ipa45( 397,j),j=1,6)/  1, 0, 1, 0, 2, 2/
      data (ipa45( 398,j),j=1,6)/  1, 0, 1, 0, 3, 1/
      data (ipa45( 399,j),j=1,6)/  1, 0, 1, 1, 0, 3/
      data (ipa45( 400,j),j=1,6)/  1, 0, 1, 1, 1, 2/
      data (ipa45( 401,j),j=1,6)/  1, 0, 1, 1, 2, 1/
      data (ipa45( 402,j),j=1,6)/  1, 0, 1, 1, 3, 0/
      data (ipa45( 403,j),j=1,6)/  1, 0, 1, 2, 0, 2/
      data (ipa45( 404,j),j=1,6)/  1, 0, 1, 2, 1, 1/
      data (ipa45( 405,j),j=1,6)/  1, 0, 1, 2, 2, 0/
      data (ipa45( 406,j),j=1,6)/  1, 0, 1, 3, 0, 1/
      data (ipa45( 407,j),j=1,6)/  1, 0, 1, 3, 1, 0/
      data (ipa45( 408,j),j=1,6)/  1, 0, 1, 4, 0, 0/
      data (ipa45( 409,j),j=1,6)/  1, 0, 2, 0, 0, 3/
      data (ipa45( 410,j),j=1,6)/  1, 0, 2, 0, 1, 2/
      data (ipa45( 411,j),j=1,6)/  1, 0, 2, 0, 2, 1/
      data (ipa45( 412,j),j=1,6)/  1, 0, 2, 1, 0, 2/
      data (ipa45( 413,j),j=1,6)/  1, 0, 2, 1, 1, 1/
      data (ipa45( 414,j),j=1,6)/  1, 0, 2, 1, 2, 0/
      data (ipa45( 415,j),j=1,6)/  1, 0, 2, 2, 0, 1/
      data (ipa45( 416,j),j=1,6)/  1, 0, 2, 2, 1, 0/
      data (ipa45( 417,j),j=1,6)/  1, 0, 2, 3, 0, 0/
      data (ipa45( 418,j),j=1,6)/  1, 0, 3, 0, 0, 2/
      data (ipa45( 419,j),j=1,6)/  1, 0, 3, 0, 1, 1/
      data (ipa45( 420,j),j=1,6)/  1, 0, 3, 1, 0, 1/
      data (ipa45( 421,j),j=1,6)/  1, 0, 3, 1, 1, 0/
      data (ipa45( 422,j),j=1,6)/  1, 0, 3, 2, 0, 0/
      data (ipa45( 423,j),j=1,6)/  1, 0, 4, 0, 0, 1/
      data (ipa45( 424,j),j=1,6)/  1, 0, 4, 1, 0, 0/
      data (ipa45( 425,j),j=1,6)/  1, 1, 0, 0, 0, 4/
      data (ipa45( 426,j),j=1,6)/  1, 1, 0, 0, 1, 3/
      data (ipa45( 427,j),j=1,6)/  1, 1, 0, 0, 2, 2/
      data (ipa45( 428,j),j=1,6)/  1, 1, 0, 0, 3, 1/
      data (ipa45( 429,j),j=1,6)/  1, 1, 0, 0, 4, 0/
      data (ipa45( 430,j),j=1,6)/  1, 1, 0, 1, 0, 3/
      data (ipa45( 431,j),j=1,6)/  1, 1, 0, 1, 1, 2/
      data (ipa45( 432,j),j=1,6)/  1, 1, 0, 1, 2, 1/
      data (ipa45( 433,j),j=1,6)/  1, 1, 0, 1, 3, 0/
      data (ipa45( 434,j),j=1,6)/  1, 1, 0, 2, 0, 2/
      data (ipa45( 435,j),j=1,6)/  1, 1, 0, 2, 1, 1/
      data (ipa45( 436,j),j=1,6)/  1, 1, 0, 2, 2, 0/
      data (ipa45( 437,j),j=1,6)/  1, 1, 0, 3, 0, 1/
      data (ipa45( 438,j),j=1,6)/  1, 1, 0, 3, 1, 0/
      data (ipa45( 439,j),j=1,6)/  1, 1, 1, 0, 0, 3/
      data (ipa45( 440,j),j=1,6)/  1, 1, 1, 0, 1, 2/
      data (ipa45( 441,j),j=1,6)/  1, 1, 1, 0, 2, 1/
      data (ipa45( 442,j),j=1,6)/  1, 1, 1, 0, 3, 0/
      data (ipa45( 443,j),j=1,6)/  1, 1, 1, 1, 0, 2/
      data (ipa45( 444,j),j=1,6)/  1, 1, 1, 1, 1, 1/
      data (ipa45( 445,j),j=1,6)/  1, 1, 1, 1, 2, 0/
      data (ipa45( 446,j),j=1,6)/  1, 1, 1, 2, 0, 1/
      data (ipa45( 447,j),j=1,6)/  1, 1, 1, 2, 1, 0/
      data (ipa45( 448,j),j=1,6)/  1, 1, 1, 3, 0, 0/
      data (ipa45( 449,j),j=1,6)/  1, 1, 2, 0, 0, 2/
      data (ipa45( 450,j),j=1,6)/  1, 1, 2, 0, 1, 1/
      data (ipa45( 451,j),j=1,6)/  1, 1, 2, 0, 2, 0/
      data (ipa45( 452,j),j=1,6)/  1, 1, 2, 1, 0, 1/
      data (ipa45( 453,j),j=1,6)/  1, 1, 2, 1, 1, 0/
      data (ipa45( 454,j),j=1,6)/  1, 1, 2, 2, 0, 0/
      data (ipa45( 455,j),j=1,6)/  1, 1, 3, 0, 0, 1/
      data (ipa45( 456,j),j=1,6)/  1, 1, 3, 0, 1, 0/
      data (ipa45( 457,j),j=1,6)/  1, 1, 3, 1, 0, 0/
      data (ipa45( 458,j),j=1,6)/  1, 1, 4, 0, 0, 0/
      data (ipa45( 459,j),j=1,6)/  1, 2, 0, 0, 0, 3/
      data (ipa45( 460,j),j=1,6)/  1, 2, 0, 0, 1, 2/
      data (ipa45( 461,j),j=1,6)/  1, 2, 0, 0, 2, 1/
      data (ipa45( 462,j),j=1,6)/  1, 2, 0, 0, 3, 0/
      data (ipa45( 463,j),j=1,6)/  1, 2, 0, 1, 0, 2/
      data (ipa45( 464,j),j=1,6)/  1, 2, 0, 1, 1, 1/
      data (ipa45( 465,j),j=1,6)/  1, 2, 0, 1, 2, 0/
      data (ipa45( 466,j),j=1,6)/  1, 2, 0, 2, 0, 1/
      data (ipa45( 467,j),j=1,6)/  1, 2, 0, 2, 1, 0/
      data (ipa45( 468,j),j=1,6)/  1, 2, 1, 0, 0, 2/
      data (ipa45( 469,j),j=1,6)/  1, 2, 1, 0, 1, 1/
      data (ipa45( 470,j),j=1,6)/  1, 2, 1, 0, 2, 0/
      data (ipa45( 471,j),j=1,6)/  1, 2, 1, 1, 0, 1/
      data (ipa45( 472,j),j=1,6)/  1, 2, 1, 1, 1, 0/
      data (ipa45( 473,j),j=1,6)/  1, 2, 1, 2, 0, 0/
      data (ipa45( 474,j),j=1,6)/  1, 2, 2, 0, 0, 1/
      data (ipa45( 475,j),j=1,6)/  1, 2, 2, 0, 1, 0/
      data (ipa45( 476,j),j=1,6)/  1, 2, 2, 1, 0, 0/
      data (ipa45( 477,j),j=1,6)/  1, 2, 3, 0, 0, 0/
      data (ipa45( 478,j),j=1,6)/  1, 3, 0, 0, 0, 2/
      data (ipa45( 479,j),j=1,6)/  1, 3, 0, 0, 1, 1/
      data (ipa45( 480,j),j=1,6)/  1, 3, 0, 0, 2, 0/
      data (ipa45( 481,j),j=1,6)/  1, 3, 0, 1, 0, 1/
      data (ipa45( 482,j),j=1,6)/  1, 3, 0, 1, 1, 0/
      data (ipa45( 483,j),j=1,6)/  1, 3, 1, 0, 0, 1/
      data (ipa45( 484,j),j=1,6)/  1, 3, 1, 0, 1, 0/
      data (ipa45( 485,j),j=1,6)/  1, 3, 1, 1, 0, 0/
      data (ipa45( 486,j),j=1,6)/  1, 3, 2, 0, 0, 0/
      data (ipa45( 487,j),j=1,6)/  1, 4, 0, 0, 0, 1/
      data (ipa45( 488,j),j=1,6)/  1, 4, 0, 0, 1, 0/
      data (ipa45( 489,j),j=1,6)/  1, 4, 1, 0, 0, 0/
      data (ipa45( 490,j),j=1,6)/  2, 0, 0, 0, 1, 3/
      data (ipa45( 491,j),j=1,6)/  2, 0, 0, 0, 2, 2/
      data (ipa45( 492,j),j=1,6)/  2, 0, 0, 0, 3, 1/
      data (ipa45( 493,j),j=1,6)/  2, 0, 0, 1, 0, 3/
      data (ipa45( 494,j),j=1,6)/  2, 0, 0, 1, 1, 2/
      data (ipa45( 495,j),j=1,6)/  2, 0, 0, 1, 2, 1/
      data (ipa45( 496,j),j=1,6)/  2, 0, 0, 1, 3, 0/
      data (ipa45( 497,j),j=1,6)/  2, 0, 0, 2, 0, 2/
      data (ipa45( 498,j),j=1,6)/  2, 0, 0, 2, 1, 1/
      data (ipa45( 499,j),j=1,6)/  2, 0, 0, 2, 2, 0/
      data (ipa45( 500,j),j=1,6)/  2, 0, 0, 3, 0, 1/
      data (ipa45( 501,j),j=1,6)/  2, 0, 0, 3, 1, 0/
      data (ipa45( 502,j),j=1,6)/  2, 0, 1, 0, 0, 3/
      data (ipa45( 503,j),j=1,6)/  2, 0, 1, 0, 1, 2/
      data (ipa45( 504,j),j=1,6)/  2, 0, 1, 0, 2, 1/
      data (ipa45( 505,j),j=1,6)/  2, 0, 1, 1, 0, 2/
      data (ipa45( 506,j),j=1,6)/  2, 0, 1, 1, 1, 1/
      data (ipa45( 507,j),j=1,6)/  2, 0, 1, 1, 2, 0/
      data (ipa45( 508,j),j=1,6)/  2, 0, 1, 2, 0, 1/
      data (ipa45( 509,j),j=1,6)/  2, 0, 1, 2, 1, 0/
      data (ipa45( 510,j),j=1,6)/  2, 0, 1, 3, 0, 0/
      data (ipa45( 511,j),j=1,6)/  2, 0, 2, 0, 0, 2/
      data (ipa45( 512,j),j=1,6)/  2, 0, 2, 0, 1, 1/
      data (ipa45( 513,j),j=1,6)/  2, 0, 2, 1, 0, 1/
      data (ipa45( 514,j),j=1,6)/  2, 0, 2, 1, 1, 0/
      data (ipa45( 515,j),j=1,6)/  2, 0, 2, 2, 0, 0/
      data (ipa45( 516,j),j=1,6)/  2, 0, 3, 0, 0, 1/
      data (ipa45( 517,j),j=1,6)/  2, 0, 3, 1, 0, 0/
      data (ipa45( 518,j),j=1,6)/  2, 1, 0, 0, 0, 3/
      data (ipa45( 519,j),j=1,6)/  2, 1, 0, 0, 1, 2/
      data (ipa45( 520,j),j=1,6)/  2, 1, 0, 0, 2, 1/
      data (ipa45( 521,j),j=1,6)/  2, 1, 0, 0, 3, 0/
      data (ipa45( 522,j),j=1,6)/  2, 1, 0, 1, 0, 2/
      data (ipa45( 523,j),j=1,6)/  2, 1, 0, 1, 1, 1/
      data (ipa45( 524,j),j=1,6)/  2, 1, 0, 1, 2, 0/
      data (ipa45( 525,j),j=1,6)/  2, 1, 0, 2, 0, 1/
      data (ipa45( 526,j),j=1,6)/  2, 1, 0, 2, 1, 0/
      data (ipa45( 527,j),j=1,6)/  2, 1, 1, 0, 0, 2/
      data (ipa45( 528,j),j=1,6)/  2, 1, 1, 0, 1, 1/
      data (ipa45( 529,j),j=1,6)/  2, 1, 1, 0, 2, 0/
      data (ipa45( 530,j),j=1,6)/  2, 1, 1, 1, 0, 1/
      data (ipa45( 531,j),j=1,6)/  2, 1, 1, 1, 1, 0/
      data (ipa45( 532,j),j=1,6)/  2, 1, 1, 2, 0, 0/
      data (ipa45( 533,j),j=1,6)/  2, 1, 2, 0, 0, 1/
      data (ipa45( 534,j),j=1,6)/  2, 1, 2, 0, 1, 0/
      data (ipa45( 535,j),j=1,6)/  2, 1, 2, 1, 0, 0/
      data (ipa45( 536,j),j=1,6)/  2, 1, 3, 0, 0, 0/
      data (ipa45( 537,j),j=1,6)/  2, 2, 0, 0, 0, 2/
      data (ipa45( 538,j),j=1,6)/  2, 2, 0, 0, 1, 1/
      data (ipa45( 539,j),j=1,6)/  2, 2, 0, 0, 2, 0/
      data (ipa45( 540,j),j=1,6)/  2, 2, 0, 1, 0, 1/
      data (ipa45( 541,j),j=1,6)/  2, 2, 0, 1, 1, 0/
      data (ipa45( 542,j),j=1,6)/  2, 2, 1, 0, 0, 1/
      data (ipa45( 543,j),j=1,6)/  2, 2, 1, 0, 1, 0/
      data (ipa45( 544,j),j=1,6)/  2, 2, 1, 1, 0, 0/
      data (ipa45( 545,j),j=1,6)/  2, 2, 2, 0, 0, 0/
      data (ipa45( 546,j),j=1,6)/  2, 3, 0, 0, 0, 1/
      data (ipa45( 547,j),j=1,6)/  2, 3, 0, 0, 1, 0/
      data (ipa45( 548,j),j=1,6)/  2, 3, 1, 0, 0, 0/
      data (ipa45( 549,j),j=1,6)/  3, 0, 0, 0, 1, 2/
      data (ipa45( 550,j),j=1,6)/  3, 0, 0, 0, 2, 1/
      data (ipa45( 551,j),j=1,6)/  3, 0, 0, 1, 0, 2/
      data (ipa45( 552,j),j=1,6)/  3, 0, 0, 1, 1, 1/
      data (ipa45( 553,j),j=1,6)/  3, 0, 0, 1, 2, 0/
      data (ipa45( 554,j),j=1,6)/  3, 0, 0, 2, 0, 1/
      data (ipa45( 555,j),j=1,6)/  3, 0, 0, 2, 1, 0/
      data (ipa45( 556,j),j=1,6)/  3, 0, 1, 0, 0, 2/
      data (ipa45( 557,j),j=1,6)/  3, 0, 1, 0, 1, 1/
      data (ipa45( 558,j),j=1,6)/  3, 0, 1, 1, 0, 1/
      data (ipa45( 559,j),j=1,6)/  3, 0, 1, 1, 1, 0/
      data (ipa45( 560,j),j=1,6)/  3, 0, 1, 2, 0, 0/
      data (ipa45( 561,j),j=1,6)/  3, 0, 2, 0, 0, 1/
      data (ipa45( 562,j),j=1,6)/  3, 0, 2, 1, 0, 0/
      data (ipa45( 563,j),j=1,6)/  3, 1, 0, 0, 0, 2/
      data (ipa45( 564,j),j=1,6)/  3, 1, 0, 0, 1, 1/
      data (ipa45( 565,j),j=1,6)/  3, 1, 0, 0, 2, 0/
      data (ipa45( 566,j),j=1,6)/  3, 1, 0, 1, 0, 1/
      data (ipa45( 567,j),j=1,6)/  3, 1, 0, 1, 1, 0/
      data (ipa45( 568,j),j=1,6)/  3, 1, 1, 0, 0, 1/
      data (ipa45( 569,j),j=1,6)/  3, 1, 1, 0, 1, 0/
      data (ipa45( 570,j),j=1,6)/  3, 1, 1, 1, 0, 0/
      data (ipa45( 571,j),j=1,6)/  3, 1, 2, 0, 0, 0/
      data (ipa45( 572,j),j=1,6)/  3, 2, 0, 0, 0, 1/
      data (ipa45( 573,j),j=1,6)/  3, 2, 0, 0, 1, 0/
      data (ipa45( 574,j),j=1,6)/  3, 2, 1, 0, 0, 0/
      data (ipa45( 575,j),j=1,6)/  4, 0, 0, 0, 1, 1/
      data (ipa45( 576,j),j=1,6)/  4, 0, 0, 1, 0, 1/
      data (ipa45( 577,j),j=1,6)/  4, 0, 0, 1, 1, 0/
      data (ipa45( 578,j),j=1,6)/  4, 0, 1, 0, 0, 1/
      data (ipa45( 579,j),j=1,6)/  4, 0, 1, 1, 0, 0/
      data (ipa45( 580,j),j=1,6)/  4, 1, 0, 0, 0, 1/
      data (ipa45( 581,j),j=1,6)/  4, 1, 0, 0, 1, 0/
      data (ipa45( 582,j),j=1,6)/  4, 1, 1, 0, 0, 0/
*********** 6 341*******************************
      data (nfas45(i),i=1,6)/0,0,16,79,241,582/
      end
