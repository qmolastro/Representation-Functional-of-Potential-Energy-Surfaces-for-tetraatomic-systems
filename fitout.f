      subroutine fit6d(r12,r13,r14,r23,r24,r34,ener,der,iop)
      implicit  real * 8 (a-h,o-z)
      dimension d123(3),d124(3),d134(3),d234(3),d1234(6),der(6)
      call diat12(r12,e12,d12,iop)
      call diat13(r13,e13,d13,iop)
      call diat14(r14,e14,d14,iop)
      call diat23(r23,e23,d23,iop)
      call diat24(r24,e24,d24,iop)
      call diat34(r34,e34,d34,iop)
      call tri123(r12,r13,r23,e123,d123,iop)
      call tri124(r12,r14,r24,e124,d124,iop)
      call tri134(r13,r14,r34,e134,d134,iop)
      call tri234(r23,r24,r34,e234,d234,iop)
      call abcd(r12,r13,r14,r23,r24,r34,e1234,d1234,iop)
      ener=e12+e13+e14+e23+e24+e34+e123+e124+e134+e234+e1234
      der(1)=d12+d123(1)+d124(1)+d1234(1) 
      der(2)=d13+d123(2)+d134(1)+d1234(2) 
      der(3)=d14+d124(2)+d134(2)+d1234(3) 
      der(4)=d23+d123(3)+d234(1)+d1234(4) 
      der(5)=d24+d124(3)+d234(2)+d1234(5) 
      der(6)=d34+d134(3)+d234(3)+d1234(6) 
      end
************************************************************************
      subroutine diat12(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
************************************************************************
      subroutine diat13(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
************************************************************************
      subroutine diat14(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
************************************************************************
      subroutine diat23(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
************************************************************************
      subroutine diat24(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
************************************************************************
      subroutine diat34(r,ener,der,iop)
************************************************************************
*     This subroutine computes the energies of a diatomic potential 
*     fitted to    87 points
*     rms =      0.12232259 kcal/mol
*     emax =      0.22786545 kcal/mol
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(mt=  6,e0= 0.0000000D+00)
      parameter(vex1=    0.165794595000D+01,vex2=    0.192681240000D+01)
      dimension cf(mt)
      data cf(  1)/    0.100582864185D+01/
      data cf(  2)/   -0.328787241355D+01/
      data cf(  3)/    0.270836301221D+02/
      data cf(  4)/   -0.180246256518D+03/
      data cf(  5)/    0.647540308481D+03/
      data cf(  6)/   -0.933018533257D+03/
      aux = 1.d0/r
      bux = dexp(-vex2*r)*aux
      cux = dexp(-vex1*r)
      ener=e0+cf(1)*bux
      dux=1.d0
      eux=r*cux
      do 1 i=2,mt
         if (iop.eq.1) der=der+(i-1)*cf(i)*dux
         dux=dux*eux
         ener=ener+cf(i)*dux
    1 continue
      if (iop.eq.1) then
         der=der*(1.d0-vex1*r)*cux
         der=der-cf(1)*(vex2+aux)*bux
      endif
      return
      end
*************************************************************
      subroutine tri123(r12,r13,r23,ener,der,iop)          
*************************************************************
*     This subroutine computes the energies of a 3D PES      
*     for the ABC system class fitted to  772 points       
*     rms =      0.56497363 kcal/mol                               
*     emax=      5.20498557 kcal/mol                               
*************************************************************
      implicit real*8(a-h,o-z)                               
      parameter(is=   40, mt_1= 4)                      
      parameter(vex1=    0.156190000000D+00)                             
      parameter(vex2=    0.150000100000D+01)                             
      parameter(vex3=    0.124744700000D+01)                             
      dimension der(3)                                       
      dimension i1(is),i2(is),i3(is),cf(is)                  
      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       
      data der12/0.d0/,der13/0.d0/,der23/0.d0/               
      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         
      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            
      data cf(  1)/    0.379652668089D+02/
      data i1(  1)/ 0/,i2(  1)/ 1/,i3(  1)/ 1/
      data cf(  2)/    0.109161179953D+01/
      data i1(  2)/ 1/,i2(  2)/ 0/,i3(  2)/ 1/
      data cf(  3)/    0.647458872025D+01/
      data i1(  3)/ 1/,i2(  3)/ 1/,i3(  3)/ 0/
      data cf(  4)/   -0.122207741460D+03/
      data i1(  4)/ 0/,i2(  4)/ 1/,i3(  4)/ 2/
      data cf(  5)/   -0.126473537977D+03/
      data i1(  5)/ 0/,i2(  5)/ 2/,i3(  5)/ 1/
      data cf(  6)/    0.277415510727D+02/
      data i1(  6)/ 1/,i2(  6)/ 0/,i3(  6)/ 2/
      data cf(  7)/   -0.828077993379D+02/
      data i1(  7)/ 1/,i2(  7)/ 1/,i3(  7)/ 1/
      data cf(  8)/   -0.365242295454D+01/
      data i1(  8)/ 1/,i2(  8)/ 2/,i3(  8)/ 0/
      data cf(  9)/   -0.876852493490D+00/
      data i1(  9)/ 2/,i2(  9)/ 0/,i3(  9)/ 1/
      data cf( 10)/   -0.835290191933D+01/
      data i1( 10)/ 2/,i2( 10)/ 1/,i3( 10)/ 0/
      data cf( 11)/    0.694723171102D+03/
      data i1( 11)/ 0/,i2( 11)/ 1/,i3( 11)/ 3/
      data cf( 12)/   -0.802353835651D+03/
      data i1( 12)/ 0/,i2( 12)/ 2/,i3( 12)/ 2/
      data cf( 13)/    0.959983210877D+03/
      data i1( 13)/ 0/,i2( 13)/ 3/,i3( 13)/ 1/
      data cf( 14)/   -0.127052191213D+03/
      data i1( 14)/ 1/,i2( 14)/ 0/,i3( 14)/ 3/
      data cf( 15)/    0.158644333845D+03/
      data i1( 15)/ 1/,i2( 15)/ 1/,i3( 15)/ 2/
      data cf( 16)/    0.209256144232D+03/
      data i1( 16)/ 1/,i2( 16)/ 2/,i3( 16)/ 1/
      data cf( 17)/   -0.112445084995D+03/
      data i1( 17)/ 1/,i2( 17)/ 3/,i3( 17)/ 0/
      data cf( 18)/   -0.216564410045D+02/
      data i1( 18)/ 2/,i2( 18)/ 0/,i3( 18)/ 2/
      data cf( 19)/    0.615761512277D+02/
      data i1( 19)/ 2/,i2( 19)/ 1/,i3( 19)/ 1/
      data cf( 20)/    0.103847737925D+02/
      data i1( 20)/ 2/,i2( 20)/ 2/,i3( 20)/ 0/
      data cf( 21)/    0.223268094402D+00/
      data i1( 21)/ 3/,i2( 21)/ 0/,i3( 21)/ 1/
      data cf( 22)/    0.464511223351D+01/
      data i1( 22)/ 3/,i2( 22)/ 1/,i3( 22)/ 0/
      data cf( 23)/   -0.258801534845D+03/
      data i1( 23)/ 0/,i2( 23)/ 1/,i3( 23)/ 4/
      data cf( 24)/   -0.169242354541D+04/
      data i1( 24)/ 0/,i2( 24)/ 2/,i3( 24)/ 3/
      data cf( 25)/    0.306755920001D+04/
      data i1( 25)/ 0/,i2( 25)/ 3/,i3( 25)/ 2/
      data cf( 26)/   -0.237466134143D+04/
      data i1( 26)/ 0/,i2( 26)/ 4/,i3( 26)/ 1/
      data cf( 27)/    0.324207797105D+02/
      data i1( 27)/ 1/,i2( 27)/ 0/,i3( 27)/ 4/
      data cf( 28)/    0.151708332767D+03/
      data i1( 28)/ 1/,i2( 28)/ 1/,i3( 28)/ 3/
      data cf( 29)/   -0.177820056390D+03/
      data i1( 29)/ 1/,i2( 29)/ 2/,i3( 29)/ 2/
      data cf( 30)/   -0.304686678991D+03/
      data i1( 30)/ 1/,i2( 30)/ 3/,i3( 30)/ 1/
      data cf( 31)/    0.228158098368D+03/
      data i1( 31)/ 1/,i2( 31)/ 4/,i3( 31)/ 0/
      data cf( 32)/    0.411570815014D+02/
      data i1( 32)/ 2/,i2( 32)/ 0/,i3( 32)/ 3/
      data cf( 33)/   -0.454749485386D+02/
      data i1( 33)/ 2/,i2( 33)/ 1/,i3( 33)/ 2/
      data cf( 34)/   -0.981048969167D+02/
      data i1( 34)/ 2/,i2( 34)/ 2/,i3( 34)/ 1/
      data cf( 35)/    0.423301056810D+02/
      data i1( 35)/ 2/,i2( 35)/ 3/,i3( 35)/ 0/
      data cf( 36)/    0.554900277613D+01/
      data i1( 36)/ 3/,i2( 36)/ 0/,i3( 36)/ 2/
      data cf( 37)/   -0.134271735248D+02/
      data i1( 37)/ 3/,i2( 37)/ 1/,i3( 37)/ 1/
      data cf( 38)/   -0.639347965902D+01/
      data i1( 38)/ 3/,i2( 38)/ 2/,i3( 38)/ 0/
      data cf( 39)/    0.157063047025D-02/
      data i1( 39)/ 4/,i2( 39)/ 0/,i3( 39)/ 1/
      data cf( 40)/   -0.946258930769D+00/
      data i1( 40)/ 4/,i2( 40)/ 1/,i3( 40)/ 0/
      ener=0.d0
      pux12=vex1*r12
      pux13=vex2*r13
      pux23=vex3*r23
      qux12=dexp(-pux12)
      qux13=dexp(-pux13)
      qux23=dexp(-pux23)
      bux12=r12*qux12
      bux13=r13*qux13
      bux23=r23*qux23
      do 1 i=1,mt_1
         f12(i)=f12(i-1)*bux12
         f13(i)=f13(i-1)*bux13
         f23(i)=f23(i-1)*bux23
1     continue

      do 2 l=1,is
      aux=f12(i1(l))*f13(i2(l))*f23(i3(l))           
      ener=ener+cf(l)*aux                         
      if (iop.eq.1) then                          
         dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))   
         dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))   
         dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)   
         der12=der12+cf(l)*dux12                          
         der13=der13+cf(l)*dux13                          
         der23=der23+cf(l)*dux23                          
      endif                         
    2 continue                                     
      if (iop.eq.1) then                          
         der(1)=der12*(1.d0-pux12)*qux12        
         der(2)=der13*(1.d0-pux13)*qux13        
         der(3)=der23*(1.d0-pux23)*qux23        
      endif                         
      return
      end
*************************************************************
      subroutine tri124(r12,r13,r23,ener,der,iop)          
*************************************************************
*     This subroutine computes the energies of a 3D PES      
*     for the ABC system class fitted to  772 points       
*     rms =      0.56497363 kcal/mol                               
*     emax=      5.20498557 kcal/mol                               
*************************************************************
      implicit real*8(a-h,o-z)                               
      parameter(is=   40, mt_1= 4)                      
      parameter(vex1=    0.156190000000D+00)                             
      parameter(vex2=    0.150000100000D+01)                             
      parameter(vex3=    0.124744700000D+01)                             
      dimension der(3)                                       
      dimension i1(is),i2(is),i3(is),cf(is)                  
      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       
      data der12/0.d0/,der13/0.d0/,der23/0.d0/               
      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         
      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            
      data cf(  1)/    0.379652668089D+02/
      data i1(  1)/ 0/,i2(  1)/ 1/,i3(  1)/ 1/
      data cf(  2)/    0.109161179953D+01/
      data i1(  2)/ 1/,i2(  2)/ 0/,i3(  2)/ 1/
      data cf(  3)/    0.647458872025D+01/
      data i1(  3)/ 1/,i2(  3)/ 1/,i3(  3)/ 0/
      data cf(  4)/   -0.122207741460D+03/
      data i1(  4)/ 0/,i2(  4)/ 1/,i3(  4)/ 2/
      data cf(  5)/   -0.126473537977D+03/
      data i1(  5)/ 0/,i2(  5)/ 2/,i3(  5)/ 1/
      data cf(  6)/    0.277415510727D+02/
      data i1(  6)/ 1/,i2(  6)/ 0/,i3(  6)/ 2/
      data cf(  7)/   -0.828077993379D+02/
      data i1(  7)/ 1/,i2(  7)/ 1/,i3(  7)/ 1/
      data cf(  8)/   -0.365242295454D+01/
      data i1(  8)/ 1/,i2(  8)/ 2/,i3(  8)/ 0/
      data cf(  9)/   -0.876852493490D+00/
      data i1(  9)/ 2/,i2(  9)/ 0/,i3(  9)/ 1/
      data cf( 10)/   -0.835290191933D+01/
      data i1( 10)/ 2/,i2( 10)/ 1/,i3( 10)/ 0/
      data cf( 11)/    0.694723171102D+03/
      data i1( 11)/ 0/,i2( 11)/ 1/,i3( 11)/ 3/
      data cf( 12)/   -0.802353835651D+03/
      data i1( 12)/ 0/,i2( 12)/ 2/,i3( 12)/ 2/
      data cf( 13)/    0.959983210877D+03/
      data i1( 13)/ 0/,i2( 13)/ 3/,i3( 13)/ 1/
      data cf( 14)/   -0.127052191213D+03/
      data i1( 14)/ 1/,i2( 14)/ 0/,i3( 14)/ 3/
      data cf( 15)/    0.158644333845D+03/
      data i1( 15)/ 1/,i2( 15)/ 1/,i3( 15)/ 2/
      data cf( 16)/    0.209256144232D+03/
      data i1( 16)/ 1/,i2( 16)/ 2/,i3( 16)/ 1/
      data cf( 17)/   -0.112445084995D+03/
      data i1( 17)/ 1/,i2( 17)/ 3/,i3( 17)/ 0/
      data cf( 18)/   -0.216564410045D+02/
      data i1( 18)/ 2/,i2( 18)/ 0/,i3( 18)/ 2/
      data cf( 19)/    0.615761512277D+02/
      data i1( 19)/ 2/,i2( 19)/ 1/,i3( 19)/ 1/
      data cf( 20)/    0.103847737925D+02/
      data i1( 20)/ 2/,i2( 20)/ 2/,i3( 20)/ 0/
      data cf( 21)/    0.223268094402D+00/
      data i1( 21)/ 3/,i2( 21)/ 0/,i3( 21)/ 1/
      data cf( 22)/    0.464511223351D+01/
      data i1( 22)/ 3/,i2( 22)/ 1/,i3( 22)/ 0/
      data cf( 23)/   -0.258801534845D+03/
      data i1( 23)/ 0/,i2( 23)/ 1/,i3( 23)/ 4/
      data cf( 24)/   -0.169242354541D+04/
      data i1( 24)/ 0/,i2( 24)/ 2/,i3( 24)/ 3/
      data cf( 25)/    0.306755920001D+04/
      data i1( 25)/ 0/,i2( 25)/ 3/,i3( 25)/ 2/
      data cf( 26)/   -0.237466134143D+04/
      data i1( 26)/ 0/,i2( 26)/ 4/,i3( 26)/ 1/
      data cf( 27)/    0.324207797105D+02/
      data i1( 27)/ 1/,i2( 27)/ 0/,i3( 27)/ 4/
      data cf( 28)/    0.151708332767D+03/
      data i1( 28)/ 1/,i2( 28)/ 1/,i3( 28)/ 3/
      data cf( 29)/   -0.177820056390D+03/
      data i1( 29)/ 1/,i2( 29)/ 2/,i3( 29)/ 2/
      data cf( 30)/   -0.304686678991D+03/
      data i1( 30)/ 1/,i2( 30)/ 3/,i3( 30)/ 1/
      data cf( 31)/    0.228158098368D+03/
      data i1( 31)/ 1/,i2( 31)/ 4/,i3( 31)/ 0/
      data cf( 32)/    0.411570815014D+02/
      data i1( 32)/ 2/,i2( 32)/ 0/,i3( 32)/ 3/
      data cf( 33)/   -0.454749485386D+02/
      data i1( 33)/ 2/,i2( 33)/ 1/,i3( 33)/ 2/
      data cf( 34)/   -0.981048969167D+02/
      data i1( 34)/ 2/,i2( 34)/ 2/,i3( 34)/ 1/
      data cf( 35)/    0.423301056810D+02/
      data i1( 35)/ 2/,i2( 35)/ 3/,i3( 35)/ 0/
      data cf( 36)/    0.554900277613D+01/
      data i1( 36)/ 3/,i2( 36)/ 0/,i3( 36)/ 2/
      data cf( 37)/   -0.134271735248D+02/
      data i1( 37)/ 3/,i2( 37)/ 1/,i3( 37)/ 1/
      data cf( 38)/   -0.639347965902D+01/
      data i1( 38)/ 3/,i2( 38)/ 2/,i3( 38)/ 0/
      data cf( 39)/    0.157063047025D-02/
      data i1( 39)/ 4/,i2( 39)/ 0/,i3( 39)/ 1/
      data cf( 40)/   -0.946258930769D+00/
      data i1( 40)/ 4/,i2( 40)/ 1/,i3( 40)/ 0/
      ener=0.d0
      pux12=vex1*r12
      pux13=vex2*r13
      pux23=vex3*r23
      qux12=dexp(-pux12)
      qux13=dexp(-pux13)
      qux23=dexp(-pux23)
      bux12=r12*qux12
      bux13=r13*qux13
      bux23=r23*qux23
      do 1 i=1,mt_1
         f12(i)=f12(i-1)*bux12
         f13(i)=f13(i-1)*bux13
         f23(i)=f23(i-1)*bux23
1     continue

      do 2 l=1,is
      aux=f12(i1(l))*f13(i2(l))*f23(i3(l))           
      ener=ener+cf(l)*aux                         
      if (iop.eq.1) then                          
         dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))   
         dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))   
         dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)   
         der12=der12+cf(l)*dux12                          
         der13=der13+cf(l)*dux13                          
         der23=der23+cf(l)*dux23                          
      endif                         
    2 continue                                     
      if (iop.eq.1) then                          
         der(1)=der12*(1.d0-pux12)*qux12        
         der(2)=der13*(1.d0-pux13)*qux13        
         der(3)=der23*(1.d0-pux23)*qux23        
      endif                         
      return
      end
*************************************************************
      subroutine tri134(r12,r13,r23,ener,der,iop)          
*************************************************************
*     This subroutine computes the energies of a 3D PES      
*     for the ABC system class fitted to  772 points       
*     rms =      0.56497363 kcal/mol                               
*     emax=      5.20498557 kcal/mol                               
*************************************************************
      implicit real*8(a-h,o-z)                               
      parameter(is=   40, mt_1= 4)                      
      parameter(vex1=    0.156190000000D+00)                             
      parameter(vex2=    0.150000100000D+01)                             
      parameter(vex3=    0.124744700000D+01)                             
      dimension der(3)                                       
      dimension i1(is),i2(is),i3(is),cf(is)                  
      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       
      data der12/0.d0/,der13/0.d0/,der23/0.d0/               
      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         
      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            
      data cf(  1)/    0.379652668089D+02/
      data i1(  1)/ 0/,i2(  1)/ 1/,i3(  1)/ 1/
      data cf(  2)/    0.109161179953D+01/
      data i1(  2)/ 1/,i2(  2)/ 0/,i3(  2)/ 1/
      data cf(  3)/    0.647458872025D+01/
      data i1(  3)/ 1/,i2(  3)/ 1/,i3(  3)/ 0/
      data cf(  4)/   -0.122207741460D+03/
      data i1(  4)/ 0/,i2(  4)/ 1/,i3(  4)/ 2/
      data cf(  5)/   -0.126473537977D+03/
      data i1(  5)/ 0/,i2(  5)/ 2/,i3(  5)/ 1/
      data cf(  6)/    0.277415510727D+02/
      data i1(  6)/ 1/,i2(  6)/ 0/,i3(  6)/ 2/
      data cf(  7)/   -0.828077993379D+02/
      data i1(  7)/ 1/,i2(  7)/ 1/,i3(  7)/ 1/
      data cf(  8)/   -0.365242295454D+01/
      data i1(  8)/ 1/,i2(  8)/ 2/,i3(  8)/ 0/
      data cf(  9)/   -0.876852493490D+00/
      data i1(  9)/ 2/,i2(  9)/ 0/,i3(  9)/ 1/
      data cf( 10)/   -0.835290191933D+01/
      data i1( 10)/ 2/,i2( 10)/ 1/,i3( 10)/ 0/
      data cf( 11)/    0.694723171102D+03/
      data i1( 11)/ 0/,i2( 11)/ 1/,i3( 11)/ 3/
      data cf( 12)/   -0.802353835651D+03/
      data i1( 12)/ 0/,i2( 12)/ 2/,i3( 12)/ 2/
      data cf( 13)/    0.959983210877D+03/
      data i1( 13)/ 0/,i2( 13)/ 3/,i3( 13)/ 1/
      data cf( 14)/   -0.127052191213D+03/
      data i1( 14)/ 1/,i2( 14)/ 0/,i3( 14)/ 3/
      data cf( 15)/    0.158644333845D+03/
      data i1( 15)/ 1/,i2( 15)/ 1/,i3( 15)/ 2/
      data cf( 16)/    0.209256144232D+03/
      data i1( 16)/ 1/,i2( 16)/ 2/,i3( 16)/ 1/
      data cf( 17)/   -0.112445084995D+03/
      data i1( 17)/ 1/,i2( 17)/ 3/,i3( 17)/ 0/
      data cf( 18)/   -0.216564410045D+02/
      data i1( 18)/ 2/,i2( 18)/ 0/,i3( 18)/ 2/
      data cf( 19)/    0.615761512277D+02/
      data i1( 19)/ 2/,i2( 19)/ 1/,i3( 19)/ 1/
      data cf( 20)/    0.103847737925D+02/
      data i1( 20)/ 2/,i2( 20)/ 2/,i3( 20)/ 0/
      data cf( 21)/    0.223268094402D+00/
      data i1( 21)/ 3/,i2( 21)/ 0/,i3( 21)/ 1/
      data cf( 22)/    0.464511223351D+01/
      data i1( 22)/ 3/,i2( 22)/ 1/,i3( 22)/ 0/
      data cf( 23)/   -0.258801534845D+03/
      data i1( 23)/ 0/,i2( 23)/ 1/,i3( 23)/ 4/
      data cf( 24)/   -0.169242354541D+04/
      data i1( 24)/ 0/,i2( 24)/ 2/,i3( 24)/ 3/
      data cf( 25)/    0.306755920001D+04/
      data i1( 25)/ 0/,i2( 25)/ 3/,i3( 25)/ 2/
      data cf( 26)/   -0.237466134143D+04/
      data i1( 26)/ 0/,i2( 26)/ 4/,i3( 26)/ 1/
      data cf( 27)/    0.324207797105D+02/
      data i1( 27)/ 1/,i2( 27)/ 0/,i3( 27)/ 4/
      data cf( 28)/    0.151708332767D+03/
      data i1( 28)/ 1/,i2( 28)/ 1/,i3( 28)/ 3/
      data cf( 29)/   -0.177820056390D+03/
      data i1( 29)/ 1/,i2( 29)/ 2/,i3( 29)/ 2/
      data cf( 30)/   -0.304686678991D+03/
      data i1( 30)/ 1/,i2( 30)/ 3/,i3( 30)/ 1/
      data cf( 31)/    0.228158098368D+03/
      data i1( 31)/ 1/,i2( 31)/ 4/,i3( 31)/ 0/
      data cf( 32)/    0.411570815014D+02/
      data i1( 32)/ 2/,i2( 32)/ 0/,i3( 32)/ 3/
      data cf( 33)/   -0.454749485386D+02/
      data i1( 33)/ 2/,i2( 33)/ 1/,i3( 33)/ 2/
      data cf( 34)/   -0.981048969167D+02/
      data i1( 34)/ 2/,i2( 34)/ 2/,i3( 34)/ 1/
      data cf( 35)/    0.423301056810D+02/
      data i1( 35)/ 2/,i2( 35)/ 3/,i3( 35)/ 0/
      data cf( 36)/    0.554900277613D+01/
      data i1( 36)/ 3/,i2( 36)/ 0/,i3( 36)/ 2/
      data cf( 37)/   -0.134271735248D+02/
      data i1( 37)/ 3/,i2( 37)/ 1/,i3( 37)/ 1/
      data cf( 38)/   -0.639347965902D+01/
      data i1( 38)/ 3/,i2( 38)/ 2/,i3( 38)/ 0/
      data cf( 39)/    0.157063047025D-02/
      data i1( 39)/ 4/,i2( 39)/ 0/,i3( 39)/ 1/
      data cf( 40)/   -0.946258930769D+00/
      data i1( 40)/ 4/,i2( 40)/ 1/,i3( 40)/ 0/
      ener=0.d0
      pux12=vex1*r12
      pux13=vex2*r13
      pux23=vex3*r23
      qux12=dexp(-pux12)
      qux13=dexp(-pux13)
      qux23=dexp(-pux23)
      bux12=r12*qux12
      bux13=r13*qux13
      bux23=r23*qux23
      do 1 i=1,mt_1
         f12(i)=f12(i-1)*bux12
         f13(i)=f13(i-1)*bux13
         f23(i)=f23(i-1)*bux23
1     continue

      do 2 l=1,is
      aux=f12(i1(l))*f13(i2(l))*f23(i3(l))           
      ener=ener+cf(l)*aux                         
      if (iop.eq.1) then                          
         dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))   
         dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))   
         dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)   
         der12=der12+cf(l)*dux12                          
         der13=der13+cf(l)*dux13                          
         der23=der23+cf(l)*dux23                          
      endif                         
    2 continue                                     
      if (iop.eq.1) then                          
         der(1)=der12*(1.d0-pux12)*qux12        
         der(2)=der13*(1.d0-pux13)*qux13        
         der(3)=der23*(1.d0-pux23)*qux23        
      endif                         
      return
      end
*************************************************************
      subroutine tri234(r12,r13,r23,ener,der,iop)          
*************************************************************
*     This subroutine computes the energies of a 3D PES      
*     for the ABC system class fitted to  772 points       
*     rms =      0.56497363 kcal/mol                               
*     emax=      5.20498557 kcal/mol                               
*************************************************************
      implicit real*8(a-h,o-z)                               
      parameter(is=   40, mt_1= 4)                      
      parameter(vex1=    0.156190000000D+00)                             
      parameter(vex2=    0.150000100000D+01)                             
      parameter(vex3=    0.124744700000D+01)                             
      dimension der(3)                                       
      dimension i1(is),i2(is),i3(is),cf(is)                  
      dimension f12(-1:mt_1),f13(-1:mt_1),f23(-1:mt_1)       
      data der12/0.d0/,der13/0.d0/,der23/0.d0/               
      data f12(-1)/0.d0/,f13(-1)/0.d0/,f23(-1)/0.d0/         
      data f12(0)/1.d0/,f13(0)/1.d0/,f23(0)/1.d0/            
      data cf(  1)/    0.379652668089D+02/
      data i1(  1)/ 0/,i2(  1)/ 1/,i3(  1)/ 1/
      data cf(  2)/    0.109161179953D+01/
      data i1(  2)/ 1/,i2(  2)/ 0/,i3(  2)/ 1/
      data cf(  3)/    0.647458872025D+01/
      data i1(  3)/ 1/,i2(  3)/ 1/,i3(  3)/ 0/
      data cf(  4)/   -0.122207741460D+03/
      data i1(  4)/ 0/,i2(  4)/ 1/,i3(  4)/ 2/
      data cf(  5)/   -0.126473537977D+03/
      data i1(  5)/ 0/,i2(  5)/ 2/,i3(  5)/ 1/
      data cf(  6)/    0.277415510727D+02/
      data i1(  6)/ 1/,i2(  6)/ 0/,i3(  6)/ 2/
      data cf(  7)/   -0.828077993379D+02/
      data i1(  7)/ 1/,i2(  7)/ 1/,i3(  7)/ 1/
      data cf(  8)/   -0.365242295454D+01/
      data i1(  8)/ 1/,i2(  8)/ 2/,i3(  8)/ 0/
      data cf(  9)/   -0.876852493490D+00/
      data i1(  9)/ 2/,i2(  9)/ 0/,i3(  9)/ 1/
      data cf( 10)/   -0.835290191933D+01/
      data i1( 10)/ 2/,i2( 10)/ 1/,i3( 10)/ 0/
      data cf( 11)/    0.694723171102D+03/
      data i1( 11)/ 0/,i2( 11)/ 1/,i3( 11)/ 3/
      data cf( 12)/   -0.802353835651D+03/
      data i1( 12)/ 0/,i2( 12)/ 2/,i3( 12)/ 2/
      data cf( 13)/    0.959983210877D+03/
      data i1( 13)/ 0/,i2( 13)/ 3/,i3( 13)/ 1/
      data cf( 14)/   -0.127052191213D+03/
      data i1( 14)/ 1/,i2( 14)/ 0/,i3( 14)/ 3/
      data cf( 15)/    0.158644333845D+03/
      data i1( 15)/ 1/,i2( 15)/ 1/,i3( 15)/ 2/
      data cf( 16)/    0.209256144232D+03/
      data i1( 16)/ 1/,i2( 16)/ 2/,i3( 16)/ 1/
      data cf( 17)/   -0.112445084995D+03/
      data i1( 17)/ 1/,i2( 17)/ 3/,i3( 17)/ 0/
      data cf( 18)/   -0.216564410045D+02/
      data i1( 18)/ 2/,i2( 18)/ 0/,i3( 18)/ 2/
      data cf( 19)/    0.615761512277D+02/
      data i1( 19)/ 2/,i2( 19)/ 1/,i3( 19)/ 1/
      data cf( 20)/    0.103847737925D+02/
      data i1( 20)/ 2/,i2( 20)/ 2/,i3( 20)/ 0/
      data cf( 21)/    0.223268094402D+00/
      data i1( 21)/ 3/,i2( 21)/ 0/,i3( 21)/ 1/
      data cf( 22)/    0.464511223351D+01/
      data i1( 22)/ 3/,i2( 22)/ 1/,i3( 22)/ 0/
      data cf( 23)/   -0.258801534845D+03/
      data i1( 23)/ 0/,i2( 23)/ 1/,i3( 23)/ 4/
      data cf( 24)/   -0.169242354541D+04/
      data i1( 24)/ 0/,i2( 24)/ 2/,i3( 24)/ 3/
      data cf( 25)/    0.306755920001D+04/
      data i1( 25)/ 0/,i2( 25)/ 3/,i3( 25)/ 2/
      data cf( 26)/   -0.237466134143D+04/
      data i1( 26)/ 0/,i2( 26)/ 4/,i3( 26)/ 1/
      data cf( 27)/    0.324207797105D+02/
      data i1( 27)/ 1/,i2( 27)/ 0/,i3( 27)/ 4/
      data cf( 28)/    0.151708332767D+03/
      data i1( 28)/ 1/,i2( 28)/ 1/,i3( 28)/ 3/
      data cf( 29)/   -0.177820056390D+03/
      data i1( 29)/ 1/,i2( 29)/ 2/,i3( 29)/ 2/
      data cf( 30)/   -0.304686678991D+03/
      data i1( 30)/ 1/,i2( 30)/ 3/,i3( 30)/ 1/
      data cf( 31)/    0.228158098368D+03/
      data i1( 31)/ 1/,i2( 31)/ 4/,i3( 31)/ 0/
      data cf( 32)/    0.411570815014D+02/
      data i1( 32)/ 2/,i2( 32)/ 0/,i3( 32)/ 3/
      data cf( 33)/   -0.454749485386D+02/
      data i1( 33)/ 2/,i2( 33)/ 1/,i3( 33)/ 2/
      data cf( 34)/   -0.981048969167D+02/
      data i1( 34)/ 2/,i2( 34)/ 2/,i3( 34)/ 1/
      data cf( 35)/    0.423301056810D+02/
      data i1( 35)/ 2/,i2( 35)/ 3/,i3( 35)/ 0/
      data cf( 36)/    0.554900277613D+01/
      data i1( 36)/ 3/,i2( 36)/ 0/,i3( 36)/ 2/
      data cf( 37)/   -0.134271735248D+02/
      data i1( 37)/ 3/,i2( 37)/ 1/,i3( 37)/ 1/
      data cf( 38)/   -0.639347965902D+01/
      data i1( 38)/ 3/,i2( 38)/ 2/,i3( 38)/ 0/
      data cf( 39)/    0.157063047025D-02/
      data i1( 39)/ 4/,i2( 39)/ 0/,i3( 39)/ 1/
      data cf( 40)/   -0.946258930769D+00/
      data i1( 40)/ 4/,i2( 40)/ 1/,i3( 40)/ 0/
      ener=0.d0
      pux12=vex1*r12
      pux13=vex2*r13
      pux23=vex3*r23
      qux12=dexp(-pux12)
      qux13=dexp(-pux13)
      qux23=dexp(-pux23)
      bux12=r12*qux12
      bux13=r13*qux13
      bux23=r23*qux23
      do 1 i=1,mt_1
         f12(i)=f12(i-1)*bux12
         f13(i)=f13(i-1)*bux13
         f23(i)=f23(i-1)*bux23
1     continue

      do 2 l=1,is
      aux=f12(i1(l))*f13(i2(l))*f23(i3(l))           
      ener=ener+cf(l)*aux                         
      if (iop.eq.1) then                          
         dux12=i1(l)*f12(i1(l)-1)*f13(i2(l))*f23(i3(l))   
         dux13=i2(l)*f12(i1(l))*f13(i2(l)-1)*f23(i3(l))   
         dux23=i3(l)*f12(i1(l))*f13(i2(l))*f23(i3(l)-1)   
         der12=der12+cf(l)*dux12                          
         der13=der13+cf(l)*dux13                          
         der23=der23+cf(l)*dux23                          
      endif                         
    2 continue                                     
      if (iop.eq.1) then                          
         der(1)=der12*(1.d0-pux12)*qux12        
         der(2)=der13*(1.d0-pux13)*qux13        
         der(3)=der23*(1.d0-pux23)*qux23        
      endif                         
      return
      end
*************************************************************
      subroutine abcd  (r12,r13,r14,r23,r24,r34,ener,der,iop)
*************************************************************
*     This subroutine computes the energies of a 6D PES      
*     for the ABCD system class fitted to  838 points      
*     rms =      0.19433131 kcal/mol                               
*     emax=      0.99925361 kcal/mol                               
*************************************************************
      implicit real*8(a-h,o-z)                               
      parameter(is=  241, mt_2= 3)                      
      dimension vex6(6)                                      
      dimension der(6)                                       
      dimension i1(is),i2(is),i3(is),i4(is),i5(is),i6(is),cf(is)        
      common/fun4/f12(-1:mt_2),f13(-1:mt_2),f14(-1:mt_2),    
     ,            f23(-1:mt_2),f24(-1:mt_2),f34(-1:mt_2)     
      data vex6(1)/    0.484131000000D+00/                               
      data vex6(2)/    0.753610000000D-01/                               
      data vex6(3)/    0.153391452121D+00/                               
      data vex6(4)/    0.108619100000D+01/                               
      data vex6(5)/    0.104325017226D+01/                               
      data vex6(6)/    0.141360000000D+01/                               
      data der12/0.d0/,der13/0.d0/,der14/0.d0/               
      data der23/0.d0/,der24/0.d0/,der34/0.d0/               
      data f12(-1)/0.d0/,f13(-1)/0.d0/,f14(-1)/0.d0/         
      data f23(-1)/0.d0/,f24(-1)/0.d0/,f34(-1)/0.d0/         
      data f12(0)/1.d0/,f13(0)/1.d0/,f14(0)/1.d0/            
      data f23(0)/1.d0/,f24(0)/1.d0/,f34(0)/1.d0/            
      data cf(  1)/   -0.301201699138D+03/                           
      data i1(  1)/ 0/,i2(  1)/ 0/,i3(  1)/ 1/
      data i4(  1)/ 0/,i5(  1)/ 1/,i6(  1)/ 1/
      data cf(  2)/    0.503085618249D+03/                           
      data i1(  2)/ 0/,i2(  2)/ 0/,i3(  2)/ 1/
      data i4(  2)/ 1/,i5(  2)/ 0/,i6(  2)/ 1/
      data cf(  3)/   -0.308129901665D+03/                           
      data i1(  3)/ 0/,i2(  3)/ 0/,i3(  3)/ 1/
      data i4(  3)/ 1/,i5(  3)/ 1/,i6(  3)/ 0/
      data cf(  4)/   -0.491778725023D+03/                           
      data i1(  4)/ 0/,i2(  4)/ 1/,i3(  4)/ 0/
      data i4(  4)/ 0/,i5(  4)/ 1/,i6(  4)/ 1/
      data cf(  5)/    0.620161528534D+02/                           
      data i1(  5)/ 0/,i2(  5)/ 1/,i3(  5)/ 0/
      data i4(  5)/ 1/,i5(  5)/ 0/,i6(  5)/ 1/
      data cf(  6)/    0.421995079118D+03/                           
      data i1(  6)/ 0/,i2(  6)/ 1/,i3(  6)/ 0/
      data i4(  6)/ 1/,i5(  6)/ 1/,i6(  6)/ 0/
      data cf(  7)/    0.200844176878D+02/                           
      data i1(  7)/ 0/,i2(  7)/ 1/,i3(  7)/ 1/
      data i4(  7)/ 0/,i5(  7)/ 1/,i6(  7)/ 0/
      data cf(  8)/   -0.291913776887D+02/                           
      data i1(  8)/ 0/,i2(  8)/ 1/,i3(  8)/ 1/
      data i4(  8)/ 1/,i5(  8)/ 0/,i6(  8)/ 0/
      data cf(  9)/    0.267437488347D+04/                           
      data i1(  9)/ 1/,i2(  9)/ 0/,i3(  9)/ 0/
      data i4(  9)/ 0/,i5(  9)/ 1/,i6(  9)/ 1/
      data cf( 10)/   -0.141626921759D+04/                           
      data i1( 10)/ 1/,i2( 10)/ 0/,i3( 10)/ 0/
      data i4( 10)/ 1/,i5( 10)/ 0/,i6( 10)/ 1/
      data cf( 11)/   -0.485824402595D+03/                           
      data i1( 11)/ 1/,i2( 11)/ 0/,i3( 11)/ 0/
      data i4( 11)/ 1/,i5( 11)/ 1/,i6( 11)/ 0/
      data cf( 12)/   -0.184911602604D+03/                           
      data i1( 12)/ 1/,i2( 12)/ 0/,i3( 12)/ 1/
      data i4( 12)/ 0/,i5( 12)/ 0/,i6( 12)/ 1/
      data cf( 13)/    0.872551143871D+02/                           
      data i1( 13)/ 1/,i2( 13)/ 0/,i3( 13)/ 1/
      data i4( 13)/ 1/,i5( 13)/ 0/,i6( 13)/ 0/
      data cf( 14)/    0.101365172019D+03/                           
      data i1( 14)/ 1/,i2( 14)/ 1/,i3( 14)/ 0/
      data i4( 14)/ 0/,i5( 14)/ 0/,i6( 14)/ 1/
      data cf( 15)/   -0.665044687548D+02/                           
      data i1( 15)/ 1/,i2( 15)/ 1/,i3( 15)/ 0/
      data i4( 15)/ 0/,i5( 15)/ 1/,i6( 15)/ 0/
      data cf( 16)/    0.145900405132D+01/                           
      data i1( 16)/ 1/,i2( 16)/ 1/,i3( 16)/ 1/
      data i4( 16)/ 0/,i5( 16)/ 0/,i6( 16)/ 0/
      data cf( 17)/   -0.582851974120D+03/                           
      data i1( 17)/ 0/,i2( 17)/ 0/,i3( 17)/ 1/
      data i4( 17)/ 0/,i5( 17)/ 1/,i6( 17)/ 2/
      data cf( 18)/    0.250347985393D+03/                           
      data i1( 18)/ 0/,i2( 18)/ 0/,i3( 18)/ 1/
      data i4( 18)/ 0/,i5( 18)/ 2/,i6( 18)/ 1/
      data cf( 19)/    0.500297835548D+03/                           
      data i1( 19)/ 0/,i2( 19)/ 0/,i3( 19)/ 1/
      data i4( 19)/ 1/,i5( 19)/ 0/,i6( 19)/ 2/
      data cf( 20)/   -0.799680184273D+01/                           
      data i1( 20)/ 0/,i2( 20)/ 0/,i3( 20)/ 1/
      data i4( 20)/ 1/,i5( 20)/ 1/,i6( 20)/ 1/
      data cf( 21)/    0.735389119101D+03/                           
      data i1( 21)/ 0/,i2( 21)/ 0/,i3( 21)/ 1/
      data i4( 21)/ 1/,i5( 21)/ 2/,i6( 21)/ 0/
      data cf( 22)/   -0.777074454824D+02/                           
      data i1( 22)/ 0/,i2( 22)/ 0/,i3( 22)/ 1/
      data i4( 22)/ 2/,i5( 22)/ 0/,i6( 22)/ 1/
      data cf( 23)/    0.658463276817D+02/                           
      data i1( 23)/ 0/,i2( 23)/ 0/,i3( 23)/ 1/
      data i4( 23)/ 2/,i5( 23)/ 1/,i6( 23)/ 0/
      data cf( 24)/    0.285587444375D+03/                           
      data i1( 24)/ 0/,i2( 24)/ 0/,i3( 24)/ 2/
      data i4( 24)/ 0/,i5( 24)/ 1/,i6( 24)/ 1/
      data cf( 25)/   -0.392225772595D+03/                           
      data i1( 25)/ 0/,i2( 25)/ 0/,i3( 25)/ 2/
      data i4( 25)/ 1/,i5( 25)/ 0/,i6( 25)/ 1/
      data cf( 26)/    0.215672909015D+03/                           
      data i1( 26)/ 0/,i2( 26)/ 0/,i3( 26)/ 2/
      data i4( 26)/ 1/,i5( 26)/ 1/,i6( 26)/ 0/
      data cf( 27)/    0.184036659999D+03/                           
      data i1( 27)/ 0/,i2( 27)/ 1/,i3( 27)/ 0/
      data i4( 27)/ 0/,i5( 27)/ 1/,i6( 27)/ 2/
      data cf( 28)/    0.620582042956D+03/                           
      data i1( 28)/ 0/,i2( 28)/ 1/,i3( 28)/ 0/
      data i4( 28)/ 0/,i5( 28)/ 2/,i6( 28)/ 1/
      data cf( 29)/    0.313548437956D+03/                           
      data i1( 29)/ 0/,i2( 29)/ 1/,i3( 29)/ 0/
      data i4( 29)/ 1/,i5( 29)/ 0/,i6( 29)/ 2/
      data cf( 30)/   -0.725059285802D+03/                           
      data i1( 30)/ 0/,i2( 30)/ 1/,i3( 30)/ 0/
      data i4( 30)/ 1/,i5( 30)/ 1/,i6( 30)/ 1/
      data cf( 31)/   -0.589573119764D+03/                           
      data i1( 31)/ 0/,i2( 31)/ 1/,i3( 31)/ 0/
      data i4( 31)/ 1/,i5( 31)/ 2/,i6( 31)/ 0/
      data cf( 32)/   -0.190846170371D+03/                           
      data i1( 32)/ 0/,i2( 32)/ 1/,i3( 32)/ 0/
      data i4( 32)/ 2/,i5( 32)/ 0/,i6( 32)/ 1/
      data cf( 33)/    0.410936937085D+03/                           
      data i1( 33)/ 0/,i2( 33)/ 1/,i3( 33)/ 0/
      data i4( 33)/ 2/,i5( 33)/ 1/,i6( 33)/ 0/
      data cf( 34)/    0.270399834664D+03/                           
      data i1( 34)/ 0/,i2( 34)/ 1/,i3( 34)/ 1/
      data i4( 34)/ 0/,i5( 34)/ 1/,i6( 34)/ 1/
      data cf( 35)/   -0.575724105910D+02/                           
      data i1( 35)/ 0/,i2( 35)/ 1/,i3( 35)/ 1/
      data i4( 35)/ 0/,i5( 35)/ 2/,i6( 35)/ 0/
      data cf( 36)/    0.415778502457D+02/                           
      data i1( 36)/ 0/,i2( 36)/ 1/,i3( 36)/ 1/
      data i4( 36)/ 1/,i5( 36)/ 0/,i6( 36)/ 1/
      data cf( 37)/   -0.355355596670D+03/                           
      data i1( 37)/ 0/,i2( 37)/ 1/,i3( 37)/ 1/
      data i4( 37)/ 1/,i5( 37)/ 1/,i6( 37)/ 0/
      data cf( 38)/    0.888099151454D+01/                           
      data i1( 38)/ 0/,i2( 38)/ 1/,i3( 38)/ 1/
      data i4( 38)/ 2/,i5( 38)/ 0/,i6( 38)/ 0/
      data cf( 39)/   -0.216832537486D+02/                           
      data i1( 39)/ 0/,i2( 39)/ 1/,i3( 39)/ 2/
      data i4( 39)/ 0/,i5( 39)/ 1/,i6( 39)/ 0/
      data cf( 40)/    0.192496494584D+02/                           
      data i1( 40)/ 0/,i2( 40)/ 1/,i3( 40)/ 2/
      data i4( 40)/ 1/,i5( 40)/ 0/,i6( 40)/ 0/
      data cf( 41)/   -0.674626355474D+02/                           
      data i1( 41)/ 0/,i2( 41)/ 2/,i3( 41)/ 0/
      data i4( 41)/ 0/,i5( 41)/ 1/,i6( 41)/ 1/
      data cf( 42)/    0.368520144159D+02/                           
      data i1( 42)/ 0/,i2( 42)/ 2/,i3( 42)/ 0/
      data i4( 42)/ 1/,i5( 42)/ 0/,i6( 42)/ 1/
      data cf( 43)/    0.777851980033D+02/                           
      data i1( 43)/ 0/,i2( 43)/ 2/,i3( 43)/ 0/
      data i4( 43)/ 1/,i5( 43)/ 1/,i6( 43)/ 0/
      data cf( 44)/    0.653680556725D+01/                           
      data i1( 44)/ 0/,i2( 44)/ 2/,i3( 44)/ 1/
      data i4( 44)/ 0/,i5( 44)/ 1/,i6( 44)/ 0/
      data cf( 45)/   -0.207387331471D+01/                           
      data i1( 45)/ 0/,i2( 45)/ 2/,i3( 45)/ 1/
      data i4( 45)/ 1/,i5( 45)/ 0/,i6( 45)/ 0/
      data cf( 46)/    0.797023433813D+03/                           
      data i1( 46)/ 1/,i2( 46)/ 0/,i3( 46)/ 0/
      data i4( 46)/ 0/,i5( 46)/ 1/,i6( 46)/ 2/
      data cf( 47)/   -0.309402008424D+04/                           
      data i1( 47)/ 1/,i2( 47)/ 0/,i3( 47)/ 0/
      data i4( 47)/ 0/,i5( 47)/ 2/,i6( 47)/ 1/
      data cf( 48)/   -0.278552713331D+04/                           
      data i1( 48)/ 1/,i2( 48)/ 0/,i3( 48)/ 0/
      data i4( 48)/ 1/,i5( 48)/ 0/,i6( 48)/ 2/
      data cf( 49)/    0.230249148462D+04/                           
      data i1( 49)/ 1/,i2( 49)/ 0/,i3( 49)/ 0/
      data i4( 49)/ 1/,i5( 49)/ 1/,i6( 49)/ 1/
      data cf( 50)/   -0.259555920681D+03/                           
      data i1( 50)/ 1/,i2( 50)/ 0/,i3( 50)/ 0/
      data i4( 50)/ 1/,i5( 50)/ 2/,i6( 50)/ 0/
      data cf( 51)/    0.614933210208D+03/                           
      data i1( 51)/ 1/,i2( 51)/ 0/,i3( 51)/ 0/
      data i4( 51)/ 2/,i5( 51)/ 0/,i6( 51)/ 1/
      data cf( 52)/   -0.165239581014D+04/                           
      data i1( 52)/ 1/,i2( 52)/ 0/,i3( 52)/ 0/
      data i4( 52)/ 2/,i5( 52)/ 1/,i6( 52)/ 0/
      data cf( 53)/    0.237183053979D+03/                           
      data i1( 53)/ 1/,i2( 53)/ 0/,i3( 53)/ 1/
      data i4( 53)/ 0/,i5( 53)/ 0/,i6( 53)/ 2/
      data cf( 54)/   -0.155315513955D+04/                           
      data i1( 54)/ 1/,i2( 54)/ 0/,i3( 54)/ 1/
      data i4( 54)/ 0/,i5( 54)/ 1/,i6( 54)/ 1/
      data cf( 55)/    0.559039909252D+03/                           
      data i1( 55)/ 1/,i2( 55)/ 0/,i3( 55)/ 1/
      data i4( 55)/ 1/,i5( 55)/ 0/,i6( 55)/ 1/
      data cf( 56)/    0.495849411797D+03/                           
      data i1( 56)/ 1/,i2( 56)/ 0/,i3( 56)/ 1/
      data i4( 56)/ 1/,i5( 56)/ 1/,i6( 56)/ 0/
      data cf( 57)/   -0.431983769097D+02/                           
      data i1( 57)/ 1/,i2( 57)/ 0/,i3( 57)/ 1/
      data i4( 57)/ 2/,i5( 57)/ 0/,i6( 57)/ 0/
      data cf( 58)/    0.123730578684D+03/                           
      data i1( 58)/ 1/,i2( 58)/ 0/,i3( 58)/ 2/
      data i4( 58)/ 0/,i5( 58)/ 0/,i6( 58)/ 1/
      data cf( 59)/   -0.468635303084D+02/                           
      data i1( 59)/ 1/,i2( 59)/ 0/,i3( 59)/ 2/
      data i4( 59)/ 1/,i5( 59)/ 0/,i6( 59)/ 0/
      data cf( 60)/   -0.126758283940D+03/                           
      data i1( 60)/ 1/,i2( 60)/ 1/,i3( 60)/ 0/
      data i4( 60)/ 0/,i5( 60)/ 0/,i6( 60)/ 2/
      data cf( 61)/    0.920555511588D+03/                           
      data i1( 61)/ 1/,i2( 61)/ 1/,i3( 61)/ 0/
      data i4( 61)/ 0/,i5( 61)/ 1/,i6( 61)/ 1/
      data cf( 62)/    0.179226853732D+03/                           
      data i1( 62)/ 1/,i2( 62)/ 1/,i3( 62)/ 0/
      data i4( 62)/ 0/,i5( 62)/ 2/,i6( 62)/ 0/
      data cf( 63)/   -0.508236074336D+03/                           
      data i1( 63)/ 1/,i2( 63)/ 1/,i3( 63)/ 0/
      data i4( 63)/ 1/,i5( 63)/ 0/,i6( 63)/ 1/
      data cf( 64)/   -0.501761193027D+03/                           
      data i1( 64)/ 1/,i2( 64)/ 1/,i3( 64)/ 0/
      data i4( 64)/ 1/,i5( 64)/ 1/,i6( 64)/ 0/
      data cf( 65)/   -0.693673531433D+02/                           
      data i1( 65)/ 1/,i2( 65)/ 1/,i3( 65)/ 1/
      data i4( 65)/ 0/,i5( 65)/ 0/,i6( 65)/ 1/
      data cf( 66)/    0.643586325783D+02/                           
      data i1( 66)/ 1/,i2( 66)/ 1/,i3( 66)/ 1/
      data i4( 66)/ 0/,i5( 66)/ 1/,i6( 66)/ 0/
      data cf( 67)/    0.254145284262D+02/                           
      data i1( 67)/ 1/,i2( 67)/ 1/,i3( 67)/ 1/
      data i4( 67)/ 1/,i5( 67)/ 0/,i6( 67)/ 0/
      data cf( 68)/   -0.969049668679D+00/                           
      data i1( 68)/ 1/,i2( 68)/ 1/,i3( 68)/ 2/
      data i4( 68)/ 0/,i5( 68)/ 0/,i6( 68)/ 0/
      data cf( 69)/   -0.739606491236D+01/                           
      data i1( 69)/ 1/,i2( 69)/ 2/,i3( 69)/ 0/
      data i4( 69)/ 0/,i5( 69)/ 0/,i6( 69)/ 1/
      data cf( 70)/   -0.216181659296D+02/                           
      data i1( 70)/ 1/,i2( 70)/ 2/,i3( 70)/ 0/
      data i4( 70)/ 0/,i5( 70)/ 1/,i6( 70)/ 0/
      data cf( 71)/    0.113301215517D+00/                           
      data i1( 71)/ 1/,i2( 71)/ 2/,i3( 71)/ 1/
      data i4( 71)/ 0/,i5( 71)/ 0/,i6( 71)/ 0/
      data cf( 72)/   -0.316237980630D+04/                           
      data i1( 72)/ 2/,i2( 72)/ 0/,i3( 72)/ 0/
      data i4( 72)/ 0/,i5( 72)/ 1/,i6( 72)/ 1/
      data cf( 73)/    0.220256962488D+04/                           
      data i1( 73)/ 2/,i2( 73)/ 0/,i3( 73)/ 0/
      data i4( 73)/ 1/,i5( 73)/ 0/,i6( 73)/ 1/
      data cf( 74)/    0.996923321660D+03/                           
      data i1( 74)/ 2/,i2( 74)/ 0/,i3( 74)/ 0/
      data i4( 74)/ 1/,i5( 74)/ 1/,i6( 74)/ 0/
      data cf( 75)/    0.261067569506D+02/                           
      data i1( 75)/ 2/,i2( 75)/ 0/,i3( 75)/ 1/
      data i4( 75)/ 0/,i5( 75)/ 0/,i6( 75)/ 1/
      data cf( 76)/   -0.869406806611D+02/                           
      data i1( 76)/ 2/,i2( 76)/ 0/,i3( 76)/ 1/
      data i4( 76)/ 1/,i5( 76)/ 0/,i6( 76)/ 0/
      data cf( 77)/    0.121756763089D+02/                           
      data i1( 77)/ 2/,i2( 77)/ 1/,i3( 77)/ 0/
      data i4( 77)/ 0/,i5( 77)/ 0/,i6( 77)/ 1/
      data cf( 78)/    0.205900179297D+02/                           
      data i1( 78)/ 2/,i2( 78)/ 1/,i3( 78)/ 0/
      data i4( 78)/ 0/,i5( 78)/ 1/,i6( 78)/ 0/
      data cf( 79)/   -0.103115817598D+01/                           
      data i1( 79)/ 2/,i2( 79)/ 1/,i3( 79)/ 1/
      data i4( 79)/ 0/,i5( 79)/ 0/,i6( 79)/ 0/
      data cf( 80)/    0.142851263704D+04/                           
      data i1( 80)/ 0/,i2( 80)/ 0/,i3( 80)/ 1/
      data i4( 80)/ 0/,i5( 80)/ 1/,i6( 80)/ 3/
      data cf( 81)/   -0.751225580976D+03/                           
      data i1( 81)/ 0/,i2( 81)/ 0/,i3( 81)/ 1/
      data i4( 81)/ 0/,i5( 81)/ 2/,i6( 81)/ 2/
      data cf( 82)/    0.442187165030D+03/                           
      data i1( 82)/ 0/,i2( 82)/ 0/,i3( 82)/ 1/
      data i4( 82)/ 0/,i5( 82)/ 3/,i6( 82)/ 1/
      data cf( 83)/   -0.871969713013D+03/                           
      data i1( 83)/ 0/,i2( 83)/ 0/,i3( 83)/ 1/
      data i4( 83)/ 1/,i5( 83)/ 0/,i6( 83)/ 3/
      data cf( 84)/   -0.568130029804D+03/                           
      data i1( 84)/ 0/,i2( 84)/ 0/,i3( 84)/ 1/
      data i4( 84)/ 1/,i5( 84)/ 1/,i6( 84)/ 2/
      data cf( 85)/    0.279786741204D+03/                           
      data i1( 85)/ 0/,i2( 85)/ 0/,i3( 85)/ 1/
      data i4( 85)/ 1/,i5( 85)/ 2/,i6( 85)/ 1/
      data cf( 86)/   -0.123688831099D+03/                           
      data i1( 86)/ 0/,i2( 86)/ 0/,i3( 86)/ 1/
      data i4( 86)/ 1/,i5( 86)/ 3/,i6( 86)/ 0/
      data cf( 87)/    0.230191734525D+03/                           
      data i1( 87)/ 0/,i2( 87)/ 0/,i3( 87)/ 1/
      data i4( 87)/ 2/,i5( 87)/ 0/,i6( 87)/ 2/
      data cf( 88)/    0.103698495128D+03/                           
      data i1( 88)/ 0/,i2( 88)/ 0/,i3( 88)/ 1/
      data i4( 88)/ 2/,i5( 88)/ 1/,i6( 88)/ 1/
      data cf( 89)/   -0.778028946643D+03/                           
      data i1( 89)/ 0/,i2( 89)/ 0/,i3( 89)/ 1/
      data i4( 89)/ 2/,i5( 89)/ 2/,i6( 89)/ 0/
      data cf( 90)/   -0.487606032503D+03/                           
      data i1( 90)/ 0/,i2( 90)/ 0/,i3( 90)/ 1/
      data i4( 90)/ 3/,i5( 90)/ 0/,i6( 90)/ 1/
      data cf( 91)/    0.538870961672D+03/                           
      data i1( 91)/ 0/,i2( 91)/ 0/,i3( 91)/ 1/
      data i4( 91)/ 3/,i5( 91)/ 1/,i6( 91)/ 0/
      data cf( 92)/   -0.959523723116D+02/                           
      data i1( 92)/ 0/,i2( 92)/ 0/,i3( 92)/ 2/
      data i4( 92)/ 0/,i5( 92)/ 1/,i6( 92)/ 2/
      data cf( 93)/   -0.501318653907D+02/                           
      data i1( 93)/ 0/,i2( 93)/ 0/,i3( 93)/ 2/
      data i4( 93)/ 0/,i5( 93)/ 2/,i6( 93)/ 1/
      data cf( 94)/    0.572963536518D+02/                           
      data i1( 94)/ 0/,i2( 94)/ 0/,i3( 94)/ 2/
      data i4( 94)/ 1/,i5( 94)/ 0/,i6( 94)/ 2/
      data cf( 95)/    0.288250201934D+03/                           
      data i1( 95)/ 0/,i2( 95)/ 0/,i3( 95)/ 2/
      data i4( 95)/ 1/,i5( 95)/ 1/,i6( 95)/ 1/
      data cf( 96)/   -0.304627556347D+03/                           
      data i1( 96)/ 0/,i2( 96)/ 0/,i3( 96)/ 2/
      data i4( 96)/ 1/,i5( 96)/ 2/,i6( 96)/ 0/
      data cf( 97)/    0.100680385451D+03/                           
      data i1( 97)/ 0/,i2( 97)/ 0/,i3( 97)/ 2/
      data i4( 97)/ 2/,i5( 97)/ 0/,i6( 97)/ 1/
      data cf( 98)/   -0.131480273166D+03/                           
      data i1( 98)/ 0/,i2( 98)/ 0/,i3( 98)/ 2/
      data i4( 98)/ 2/,i5( 98)/ 1/,i6( 98)/ 0/
      data cf( 99)/   -0.560717387084D+02/                           
      data i1( 99)/ 0/,i2( 99)/ 0/,i3( 99)/ 3/
      data i4( 99)/ 0/,i5( 99)/ 1/,i6( 99)/ 1/
      data cf(100)/    0.662930234847D+02/                           
      data i1(100)/ 0/,i2(100)/ 0/,i3(100)/ 3/
      data i4(100)/ 1/,i5(100)/ 0/,i6(100)/ 1/
      data cf(101)/   -0.559701510309D+02/                           
      data i1(101)/ 0/,i2(101)/ 0/,i3(101)/ 3/
      data i4(101)/ 1/,i5(101)/ 1/,i6(101)/ 0/
      data cf(102)/   -0.663307490043D+02/                           
      data i1(102)/ 0/,i2(102)/ 1/,i3(102)/ 0/
      data i4(102)/ 0/,i5(102)/ 1/,i6(102)/ 3/
      data cf(103)/    0.148849726140D+03/                           
      data i1(103)/ 0/,i2(103)/ 1/,i3(103)/ 0/
      data i4(103)/ 0/,i5(103)/ 2/,i6(103)/ 2/
      data cf(104)/   -0.427087035348D+03/                           
      data i1(104)/ 0/,i2(104)/ 1/,i3(104)/ 0/
      data i4(104)/ 0/,i5(104)/ 3/,i6(104)/ 1/
      data cf(105)/    0.109568754742D+03/                           
      data i1(105)/ 0/,i2(105)/ 1/,i3(105)/ 0/
      data i4(105)/ 1/,i5(105)/ 0/,i6(105)/ 3/
      data cf(106)/   -0.111247804332D+02/                           
      data i1(106)/ 0/,i2(106)/ 1/,i3(106)/ 0/
      data i4(106)/ 1/,i5(106)/ 1/,i6(106)/ 2/
      data cf(107)/    0.243415539861D+03/                           
      data i1(107)/ 0/,i2(107)/ 1/,i3(107)/ 0/
      data i4(107)/ 1/,i5(107)/ 2/,i6(107)/ 1/
      data cf(108)/    0.419396033354D+03/                           
      data i1(108)/ 0/,i2(108)/ 1/,i3(108)/ 0/
      data i4(108)/ 1/,i5(108)/ 3/,i6(108)/ 0/
      data cf(109)/   -0.464944550410D+03/                           
      data i1(109)/ 0/,i2(109)/ 1/,i3(109)/ 0/
      data i4(109)/ 2/,i5(109)/ 0/,i6(109)/ 2/
      data cf(110)/    0.539402870059D+03/                           
      data i1(110)/ 0/,i2(110)/ 1/,i3(110)/ 0/
      data i4(110)/ 2/,i5(110)/ 1/,i6(110)/ 1/
      data cf(111)/   -0.334638894826D+03/                           
      data i1(111)/ 0/,i2(111)/ 1/,i3(111)/ 0/
      data i4(111)/ 2/,i5(111)/ 2/,i6(111)/ 0/
      data cf(112)/    0.302989132379D+03/                           
      data i1(112)/ 0/,i2(112)/ 1/,i3(112)/ 0/
      data i4(112)/ 3/,i5(112)/ 0/,i6(112)/ 1/
      data cf(113)/   -0.337865480532D+03/                           
      data i1(113)/ 0/,i2(113)/ 1/,i3(113)/ 0/
      data i4(113)/ 3/,i5(113)/ 1/,i6(113)/ 0/
      data cf(114)/    0.300574120469D+02/                           
      data i1(114)/ 0/,i2(114)/ 1/,i3(114)/ 1/
      data i4(114)/ 0/,i5(114)/ 1/,i6(114)/ 2/
      data cf(115)/   -0.211040230639D+03/                           
      data i1(115)/ 0/,i2(115)/ 1/,i3(115)/ 1/
      data i4(115)/ 0/,i5(115)/ 2/,i6(115)/ 1/
      data cf(116)/   -0.183587807022D+02/                           
      data i1(116)/ 0/,i2(116)/ 1/,i3(116)/ 1/
      data i4(116)/ 0/,i5(116)/ 3/,i6(116)/ 0/
      data cf(117)/   -0.394306259788D+02/                           
      data i1(117)/ 0/,i2(117)/ 1/,i3(117)/ 1/
      data i4(117)/ 1/,i5(117)/ 0/,i6(117)/ 2/
      data cf(118)/    0.328329861799D+01/                           
      data i1(118)/ 0/,i2(118)/ 1/,i3(118)/ 1/
      data i4(118)/ 1/,i5(118)/ 1/,i6(118)/ 1/
      data cf(119)/    0.322275533313D+03/                           
      data i1(119)/ 0/,i2(119)/ 1/,i3(119)/ 1/
      data i4(119)/ 1/,i5(119)/ 2/,i6(119)/ 0/
      data cf(120)/   -0.160289354568D+02/                           
      data i1(120)/ 0/,i2(120)/ 1/,i3(120)/ 1/
      data i4(120)/ 2/,i5(120)/ 0/,i6(120)/ 1/
      data cf(121)/   -0.534836566747D+01/                           
      data i1(121)/ 0/,i2(121)/ 1/,i3(121)/ 1/
      data i4(121)/ 2/,i5(121)/ 1/,i6(121)/ 0/
      data cf(122)/   -0.344831289197D+01/                           
      data i1(122)/ 0/,i2(122)/ 1/,i3(122)/ 1/
      data i4(122)/ 3/,i5(122)/ 0/,i6(122)/ 0/
      data cf(123)/   -0.382638257656D+02/                           
      data i1(123)/ 0/,i2(123)/ 1/,i3(123)/ 2/
      data i4(123)/ 0/,i5(123)/ 1/,i6(123)/ 1/
      data cf(124)/    0.188525384565D+02/                           
      data i1(124)/ 0/,i2(124)/ 1/,i3(124)/ 2/
      data i4(124)/ 0/,i5(124)/ 2/,i6(124)/ 0/
      data cf(125)/   -0.173822073182D+02/                           
      data i1(125)/ 0/,i2(125)/ 1/,i3(125)/ 2/
      data i4(125)/ 1/,i5(125)/ 0/,i6(125)/ 1/
      data cf(126)/    0.821354775951D+02/                           
      data i1(126)/ 0/,i2(126)/ 1/,i3(126)/ 2/
      data i4(126)/ 1/,i5(126)/ 1/,i6(126)/ 0/
      data cf(127)/    0.528809319761D+00/                           
      data i1(127)/ 0/,i2(127)/ 1/,i3(127)/ 2/
      data i4(127)/ 2/,i5(127)/ 0/,i6(127)/ 0/
      data cf(128)/    0.530664448902D+01/                           
      data i1(128)/ 0/,i2(128)/ 1/,i3(128)/ 3/
      data i4(128)/ 0/,i5(128)/ 1/,i6(128)/ 0/
      data cf(129)/   -0.377776819510D+01/                           
      data i1(129)/ 0/,i2(129)/ 1/,i3(129)/ 3/
      data i4(129)/ 1/,i5(129)/ 0/,i6(129)/ 0/
      data cf(130)/   -0.259855504868D+01/                           
      data i1(130)/ 0/,i2(130)/ 2/,i3(130)/ 0/
      data i4(130)/ 0/,i5(130)/ 1/,i6(130)/ 2/
      data cf(131)/    0.288684849214D+02/                           
      data i1(131)/ 0/,i2(131)/ 2/,i3(131)/ 0/
      data i4(131)/ 0/,i5(131)/ 2/,i6(131)/ 1/
      data cf(132)/   -0.796958242747D+02/                           
      data i1(132)/ 0/,i2(132)/ 2/,i3(132)/ 0/
      data i4(132)/ 1/,i5(132)/ 0/,i6(132)/ 2/
      data cf(133)/    0.858401437291D+02/                           
      data i1(133)/ 0/,i2(133)/ 2/,i3(133)/ 0/
      data i4(133)/ 1/,i5(133)/ 1/,i6(133)/ 1/
      data cf(134)/   -0.322531397013D+02/                           
      data i1(134)/ 0/,i2(134)/ 2/,i3(134)/ 0/
      data i4(134)/ 1/,i5(134)/ 2/,i6(134)/ 0/
      data cf(135)/    0.552223496315D+01/                           
      data i1(135)/ 0/,i2(135)/ 2/,i3(135)/ 0/
      data i4(135)/ 2/,i5(135)/ 0/,i6(135)/ 1/
      data cf(136)/   -0.844112164944D+02/                           
      data i1(136)/ 0/,i2(136)/ 2/,i3(136)/ 0/
      data i4(136)/ 2/,i5(136)/ 1/,i6(136)/ 0/
      data cf(137)/    0.750597406673D+01/                           
      data i1(137)/ 0/,i2(137)/ 2/,i3(137)/ 1/
      data i4(137)/ 0/,i5(137)/ 1/,i6(137)/ 1/
      data cf(138)/   -0.187160231524D+00/                           
      data i1(138)/ 0/,i2(138)/ 2/,i3(138)/ 1/
      data i4(138)/ 0/,i5(138)/ 2/,i6(138)/ 0/
      data cf(139)/   -0.582489388201D+01/                           
      data i1(139)/ 0/,i2(139)/ 2/,i3(139)/ 1/
      data i4(139)/ 1/,i5(139)/ 0/,i6(139)/ 1/
      data cf(140)/   -0.174553794688D+02/                           
      data i1(140)/ 0/,i2(140)/ 2/,i3(140)/ 1/
      data i4(140)/ 1/,i5(140)/ 1/,i6(140)/ 0/
      data cf(141)/   -0.167749728558D+01/                           
      data i1(141)/ 0/,i2(141)/ 2/,i3(141)/ 1/
      data i4(141)/ 2/,i5(141)/ 0/,i6(141)/ 0/
      data cf(142)/   -0.221741374534D+01/                           
      data i1(142)/ 0/,i2(142)/ 2/,i3(142)/ 2/
      data i4(142)/ 0/,i5(142)/ 1/,i6(142)/ 0/
      data cf(143)/    0.111211601043D+01/                           
      data i1(143)/ 0/,i2(143)/ 2/,i3(143)/ 2/
      data i4(143)/ 1/,i5(143)/ 0/,i6(143)/ 0/
      data cf(144)/    0.428751102711D+01/                           
      data i1(144)/ 0/,i2(144)/ 3/,i3(144)/ 0/
      data i4(144)/ 0/,i5(144)/ 1/,i6(144)/ 1/
      data cf(145)/   -0.171798633623D+01/                           
      data i1(145)/ 0/,i2(145)/ 3/,i3(145)/ 0/
      data i4(145)/ 1/,i5(145)/ 0/,i6(145)/ 1/
      data cf(146)/   -0.496354677492D+01/                           
      data i1(146)/ 0/,i2(146)/ 3/,i3(146)/ 0/
      data i4(146)/ 1/,i5(146)/ 1/,i6(146)/ 0/
      data cf(147)/    0.115915301209D+00/                           
      data i1(147)/ 0/,i2(147)/ 3/,i3(147)/ 1/
      data i4(147)/ 0/,i5(147)/ 1/,i6(147)/ 0/
      data cf(148)/   -0.166608600732D+00/                           
      data i1(148)/ 0/,i2(148)/ 3/,i3(148)/ 1/
      data i4(148)/ 1/,i5(148)/ 0/,i6(148)/ 0/
      data cf(149)/   -0.452817339529D+04/                           
      data i1(149)/ 1/,i2(149)/ 0/,i3(149)/ 0/
      data i4(149)/ 0/,i5(149)/ 1/,i6(149)/ 3/
      data cf(150)/    0.127923663390D+04/                           
      data i1(150)/ 1/,i2(150)/ 0/,i3(150)/ 0/
      data i4(150)/ 0/,i5(150)/ 2/,i6(150)/ 2/
      data cf(151)/    0.936265037109D+03/                           
      data i1(151)/ 1/,i2(151)/ 0/,i3(151)/ 0/
      data i4(151)/ 0/,i5(151)/ 3/,i6(151)/ 1/
      data cf(152)/    0.255438290437D+04/                           
      data i1(152)/ 1/,i2(152)/ 0/,i3(152)/ 0/
      data i4(152)/ 1/,i5(152)/ 0/,i6(152)/ 3/
      data cf(153)/    0.206842157255D+04/                           
      data i1(153)/ 1/,i2(153)/ 0/,i3(153)/ 0/
      data i4(153)/ 1/,i5(153)/ 1/,i6(153)/ 2/
      data cf(154)/   -0.181951992472D+04/                           
      data i1(154)/ 1/,i2(154)/ 0/,i3(154)/ 0/
      data i4(154)/ 1/,i5(154)/ 2/,i6(154)/ 1/
      data cf(155)/   -0.116984836414D+04/                           
      data i1(155)/ 1/,i2(155)/ 0/,i3(155)/ 0/
      data i4(155)/ 1/,i5(155)/ 3/,i6(155)/ 0/
      data cf(156)/    0.394437303290D+03/                           
      data i1(156)/ 1/,i2(156)/ 0/,i3(156)/ 0/
      data i4(156)/ 2/,i5(156)/ 0/,i6(156)/ 2/
      data cf(157)/   -0.190013725929D+04/                           
      data i1(157)/ 1/,i2(157)/ 0/,i3(157)/ 0/
      data i4(157)/ 2/,i5(157)/ 1/,i6(157)/ 1/
      data cf(158)/    0.355863756452D+04/                           
      data i1(158)/ 1/,i2(158)/ 0/,i3(158)/ 0/
      data i4(158)/ 2/,i5(158)/ 2/,i6(158)/ 0/
      data cf(159)/    0.685688035406D+03/                           
      data i1(159)/ 1/,i2(159)/ 0/,i3(159)/ 0/
      data i4(159)/ 3/,i5(159)/ 0/,i6(159)/ 1/
      data cf(160)/   -0.474595115803D+03/                           
      data i1(160)/ 1/,i2(160)/ 0/,i3(160)/ 0/
      data i4(160)/ 3/,i5(160)/ 1/,i6(160)/ 0/
      data cf(161)/   -0.279570477702D+01/                           
      data i1(161)/ 1/,i2(161)/ 0/,i3(161)/ 1/
      data i4(161)/ 0/,i5(161)/ 0/,i6(161)/ 3/
      data cf(162)/    0.736294124868D+03/                           
      data i1(162)/ 1/,i2(162)/ 0/,i3(162)/ 1/
      data i4(162)/ 0/,i5(162)/ 1/,i6(162)/ 2/
      data cf(163)/    0.705413701454D+03/                           
      data i1(163)/ 1/,i2(163)/ 0/,i3(163)/ 1/
      data i4(163)/ 0/,i5(163)/ 2/,i6(163)/ 1/
      data cf(164)/   -0.408276119509D+03/                           
      data i1(164)/ 1/,i2(164)/ 0/,i3(164)/ 1/
      data i4(164)/ 1/,i5(164)/ 0/,i6(164)/ 2/
      data cf(165)/   -0.150347457011D+04/                           
      data i1(165)/ 1/,i2(165)/ 0/,i3(165)/ 1/
      data i4(165)/ 1/,i5(165)/ 1/,i6(165)/ 1/
      data cf(166)/    0.138477465722D+03/                           
      data i1(166)/ 1/,i2(166)/ 0/,i3(166)/ 1/
      data i4(166)/ 1/,i5(166)/ 2/,i6(166)/ 0/
      data cf(167)/   -0.223733190459D+03/                           
      data i1(167)/ 1/,i2(167)/ 0/,i3(167)/ 1/
      data i4(167)/ 2/,i5(167)/ 0/,i6(167)/ 1/
      data cf(168)/    0.492412839903D+03/                           
      data i1(168)/ 1/,i2(168)/ 0/,i3(168)/ 1/
      data i4(168)/ 2/,i5(168)/ 1/,i6(168)/ 0/
      data cf(169)/    0.373080591583D+02/                           
      data i1(169)/ 1/,i2(169)/ 0/,i3(169)/ 1/
      data i4(169)/ 3/,i5(169)/ 0/,i6(169)/ 0/
      data cf(170)/   -0.469912835786D+02/                           
      data i1(170)/ 1/,i2(170)/ 0/,i3(170)/ 2/
      data i4(170)/ 0/,i5(170)/ 0/,i6(170)/ 2/
      data cf(171)/    0.193874621446D+03/                           
      data i1(171)/ 1/,i2(171)/ 0/,i3(171)/ 2/
      data i4(171)/ 0/,i5(171)/ 1/,i6(171)/ 1/
      data cf(172)/   -0.592907064668D+02/                           
      data i1(172)/ 1/,i2(172)/ 0/,i3(172)/ 2/
      data i4(172)/ 1/,i5(172)/ 0/,i6(172)/ 1/
      data cf(173)/    0.229239476008D+02/                           
      data i1(173)/ 1/,i2(173)/ 0/,i3(173)/ 2/
      data i4(173)/ 1/,i5(173)/ 1/,i6(173)/ 0/
      data cf(174)/   -0.347100411016D+01/                           
      data i1(174)/ 1/,i2(174)/ 0/,i3(174)/ 2/
      data i4(174)/ 2/,i5(174)/ 0/,i6(174)/ 0/
      data cf(175)/   -0.222187119257D+02/                           
      data i1(175)/ 1/,i2(175)/ 0/,i3(175)/ 3/
      data i4(175)/ 0/,i5(175)/ 0/,i6(175)/ 1/
      data cf(176)/    0.105703529351D+02/                           
      data i1(176)/ 1/,i2(176)/ 0/,i3(176)/ 3/
      data i4(176)/ 1/,i5(176)/ 0/,i6(176)/ 0/
      data cf(177)/    0.596151542199D+01/                           
      data i1(177)/ 1/,i2(177)/ 1/,i3(177)/ 0/
      data i4(177)/ 0/,i5(177)/ 0/,i6(177)/ 3/
      data cf(178)/   -0.342999439483D+03/                           
      data i1(178)/ 1/,i2(178)/ 1/,i3(178)/ 0/
      data i4(178)/ 0/,i5(178)/ 1/,i6(178)/ 2/
      data cf(179)/   -0.315999239119D+03/                           
      data i1(179)/ 1/,i2(179)/ 1/,i3(179)/ 0/
      data i4(179)/ 0/,i5(179)/ 2/,i6(179)/ 1/
      data cf(180)/    0.529584992667D+02/                           
      data i1(180)/ 1/,i2(180)/ 1/,i3(180)/ 0/
      data i4(180)/ 0/,i5(180)/ 3/,i6(180)/ 0/
      data cf(181)/    0.418058939049D+03/                           
      data i1(181)/ 1/,i2(181)/ 1/,i3(181)/ 0/
      data i4(181)/ 1/,i5(181)/ 0/,i6(181)/ 2/
      data cf(182)/    0.434824041322D+02/                           
      data i1(182)/ 1/,i2(182)/ 1/,i3(182)/ 0/
      data i4(182)/ 1/,i5(182)/ 1/,i6(182)/ 1/
      data cf(183)/   -0.546438497066D+02/                           
      data i1(183)/ 1/,i2(183)/ 1/,i3(183)/ 0/
      data i4(183)/ 1/,i5(183)/ 2/,i6(183)/ 0/
      data cf(184)/    0.173068936654D+03/                           
      data i1(184)/ 1/,i2(184)/ 1/,i3(184)/ 0/
      data i4(184)/ 2/,i5(184)/ 0/,i6(184)/ 1/
      data cf(185)/    0.301871922789D+03/                           
      data i1(185)/ 1/,i2(185)/ 1/,i3(185)/ 0/
      data i4(185)/ 2/,i5(185)/ 1/,i6(185)/ 0/
      data cf(186)/    0.492977220240D+01/                           
      data i1(186)/ 1/,i2(186)/ 1/,i3(186)/ 1/
      data i4(186)/ 0/,i5(186)/ 0/,i6(186)/ 2/
      data cf(187)/   -0.119234202340D+03/                           
      data i1(187)/ 1/,i2(187)/ 1/,i3(187)/ 1/
      data i4(187)/ 0/,i5(187)/ 1/,i6(187)/ 1/
      data cf(188)/   -0.459103547407D+02/                           
      data i1(188)/ 1/,i2(188)/ 1/,i3(188)/ 1/
      data i4(188)/ 0/,i5(188)/ 2/,i6(188)/ 0/
      data cf(189)/    0.101760338408D+03/                           
      data i1(189)/ 1/,i2(189)/ 1/,i3(189)/ 1/
      data i4(189)/ 1/,i5(189)/ 0/,i6(189)/ 1/
      data cf(190)/    0.931935866221D+01/                           
      data i1(190)/ 1/,i2(190)/ 1/,i3(190)/ 1/
      data i4(190)/ 1/,i5(190)/ 1/,i6(190)/ 0/
      data cf(191)/    0.240509240303D+01/                           
      data i1(191)/ 1/,i2(191)/ 1/,i3(191)/ 1/
      data i4(191)/ 2/,i5(191)/ 0/,i6(191)/ 0/
      data cf(192)/    0.159077541514D+02/                           
      data i1(192)/ 1/,i2(192)/ 1/,i3(192)/ 2/
      data i4(192)/ 0/,i5(192)/ 0/,i6(192)/ 1/
      data cf(193)/   -0.149569985785D+02/                           
      data i1(193)/ 1/,i2(193)/ 1/,i3(193)/ 2/
      data i4(193)/ 0/,i5(193)/ 1/,i6(193)/ 0/
      data cf(194)/   -0.901463304374D+01/                           
      data i1(194)/ 1/,i2(194)/ 1/,i3(194)/ 2/
      data i4(194)/ 1/,i5(194)/ 0/,i6(194)/ 0/
      data cf(195)/    0.155055676376D+00/                           
      data i1(195)/ 1/,i2(195)/ 1/,i3(195)/ 3/
      data i4(195)/ 0/,i5(195)/ 0/,i6(195)/ 0/
      data cf(196)/    0.145626888246D+02/                           
      data i1(196)/ 1/,i2(196)/ 2/,i3(196)/ 0/
      data i4(196)/ 0/,i5(196)/ 0/,i6(196)/ 2/
      data cf(197)/   -0.684544564016D+01/                           
      data i1(197)/ 1/,i2(197)/ 2/,i3(197)/ 0/
      data i4(197)/ 0/,i5(197)/ 1/,i6(197)/ 1/
      data cf(198)/    0.142443320598D+01/                           
      data i1(198)/ 1/,i2(198)/ 2/,i3(198)/ 0/
      data i4(198)/ 0/,i5(198)/ 2/,i6(198)/ 0/
      data cf(199)/    0.474030903046D+01/                           
      data i1(199)/ 1/,i2(199)/ 2/,i3(199)/ 0/
      data i4(199)/ 1/,i5(199)/ 0/,i6(199)/ 1/
      data cf(200)/    0.310077170419D+02/                           
      data i1(200)/ 1/,i2(200)/ 2/,i3(200)/ 0/
      data i4(200)/ 1/,i5(200)/ 1/,i6(200)/ 0/
      data cf(201)/   -0.268344236699D+01/                           
      data i1(201)/ 1/,i2(201)/ 2/,i3(201)/ 1/
      data i4(201)/ 0/,i5(201)/ 0/,i6(201)/ 1/
      data cf(202)/    0.491101814121D+01/                           
      data i1(202)/ 1/,i2(202)/ 2/,i3(202)/ 1/
      data i4(202)/ 0/,i5(202)/ 1/,i6(202)/ 0/
      data cf(203)/    0.132295364419D+01/                           
      data i1(203)/ 1/,i2(203)/ 2/,i3(203)/ 1/
      data i4(203)/ 1/,i5(203)/ 0/,i6(203)/ 0/
      data cf(204)/   -0.293091475025D-01/                           
      data i1(204)/ 1/,i2(204)/ 2/,i3(204)/ 2/
      data i4(204)/ 0/,i5(204)/ 0/,i6(204)/ 0/
      data cf(205)/    0.153972918092D+01/                           
      data i1(205)/ 1/,i2(205)/ 3/,i3(205)/ 0/
      data i4(205)/ 0/,i5(205)/ 0/,i6(205)/ 1/
      data cf(206)/   -0.224264670551D+00/                           
      data i1(206)/ 1/,i2(206)/ 3/,i3(206)/ 0/
      data i4(206)/ 0/,i5(206)/ 1/,i6(206)/ 0/
      data cf(207)/    0.534148240518D-02/                           
      data i1(207)/ 1/,i2(207)/ 3/,i3(207)/ 1/
      data i4(207)/ 0/,i5(207)/ 0/,i6(207)/ 0/
      data cf(208)/    0.975361739111D+02/                           
      data i1(208)/ 2/,i2(208)/ 0/,i3(208)/ 0/
      data i4(208)/ 0/,i5(208)/ 1/,i6(208)/ 2/
      data cf(209)/    0.134575328462D+04/                           
      data i1(209)/ 2/,i2(209)/ 0/,i3(209)/ 0/
      data i4(209)/ 0/,i5(209)/ 2/,i6(209)/ 1/
      data cf(210)/    0.102948388187D+04/                           
      data i1(210)/ 2/,i2(210)/ 0/,i3(210)/ 0/
      data i4(210)/ 1/,i5(210)/ 0/,i6(210)/ 2/
      data cf(211)/    0.866189861908D+03/                           
      data i1(211)/ 2/,i2(211)/ 0/,i3(211)/ 0/
      data i4(211)/ 1/,i5(211)/ 1/,i6(211)/ 1/
      data cf(212)/   -0.319833290751D+03/                           
      data i1(212)/ 2/,i2(212)/ 0/,i3(212)/ 0/
      data i4(212)/ 1/,i5(212)/ 2/,i6(212)/ 0/
      data cf(213)/   -0.564190763411D+03/                           
      data i1(213)/ 2/,i2(213)/ 0/,i3(213)/ 0/
      data i4(213)/ 2/,i5(213)/ 0/,i6(213)/ 1/
      data cf(214)/   -0.242487181986D+03/                           
      data i1(214)/ 2/,i2(214)/ 0/,i3(214)/ 0/
      data i4(214)/ 2/,i5(214)/ 1/,i6(214)/ 0/
      data cf(215)/   -0.296207075527D+02/                           
      data i1(215)/ 2/,i2(215)/ 0/,i3(215)/ 1/
      data i4(215)/ 0/,i5(215)/ 0/,i6(215)/ 2/
      data cf(216)/    0.587680804951D+03/                           
      data i1(216)/ 2/,i2(216)/ 0/,i3(216)/ 1/
      data i4(216)/ 0/,i5(216)/ 1/,i6(216)/ 1/
      data cf(217)/   -0.193248365814D+03/                           
      data i1(217)/ 2/,i2(217)/ 0/,i3(217)/ 1/
      data i4(217)/ 1/,i5(217)/ 0/,i6(217)/ 1/
      data cf(218)/   -0.379517171270D+03/                           
      data i1(218)/ 2/,i2(218)/ 0/,i3(218)/ 1/
      data i4(218)/ 1/,i5(218)/ 1/,i6(218)/ 0/
      data cf(219)/    0.239263179809D+02/                           
      data i1(219)/ 2/,i2(219)/ 0/,i3(219)/ 1/
      data i4(219)/ 2/,i5(219)/ 0/,i6(219)/ 0/
      data cf(220)/   -0.174655330523D+02/                           
      data i1(220)/ 2/,i2(220)/ 0/,i3(220)/ 2/
      data i4(220)/ 0/,i5(220)/ 0/,i6(220)/ 1/
      data cf(221)/    0.972920760651D+01/                           
      data i1(221)/ 2/,i2(221)/ 0/,i3(221)/ 2/
      data i4(221)/ 1/,i5(221)/ 0/,i6(221)/ 0/
      data cf(222)/    0.304621616143D+00/                           
      data i1(222)/ 2/,i2(222)/ 1/,i3(222)/ 0/
      data i4(222)/ 0/,i5(222)/ 0/,i6(222)/ 2/
      data cf(223)/   -0.335512088414D+03/                           
      data i1(223)/ 2/,i2(223)/ 1/,i3(223)/ 0/
      data i4(223)/ 0/,i5(223)/ 1/,i6(223)/ 1/
      data cf(224)/   -0.421004985175D+02/                           
      data i1(224)/ 2/,i2(224)/ 1/,i3(224)/ 0/
      data i4(224)/ 0/,i5(224)/ 2/,i6(224)/ 0/
      data cf(225)/    0.644175892667D+02/                           
      data i1(225)/ 2/,i2(225)/ 1/,i3(225)/ 0/
      data i4(225)/ 1/,i5(225)/ 0/,i6(225)/ 1/
      data cf(226)/    0.165903812759D+03/                           
      data i1(226)/ 2/,i2(226)/ 1/,i3(226)/ 0/
      data i4(226)/ 1/,i5(226)/ 1/,i6(226)/ 0/
      data cf(227)/    0.121005994015D+02/                           
      data i1(227)/ 2/,i2(227)/ 1/,i3(227)/ 1/
      data i4(227)/ 0/,i5(227)/ 0/,i6(227)/ 1/
      data cf(228)/   -0.303101104598D+01/                           
      data i1(228)/ 2/,i2(228)/ 1/,i3(228)/ 1/
      data i4(228)/ 0/,i5(228)/ 1/,i6(228)/ 0/
      data cf(229)/   -0.619470477363D+01/                           
      data i1(229)/ 2/,i2(229)/ 1/,i3(229)/ 1/
      data i4(229)/ 1/,i5(229)/ 0/,i6(229)/ 0/
      data cf(230)/    0.369880456480D+00/                           
      data i1(230)/ 2/,i2(230)/ 1/,i3(230)/ 2/
      data i4(230)/ 0/,i5(230)/ 0/,i6(230)/ 0/
      data cf(231)/   -0.275354963426D+01/                           
      data i1(231)/ 2/,i2(231)/ 2/,i3(231)/ 0/
      data i4(231)/ 0/,i5(231)/ 0/,i6(231)/ 1/
      data cf(232)/    0.627953939226D+01/                           
      data i1(232)/ 2/,i2(232)/ 2/,i3(232)/ 0/
      data i4(232)/ 0/,i5(232)/ 1/,i6(232)/ 0/
      data cf(233)/   -0.100087490686D+00/                           
      data i1(233)/ 2/,i2(233)/ 2/,i3(233)/ 1/
      data i4(233)/ 0/,i5(233)/ 0/,i6(233)/ 0/
      data cf(234)/    0.118106910390D+04/                           
      data i1(234)/ 3/,i2(234)/ 0/,i3(234)/ 0/
      data i4(234)/ 0/,i5(234)/ 1/,i6(234)/ 1/
      data cf(235)/   -0.946818406502D+03/                           
      data i1(235)/ 3/,i2(235)/ 0/,i3(235)/ 0/
      data i4(235)/ 1/,i5(235)/ 0/,i6(235)/ 1/
      data cf(236)/   -0.314099976753D+03/                           
      data i1(236)/ 3/,i2(236)/ 0/,i3(236)/ 0/
      data i4(236)/ 1/,i5(236)/ 1/,i6(236)/ 0/
      data cf(237)/   -0.214107591756D+01/                           
      data i1(237)/ 3/,i2(237)/ 0/,i3(237)/ 1/
      data i4(237)/ 0/,i5(237)/ 0/,i6(237)/ 1/
      data cf(238)/    0.372078669606D+02/                           
      data i1(238)/ 3/,i2(238)/ 0/,i3(238)/ 1/
      data i4(238)/ 1/,i5(238)/ 0/,i6(238)/ 0/
      data cf(239)/   -0.750908536548D+01/                           
      data i1(239)/ 3/,i2(239)/ 1/,i3(239)/ 0/
      data i4(239)/ 0/,i5(239)/ 0/,i6(239)/ 1/
      data cf(240)/   -0.873202217951D+01/                           
      data i1(240)/ 3/,i2(240)/ 1/,i3(240)/ 0/
      data i4(240)/ 0/,i5(240)/ 1/,i6(240)/ 0/
      data cf(241)/    0.261533605386D+00/                           
      data i1(241)/ 3/,i2(241)/ 1/,i3(241)/ 1/
      data i4(241)/ 0/,i5(241)/ 0/,i6(241)/ 0/
      ener=0.d0          
      pux12= vex6(1)*r12 
      pux13= vex6(2)*r13 
      pux14= vex6(3)*r14 
      pux23= vex6(4)*r23 
      pux24= vex6(5)*r24 
      pux34= vex6(6)*r34 
      qux12= dexp(-pux12) 
      qux13= dexp(-pux13) 
      qux14= dexp(-pux14) 
      qux23= dexp(-pux23) 
      qux24= dexp(-pux24) 
      qux34= dexp(-pux34) 
      aux12= r12*qux12 
      aux13= r13*qux13 
      aux14= r14*qux14 
      aux23= r23*qux23 
      aux24= r24*qux24 
      aux34= r34*qux34 
      do 1 i=1,mt_2               
         f12(i) = aux12*f12(i-1)  
         f13(i) = aux13*f13(i-1)  
         f14(i) = aux14*f14(i-1)  
         f23(i) = aux23*f23(i-1)  
         f24(i) = aux24*f24(i-1)  
         f34(i) = aux34*f34(i-1)  
1     continue                 
      do 2 j=1,is                               
         f1213=f12(i1(j))*f13(i2(j))            
         f1423=f14(i3(j))*f23(i4(j))            
         f2434=f24(i5(j))*f34(i6(j))            
         au=f1213*f1423*f2434                
         ener=ener+cf(j)*au                     
         if (iop.eq.1) then                     
            f14232434=f1423*f2434               
            f12132434=f1213*f2434               
            f12131423=f1213*f1423               
            d12=i1(j)*f12(i1(j)-1)*f13(i2(j))*f14232434  
            d13=i2(j)*f12(i1(j))*f13(i2(j)-1)*f14232434  
            d14=i3(j)*f14(i3(j)-1)*f23(i4(j))*f12132434  
            d23=i4(j)*f14(i3(j))*f23(i4(j)-1)*f12132434  
            d24=i5(j)*f12131423*f24(i5(j)-1)*f34(i6(j))  
            d34=i6(j)*f12131423*f24(i5(j))*f34(i6(j)-1)  
            der12=der12+cf(j)*d12                          
            der13=der13+cf(j)*d13                          
            der14=der14+cf(j)*d14                          
            der23=der23+cf(j)*d23                          
            der24=der24+cf(j)*d24                          
            der34=der34+cf(j)*d34                          
         endif                                             
2     continue                                              
      if (iop.eq.1) then                 
         der(1)=der12*(1.d0-pux12)*qux12    
         der(2)=der13*(1.d0-pux13)*qux13    
         der(3)=der14*(1.d0-pux14)*qux14    
         der(4)=der23*(1.d0-pux23)*qux23    
         der(5)=der24*(1.d0-pux24)*qux24    
         der(6)=der34*(1.d0-pux34)*qux34    
      endif                                        
      return                                                
      end                                                   
