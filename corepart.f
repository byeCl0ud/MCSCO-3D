        subroutine corepart(
     1          outputFile,
     2          nMC,
     3          L, Ti, Tf, nTs, nrp,
     4          D, Jcoop, r, q,
     5          ConcInit )

         character output*128
         integer L, nTs, q, nrp
         integer iTstep
         real*8 D, Jcoop, r, Ti, Tf
         real*8 dT, ConcInit

         integer L_p, nTs_p, lt_p, nrp_p
       parameter( 
     1   L_p=128,
     2   nTs_p=500,
     3   lt_p=8,
     4   nrp_p=2000)   
       
       real*8  Einst( nTs_p, nrp_p )
       integer Minst( nTs_p, nrp_p )
       real*8  rMinst( nrp_p )
       real*8  Et( nrp_p )
       integer Mt( nrp_p )
       real*8  Eav( nTs_p ), dEav( nTs_p )
       real*8  E2av( nTs_p )
       real*8  avM( nTs_p ), dMav( nTs_p ) 
       real*8  Eavg, avgM, dEavg, avgdM, E2avg

c Hysteresis-leading variables
       real*8 ni, ntt, leadT
c  Internal variables:
       integer itime, M, i,j,k
       real*8  T, E
c nncoords:
       integer nnnumber, nnnum_p
       parameter( nnnum_p=8 )
       integer idx(nnnum_p), jdy(nnnum_p), kdz(nnnum_p)
c latinit:
       integer ijatom(0:L_p+1,0:L_p+1,0:L_p+1)
c functions:
       real*8  Energy, dEnergy

       call lattice(
     2       nnnumber, nnnum_p,
     3       idx, jdy, kdz )
       
       call latinit(
     1       L, L_p,
     2       ConcInit,
     3       ijatom )

      WRITE(*,*) 'Primary lattice initialization O.K.'
       
       E = Energy(
     1       L, L_p,
     2       nnnumber, nnnum_p,
     3       idx, jdy, kdz,
     4       D, Jcoop,
     5       ijatom,
     6       M )  

       ni=nint(float(nTs-q)/float(q))
       dT=(Tf-Ti)/ni
       T=Ti-dT
       
       if (q.eq.1) then
          ntt=1
       else
          ntt= 2 - nint(float(nTs)/2)*2 + nTs
       endif

       leadT = dabs((Tf-Ti)+dT+dT*ntt/2)

       open(UNIT=31, FILE='outputFile', status='UNKNOWN')
       
       WRITE(*,*) 'E=',E,' D=',D,' J=',Jcoop,' r=',r
       WRITE(*,*) 'Here we go, running the calculation...'
       WRITE(*,*) 'To preview progress of the calculation,'
       WRITE(*,*) 'run "tail -f output" in another terminal'
       
      do 100 iTstep=1, nTs
       leadT = leadT - dabs(dT)
       T = T + dT*(nint(leadT/dabs(leadT-1.d-3)))**(3-q)

       call MonteCarlo(
     1        L, L_p,
     2        nnnumber, nnnum_p,
     3        idx, jdy, kdz,
     4        D, Jcoop, r, T,
     5        E, M, 
     6        Et, Mt,
     7        ijatom,
     8        nMC, nrp, nrp_p)

       Eavg=0
       E2avg=0
       avgM=0
       do 111 itime=1,nrp
         Einst(iTstep,itime)=Et(itime)
         Minst(iTstep,itime)=Mt(itime) 
         rMinst(itime)=dfloat(Mt(itime))/dfloat(L*L*L)       
         Eavg=Eavg+Et(itime)
         E2avg=E2avg+Et(itime)*Et(itime)
         avgM=avgM+rMinst(itime)
111    continue 
       Eavg = Eavg/dfloat(nrp)
       E2avg = E2avg/dfloat(nrp)
       avgM = avgM/dfloat(nrp)
       avM( iTstep ) = avgM         
       
       dEavg=0
       avgdM=0
       do 112 itime=1,nrp
         dEavg = dEavg + (Et(itime)-Eavg)**2
         avgdM = avgdM + (rMinst(itime)-avgM)**2
112    continue 
       dEavg = dEavg/dfloat(nrp*(nrp-1))
       avgdM = avgdM/dfloat(nrp*(nrp-1))
       dEavg = sqrt( dEavg ) 
       avgdM = sqrt( avgdM )
       dMav( iTstep ) = avgdM
       Eav( iTstep ) = Eavg/dfloat(L*L*L)  
       E2av( iTstep ) = E2avg/dfloat(L*L*L*L)    
       dEav( iTstep ) = dEavg/dfloat(L*L*L) 
       
       write(31,1131) T, Eav(iTstep), avM( iTstep )
1131   format(1X,F15.9,';',F18.9,'; ',F12.9)
           
           
100   continue       

      close(35)

      return
      end
