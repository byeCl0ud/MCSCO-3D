      real*8 function dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,j,k  ) 
       integer L, L_p, nnnumber, nnnum_p
       integer idx(nnnum_p), jdy(nnnum_p), kdz(nnnum_p) 
       real*8  D, Jcoop, r, T, dE2, dE1
       integer idE2, inE2
       integer ijatom(0:L_p+1,0:L_p+1,0:L_p+1)
       integer i,j,k       
       integer o,i2,j2,k2
       
        dE1=-2*ijatom(i,j,k)*D+2*ijatom(i,j,k)*r*T
       
         inE2=0 
         do 33 o=1,nnnumber
           i2=i+idx(o)
           j2=j+jdy(o)
	   k2=k+kdz(o)
           inE2=inE2+ijatom(i2,j2,k2)
33       continue           
         idE2=-2*ijatom(i,j,k)*inE2   
        dE2=Jcoop*idE2
       dEnergy=dE1+dE2

      end
