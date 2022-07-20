      real*8 function Energy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop,      
     5       ijatom,
     6       M  ) 
       integer L, L_p,nnnumber, nnnum_p
       integer idx(nnnum_p), jdy(nnnum_p), kdz(nnnum_p)
       real*8  D, Jcoop, E2, E1
       integer iE1, iE2, inE2
       integer ijatom(0:L_p+1,0:L_p+1,0:L_p+1)
       integer M
       integer i,j,k,o,i2,j2,k2
       
       
       iE1=0       
         do 21 i=1,L
          do 22 j=1,L     
	   do 23 k=1,L
           iE1=iE1+ijatom(i,j,k)
23	   continue
22        continue
21       continue  
       M=iE1       
       E1=D*iE1
       
       iE2=0
       do 31 i=1,L
        do 32 j=1,L
	 do 33 k=1,L
         inE2=0 
         do 34 o=1,nnnumber
           i2=i+idx(o)
           j2=j+jdy(o)
	   k2=k+kdz(o)
           inE2=inE2+ijatom(i2,j2,k2)
34       continue           
         iE2=iE2+ijatom(i,j,k)*inE2
33	 continue 
32      continue
31     continue  
       E2=Jcoop*iE2
       Energy=E1+0.5*E2
      end
