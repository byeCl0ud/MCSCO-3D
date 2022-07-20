        subroutine MonteCarlo(
     1        L, L_p,
     2        nnnumber, nnnum_p,
     3        idx, jdy, kdz,
     4        D, Jcoop, r, T,
     5        E, M, 
     6        Et, Mt,
     7        ijatom,
     8        nMC, nrp, nrp_p)

       integer L, L_p, nnnumber, nnnum_p
       integer idx(nnnum_p), jdy(nnnum_p), kdz(nnnum_p)       
       real*8  D, Jcoop, r, T, E, dE, Et( nrp_p )
       integer M, Mt( nrp_p )
       integer ijatom(0:L_p+1,0:L_p+1,0:L_p+1)
       integer nMC, nrp, nrp_p 
       integer iMC, itime
      
       integer i,j,k
       real*8  beta
       real*8  random, dEnergy
       integer nreject,nacrais,naclower
       
       nattemp=0
       nreject=0
       nacrais=0
       naclower=0
       beta=1.0/T
       
       do 101 itime=1,nrp
        do 102 iMC=1, nMC 

cc   Periodic Boundary: ijatom
c    (1,1,1)=(1,1,L+1)=(1,L+1,1)=(L+1,1,1)=
c    =(1,L+1,L+1)=(L+1,1,L+1)=(L+1,L+1,1)=
c    =(L+1,L+1,L+1)

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,1,1  )   
         if(dE .le. 0)then
             ijatom(1,1,1)=-ijatom(1,1,1)
	     ijatom(1,1,L+1)=-ijatom(1,1,L+1)
	     ijatom(1,L+1,1)=-ijatom(1,L+1,1)
	     ijatom(L+1,1,1)=-ijatom(L+1,1,1)
	     ijatom(1,L+1,L+1)=-ijatom(1,L+1,L+1)
	     ijatom(L+1,1,L+1)=-ijatom(L+1,1,L+1)
	     ijatom(L+1,L+1,1)=-ijatom(L+1,L+1,1)
	     ijatom(L+1,L+1,L+1)=-ijatom(L+1,L+1,L+1)
             M=M+2*ijatom(1,1,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
             ijatom(1,1,1)=-ijatom(1,1,1)
	     ijatom(1,1,L+1)=-ijatom(1,1,L+1)
	     ijatom(1,L+1,1)=-ijatom(1,L+1,1)
	     ijatom(L+1,1,1)=-ijatom(L+1,1,1)
	     ijatom(1,L+1,L+1)=-ijatom(1,L+1,L+1)
	     ijatom(L+1,1,L+1)=-ijatom(L+1,1,L+1)
	     ijatom(L+1,L+1,1)=-ijatom(L+1,L+1,1)
	     ijatom(L+1,L+1,L+1)=-ijatom(L+1,L+1,L+1)
             M=M+2*ijatom(1,1,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif         

c    (1,1,L)=(1,1,0)=(1,L+1,0)=(L+1,1,0)=
c    (L+1,L+1,0)=(1,L+1,L)=(L+1,1,L)=
c    (L+1,L+1,L)

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,1,L  )   
         if(dE .le. 0)then
	     ijatom(1,1,L)=-ijatom(1,1,L)
	     ijatom(1,1,0)=-ijatom(1,1,0)
	     ijatom(1,L+1,0)=-ijatom(1,L+1,0)
	     ijatom(L+1,1,0)=-ijatom(L+1,1,0)
	     ijatom(L+1,L+1,0)=-ijatom(L+1,L+1,0)
	     ijatom(1,L+1,L)=-ijatom(1,L+1,L)
	     ijatom(L+1,1,L)=-ijatom(L+1,1,L)
	     ijatom(L+1,L+1,L)=-ijatom(L+1,L+1,L)
             M=M+2*ijatom(1,1,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,1,L)=-ijatom(1,1,L)
	     ijatom(1,1,0)=-ijatom(1,1,0)
	     ijatom(1,L+1,0)=-ijatom(1,L+1,0)
	     ijatom(L+1,1,0)=-ijatom(L+1,1,0)
	     ijatom(L+1,L+1,0)=-ijatom(L+1,L+1,0)
	     ijatom(1,L+1,L)=-ijatom(1,L+1,L)
	     ijatom(L+1,1,L)=-ijatom(L+1,1,L)
	     ijatom(L+1,L+1,L)=-ijatom(L+1,L+1,L)
             M=M+2*ijatom(1,1,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (1,L,1)=(1,0,1)=(L+1,0,1)=(1,0,L+1)=
c    (L+1,0,L+1)=(L+1,L,1)=(1,L,L+1)=
c    (L+1,L,L+1)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,L,1  )   
         if(dE .le. 0)then
	     ijatom(1,L,1)=-ijatom(1,L,1)
	     ijatom(1,0,1)=-ijatom(1,0,1)
	     ijatom(L+1,0,1)=-ijatom(L+1,0,1)
	     ijatom(1,0,L+1)=-ijatom(1,0,L+1)
	     ijatom(L+1,0,L+1)=-ijatom(L+1,0,L+1)
	     ijatom(L+1,L,1)=-ijatom(L+1,L,1)
	     ijatom(1,L,L+1)=-ijatom(1,L,L+1)
	     ijatom(L+1,L,L+1)=-ijatom(L+1,L,L+1)
             M=M+2*ijatom(1,L,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,L,1)=-ijatom(1,L,1)
	     ijatom(1,0,1)=-ijatom(1,0,1)
	     ijatom(L+1,0,1)=-ijatom(L+1,0,1)
	     ijatom(1,0,L+1)=-ijatom(1,0,L+1)
	     ijatom(L+1,0,L+1)=-ijatom(L+1,0,L+1)
	     ijatom(L+1,L,1)=-ijatom(L+1,L,1)
	     ijatom(1,L,L+1)=-ijatom(1,L,L+1)
	     ijatom(L+1,L,L+1)=-ijatom(L+1,L,L+1)
             M=M+2*ijatom(1,L,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (L,1,1)=(0,1,1)=(L,L+1,1)=(L,1,L+1)=
c    (L,L+1,L+1)=(0,L+1,1)=(0,1,L+1)=
c    (0,L+1,L+1)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,1,1  )   
         if(dE .le. 0)then
	     ijatom(L,1,1)=-ijatom(L,1,1)
	     ijatom(0,1,1)=-ijatom(0,1,1)
	     ijatom(L,L+1,1)=-ijatom(L,L+1,1)
	     ijatom(L,1,L+1)=-ijatom(L,1,L+1)
	     ijatom(L,L+1,L+1)=-ijatom(L,L+1,L+1)
	     ijatom(0,L+1,1)=-ijatom(0,L+1,1)
	     ijatom(0,1,L+1)=-ijatom(0,1,L+1)
	     ijatom(0,L+1,L+1)=-ijatom(0,L+1,L+1)
             M=M+2*ijatom(L,1,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(L,1,1)=-ijatom(L,1,1)
	     ijatom(0,1,1)=-ijatom(0,1,1)
	     ijatom(L,L+1,1)=-ijatom(L,L+1,1)
	     ijatom(L,1,L+1)=-ijatom(L,1,L+1)
	     ijatom(L,L+1,L+1)=-ijatom(L,L+1,L+1)
	     ijatom(0,L+1,1)=-ijatom(0,L+1,1)
	     ijatom(0,1,L+1)=-ijatom(0,1,L+1)
	     ijatom(0,L+1,L+1)=-ijatom(0,L+1,L+1)
             M=M+2*ijatom(L,1,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (1,L,L)=(1,0,0)=(L+1,L,L)=(L+1,0,0)=
c    (1,L,0)=(1,0,L)=(L+1,L,0)=(L+1,0,L)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,L,L  )   
         if(dE .le. 0)then
	     ijatom(1,L,L)=-ijatom(1,L,L)
	     ijatom(1,0,0)=-ijatom(1,0,0)
	     ijatom(L+1,L,L)=-ijatom(L+1,L,L)
	     ijatom(L+1,0,0)=-ijatom(L+1,0,0)
	     ijatom(1,L,0)=-ijatom(1,L,0)
	     ijatom(1,0,L)=-ijatom(1,0,L)
	     ijatom(L+1,L,0)=-ijatom(L+1,L,0)
	     ijatom(L+1,0,L)=-ijatom(L+1,0,L)
             M=M+2*ijatom(1,L,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,L,L)=-ijatom(1,L,L)
	     ijatom(1,0,0)=-ijatom(1,0,0)
	     ijatom(L+1,L,L)=-ijatom(L+1,L,L)
	     ijatom(L+1,0,0)=-ijatom(L+1,0,0)
	     ijatom(1,L,0)=-ijatom(1,L,0)
	     ijatom(1,0,L)=-ijatom(1,0,L)
	     ijatom(L+1,L,0)=-ijatom(L+1,L,0)
	     ijatom(L+1,0,L)=-ijatom(L+1,0,L)
             M=M+2*ijatom(1,L,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (L,1,L)=(0,1,0)=(L,L+1,L)=(0,L+1,0)=
c    (L,1,0)=(0,1,L)=(L,L+1,0)=(0,L+1,L)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,1,L  )   
         if(dE .le. 0)then
	     ijatom(L,1,L)=-ijatom(L,1,L)
	     ijatom(0,1,0)=-ijatom(0,1,0)
	     ijatom(L,L+1,L)=-ijatom(L,L+1,L)
	     ijatom(0,L+1,0)=-ijatom(0,L+1,0)
	     ijatom(L,1,0)=-ijatom(L,1,0)
	     ijatom(0,1,L)=-ijatom(0,1,L)
	     ijatom(L,L+1,0)=-ijatom(L,L+1,0)
	     ijatom(0,L+1,L)=-ijatom(0,L+1,L)
             M=M+2*ijatom(L,1,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(L,1,L)=-ijatom(L,1,L)
	     ijatom(0,1,0)=-ijatom(0,1,0)
	     ijatom(L,L+1,L)=-ijatom(L,L+1,L)
	     ijatom(0,L+1,0)=-ijatom(0,L+1,0)
	     ijatom(L,1,0)=-ijatom(L,1,0)
	     ijatom(0,1,L)=-ijatom(0,1,L)
	     ijatom(L,L+1,0)=-ijatom(L,L+1,0)
	     ijatom(0,L+1,L)=-ijatom(0,L+1,L)
             M=M+2*ijatom(L,1,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (1,L,L)=(1,0,0)=(L+1,L,L)=(L+1,0,0)=
c    (1,L,0)=(1,0,L)=(L+1,L,0)=(L+1,0,L)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,L,L  )   
         if(dE .le. 0)then
	     ijatom(1,L,L)=-ijatom(1,L,L)
	     ijatom(1,0,0)=-ijatom(1,0,0)
	     ijatom(L+1,L,L)=-ijatom(L+1,L,L)
	     ijatom(L+1,0,0)=-ijatom(L+1,0,0)
	     ijatom(1,L,0)=-ijatom(1,L,0)
	     ijatom(1,0,L)=-ijatom(1,0,L)
	     ijatom(L+1,L,0)=-ijatom(L+1,L,0)
	     ijatom(L+1,0,L)=-ijatom(L+1,0,L)
             M=M+2*ijatom(1,L,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,L,L)=-ijatom(1,L,L)
	     ijatom(1,0,0)=-ijatom(1,0,0)
	     ijatom(L+1,L,L)=-ijatom(L+1,L,L)
	     ijatom(L+1,0,0)=-ijatom(L+1,0,0)
	     ijatom(1,L,0)=-ijatom(1,L,0)
	     ijatom(1,0,L)=-ijatom(1,0,L)
	     ijatom(L+1,L,0)=-ijatom(L+1,L,0)
	     ijatom(L+1,0,L)=-ijatom(L+1,0,L)
             M=M+2*ijatom(1,L,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (L,L,L)=(0,0,0)=(L,L,0)=(L,0,L)=
c    (0,L,L)=(0,0,L)=(L,0,0)=(0,L,0)
         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,L,L  )   
         if(dE .le. 0)then
	     ijatom(L,L,L)=-ijatom(L,L,L)
	     ijatom(0,0,0)=-ijatom(0,0,0)
	     ijatom(L,L,0)=-ijatom(L,L,0)
	     ijatom(L,0,L)=-ijatom(L,0,L)
	     ijatom(0,L,L)=-ijatom(0,L,L)
	     ijatom(0,0,L)=-ijatom(0,0,L)
	     ijatom(L,0,0)=-ijatom(L,0,0)
	     ijatom(0,L,0)=-ijatom(0,L,0)
             M=M+2*ijatom(L,L,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(L,L,L)=-ijatom(L,L,L)
	     ijatom(0,0,0)=-ijatom(0,0,0)
	     ijatom(L,L,0)=-ijatom(L,L,0)
	     ijatom(L,0,L)=-ijatom(L,0,L)
	     ijatom(0,L,L)=-ijatom(0,L,L)
	     ijatom(0,0,L)=-ijatom(0,0,L)
	     ijatom(L,0,0)=-ijatom(L,0,0)
	     ijatom(0,L,0)=-ijatom(0,L,0)
             M=M+2*ijatom(L,L,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

c    (1,1,k)=(L+1,1,k)=(1,L+1,k)=(L+1,L+1,k)
       do 60 k=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,1,k  )   
         if(dE .le. 0)then
	     ijatom(1,1,k)=-ijatom(1,1,k)
	     ijatom(L+1,1,k)=-ijatom(L+1,1,k)
	     ijatom(1,L+1,k)=-ijatom(1,L+1,k)
	     ijatom(L+1,L+1,k)=-ijatom(L+1,L+1,k)             
             M=M+2*ijatom(1,1,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,1,k)=-ijatom(1,1,k)
	     ijatom(L+1,1,k)=-ijatom(L+1,1,k)
	     ijatom(1,L+1,k)=-ijatom(1,L+1,k)
	     ijatom(L+1,L+1,k)=-ijatom(L+1,L+1,k)            
             M=M+2*ijatom(1,1,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

60     continue

c    (1,L,k)=(L+1,L,k)=(1,0,k)=(L+1,0,k)
       do 61 k=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,L,k  )   
         if(dE .le. 0)then
	     ijatom(1,L,k)=-ijatom(1,L,k)
	     ijatom(L+1,L,k)=-ijatom(L+1,L,k)
	     ijatom(1,0,k)=-ijatom(1,0,k)
	     ijatom(L+1,0,k)=-ijatom(L+1,0,k)             
             M=M+2*ijatom(1,L,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,L,k)=-ijatom(1,L,k)
	     ijatom(L+1,L,k)=-ijatom(L+1,L,k)
	     ijatom(1,0,k)=-ijatom(1,0,k)
	     ijatom(L+1,0,k)=-ijatom(L+1,0,k)            
             M=M+2*ijatom(1,L,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

61     continue

c    (L,1,k)=(L,L+1,k)=(0,1,k)=(0,L+1,k)
       do 62 k=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,1,k  )   
         if(dE .le. 0)then
  	     ijatom(L,1,k)=-ijatom(L,1,k)
  	     ijatom(L,L+1,k)=-ijatom(L,L+1,k)
  	     ijatom(0,1,k)=-ijatom(0,1,k)
  	     ijatom(0,L+1,k)=-ijatom(0,L+1,k)            
             M=M+2*ijatom(L,1,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
  	     ijatom(L,1,k)=-ijatom(L,1,k)
  	     ijatom(L,L+1,k)=-ijatom(L,L+1,k)
  	     ijatom(0,1,k)=-ijatom(0,1,k)
  	     ijatom(0,L+1,k)=-ijatom(0,L+1,k)            
             M=M+2*ijatom(L,1,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

62     continue

c    (L,L,k)=(L,0,k)=(0,L,k)=(0,0,k)
       do 63 k=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,L,k  )   
         if(dE .le. 0)then
 	     ijatom(L,L,k)=-ijatom(L,L,k)
 	     ijatom(L,0,k)=-ijatom(L,0,k)
 	     ijatom(0,L,k)=-ijatom(0,L,k)
 	     ijatom(0,0,k)=-ijatom(0,0,k)            
             M=M+2*ijatom(L,L,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
 	     ijatom(L,L,k)=-ijatom(L,L,k)
 	     ijatom(L,0,k)=-ijatom(L,0,k)
 	     ijatom(0,L,k)=-ijatom(0,L,k)
 	     ijatom(0,0,k)=-ijatom(0,0,k)                 
             M=M+2*ijatom(L,L,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

63     continue

c    (1,j,1)=(1,j,L+1)=(L+1,j,1)=(L+1,j,L+1)
       do 70 j=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,j,1  )   
         if(dE .le. 0)then
             ijatom(1,j,1)=-ijatom(1,j,1)
  	     ijatom(1,j,L+1)=-ijatom(1,j,L+1)
  	     ijatom(L+1,j,1)=-ijatom(L+1,j,1)
  	     ijatom(L+1,j,L+1)=-ijatom(L+1,j,L+1)            
             M=M+2*ijatom(1,j,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
             ijatom(1,j,1)=-ijatom(1,j,1)
  	     ijatom(1,j,L+1)=-ijatom(1,j,L+1)
  	     ijatom(L+1,j,1)=-ijatom(L+1,j,1)
  	     ijatom(L+1,j,L+1)=-ijatom(L+1,j,L+1)               
             M=M+2*ijatom(1,j,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

70     continue

c    (1,j,L)=(1,j,0)=(L+1,j,L)=(L+1,j,0)
       do 71 j=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,j,L  )   
         if(dE .le. 0)then
    	     ijatom(1,j,L)=-ijatom(1,j,L)
    	     ijatom(1,j,0)=-ijatom(1,j,0)
    	     ijatom(L+1,j,L)=-ijatom(L+1,j,L)
    	     ijatom(L+1,j,0)=-ijatom(L+1,j,0)           
             M=M+2*ijatom(1,j,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
    	     ijatom(1,j,L)=-ijatom(1,j,L)
    	     ijatom(1,j,0)=-ijatom(1,j,0)
    	     ijatom(L+1,j,L)=-ijatom(L+1,j,L)
    	     ijatom(L+1,j,0)=-ijatom(L+1,j,0)             
             M=M+2*ijatom(1,j,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

71     continue

c    (L,j,1)=(L,j,L+1)=(0,j,1)=(0,j,L+1)
       do 72 j=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,j,1  )   
         if(dE .le. 0)then
     	     ijatom(L,j,1)=-ijatom(L,j,1)
    	     ijatom(L,j,L+1)=-ijatom(L,j,L+1)
     	     ijatom(0,j,1)=-ijatom(0,j,1)
     	     ijatom(0,j,L+1)=-ijatom(0,j,L+1)          
             M=M+2*ijatom(L,j,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
     	     ijatom(L,j,1)=-ijatom(L,j,1)
    	     ijatom(L,j,L+1)=-ijatom(L,j,L+1)
     	     ijatom(0,j,1)=-ijatom(0,j,1)
     	     ijatom(0,j,L+1)=-ijatom(0,j,L+1)            
             M=M+2*ijatom(L,j,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

72     continue

c    (L,j,L)=(0,j,L)=(L,j,0)=(0,j,0)
       do 73 j=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,j,L  )   
         if(dE .le. 0)then
             ijatom(L,j,L)=-ijatom(L,j,L)
             ijatom(0,j,L)=-ijatom(0,j,L)
             ijatom(L,j,0)=-ijatom(L,j,0)
             ijatom(0,j,0)=-ijatom(0,j,0)        
             M=M+2*ijatom(L,j,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
             ijatom(L,j,L)=-ijatom(L,j,L)
             ijatom(0,j,L)=-ijatom(0,j,L)
             ijatom(L,j,0)=-ijatom(L,j,0)
             ijatom(0,j,0)=-ijatom(0,j,0)           
             M=M+2*ijatom(L,j,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

73     continue

c    (i,1,1)=(i,1,L+1)=(i,L+1,1)=(i,L+1,L+1)
       do 80 i=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,1,1  )   
         if(dE .le. 0)then
	     ijatom(i,1,1)=-ijatom(i,1,1)
	     ijatom(i,1,L+1)=-ijatom(i,1,L+1)
	     ijatom(i,L+1,1)=-ijatom(i,L+1,1)
	     ijatom(i,L+1,L+1)=-ijatom(i,L+1,L+1)        
             M=M+2*ijatom(i,1,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,1,1)=-ijatom(i,1,1)
	     ijatom(i,1,L+1)=-ijatom(i,1,L+1)
	     ijatom(i,L+1,1)=-ijatom(i,L+1,1)
	     ijatom(i,L+1,L+1)=-ijatom(i,L+1,L+1)           
             M=M+2*ijatom(i,1,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

80     continue

c    (i,1,L)=(i,L+1,L)=(i,1,0)=(i,L+1,0)
       do 81 i=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,1,L  )   
         if(dE .le. 0)then
	     ijatom(i,1,L)=-ijatom(i,1,L)
	     ijatom(i,L+1,L)=-ijatom(i,L+1,L)
	     ijatom(i,1,0)=-ijatom(i,1,0)
	     ijatom(i,L+1,0)=-ijatom(i,L+1,0)       
             M=M+2*ijatom(i,1,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,1,L)=-ijatom(i,1,L)
	     ijatom(i,L+1,L)=-ijatom(i,L+1,L)
	     ijatom(i,1,0)=-ijatom(i,1,0)
	     ijatom(i,L+1,0)=-ijatom(i,L+1,0)           
             M=M+2*ijatom(i,1,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

81     continue

c    (i,L,1)=(i,L,L+1)=(i,0,1)=(i,0,L+1)
       do 82 i=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,L,1  )   
         if(dE .le. 0)then
	     ijatom(i,L,1)=-ijatom(i,L,1)
             ijatom(i,L,L+1)=-ijatom(i,L,L+1)
	     ijatom(i,0,1)=-ijatom(i,0,1)
	     ijatom(i,0,L+1)=-ijatom(i,0,L+1)       
             M=M+2*ijatom(i,L,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,L,1)=-ijatom(i,L,1)
             ijatom(i,L,L+1)=-ijatom(i,L,L+1)
	     ijatom(i,0,1)=-ijatom(i,0,1)
	     ijatom(i,0,L+1)=-ijatom(i,0,L+1)           
             M=M+2*ijatom(i,L,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

82     continue

c    (i,L,L)=(i,L,0)=(i,0,L)=(i,0,0)
       do 83 i=2,L-1 

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,L,L  )   
         if(dE .le. 0)then
	     ijatom(i,L,L)=-ijatom(i,L,L)
 	     ijatom(i,L,0)=-ijatom(i,L,0)
	     ijatom(i,0,L)=-ijatom(i,0,L)
	     ijatom(i,0,0)=-ijatom(i,0,0)       
             M=M+2*ijatom(i,L,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,L,L)=-ijatom(i,L,L)
 	     ijatom(i,L,0)=-ijatom(i,L,0)
	     ijatom(i,0,L)=-ijatom(i,0,L)
	     ijatom(i,0,0)=-ijatom(i,0,0)           
             M=M+2*ijatom(i,L,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

83     continue

c    (i,j,1)=(i,j,L+1)
       do 600 i=2,L-1 
        do 601 j=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,j,1  )   
         if(dE .le. 0)then
	     ijatom(i,j,1)=-ijatom(i,j,1)
	     ijatom(i,j,L+1)=-ijatom(i,j,L+1)      
             M=M+2*ijatom(i,j,1) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,j,1)=-ijatom(i,j,1)
	     ijatom(i,j,L+1)=-ijatom(i,j,L+1)           
             M=M+2*ijatom(i,j,1) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
601     continue
600    continue

c    (i,j,L)=(i,j,0)
       do 602 i=2,L-1 
        do 603 j=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,j,L  )   
         if(dE .le. 0)then
	     ijatom(i,j,L)=-ijatom(i,j,L)
	     ijatom(i,j,0)=-ijatom(i,j,0)     
             M=M+2*ijatom(i,j,L) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,j,L)=-ijatom(i,j,L)
	     ijatom(i,j,0)=-ijatom(i,j,0)           
             M=M+2*ijatom(i,j,L) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
603     continue
602    continue

c    (i,1,k)=(i,L+1,k)
       do 700 i=2,L-1 
        do 701 k=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,1,k  )   
         if(dE .le. 0)then
             ijatom(i,1,k)=-ijatom(i,1,k)
	     ijatom(i,L+1,k)=-ijatom(i,L+1,k)  
             M=M+2*ijatom(i,1,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
             ijatom(i,1,k)=-ijatom(i,1,k)
	     ijatom(i,L+1,k)=-ijatom(i,L+1,k)           
             M=M+2*ijatom(i,1,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
701     continue
700    continue

c    (i,L,k)=(i,0,k)
       do 702 i=2,L-1 
        do 703 k=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,L,k  )   
         if(dE .le. 0)then
	     ijatom(i,L,k)=-ijatom(i,L,k)
	     ijatom(i,0,k)=-ijatom(i,0,k)
             M=M+2*ijatom(i,L,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,L,k)=-ijatom(i,L,k)
	     ijatom(i,0,k)=-ijatom(i,0,k)           
             M=M+2*ijatom(i,L,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
703     continue
702    continue

c    (1,j,k)=(L+1,j,k)
       do 800 j=2,L-1 
        do 801 k=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,1,j,k  )   
         if(dE .le. 0)then
	     ijatom(1,j,k)=-ijatom(1,j,k)
	     ijatom(L+1,j,k)=-ijatom(L+1,j,k)
             M=M+2*ijatom(1,j,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(1,j,k)=-ijatom(1,j,k)
	     ijatom(L+1,j,k)=-ijatom(L+1,j,k)           
             M=M+2*ijatom(1,j,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
801     continue
800    continue

c    (L,j,k)=(0,j,k)
       do 802 j=2,L-1 
        do 803 k=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,L,j,k  )   
         if(dE .le. 0)then
	     ijatom(L,j,k)=-ijatom(L,j,k)
	     ijatom(0,j,k)=-ijatom(0,j,k)
             M=M+2*ijatom(L,j,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(L,j,k)=-ijatom(L,j,k)
	     ijatom(0,j,k)=-ijatom(0,j,k)           
             M=M+2*ijatom(L,j,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif
803     continue
802    continue

c    (i,j,k) has no borders
       do 900 i=2,L-1 
        do 901 j=2,L-1
	 do 902 k=2,L-1

         dE=dEnergy(
     1       L, L_p,       
     2       nnnumber, nnnum_p,       
     3       idx, jdy, kdz,
     4       D, Jcoop, r, T,    
     5       ijatom,i,j,k  )   
         if(dE .le. 0)then
	     ijatom(i,j,k)=-ijatom(i,j,k)
             M=M+2*ijatom(i,j,k) 
             E=E+dE
             naclower=naclower+1
         else if( random() .le. exp(-beta*dE) )then
	     ijatom(i,j,k)=-ijatom(i,j,k)           
             M=M+2*ijatom(i,j,k) 
             E=E+dE
             nacrais=nacrais+1 
         else     
             nreject=nreject+1     
         endif

902	 continue
901     continue
900    continue
                
102    continue

       Et(itime)=E
       Mt(itime)=M

101   continue 

      end
