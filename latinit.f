      subroutine latinit(
     1       L, L_p,
     2       ConcInit,
     1       ijatom  ) 
c
c    Initialises lattice and buffer layer.
c      ijatom(i,j) = 1 or -1
c
       implicit none 
       integer L, L_p
       real*8  ConcInit      
       integer ijatom(0:L_p+1,0:L_p+1,0:L_p+1)
       integer itype 
       integer i,j,k 
       integer nexttype
       
cc   Periodic Boundary: ijatom
c    (1,1,1)=(1,1,L+1)=(1,L+1,1)=(L+1,1,1)=
c    =(1,L+1,L+1)=(L+1,1,L+1)=(L+1,L+1,1)=
c    =(L+1,L+1,L+1)
       itype=nexttype(ConcInit)
       ijatom(1,1,1)=itype
       ijatom(1,1,L+1)=itype
       ijatom(1,L+1,1)=itype
       ijatom(L+1,1,1)=itype
       ijatom(1,L+1,L+1)=itype
       ijatom(L+1,1,L+1)=itype
       ijatom(L+1,L+1,1)=itype
       ijatom(L+1,L+1,L+1)=itype
c    (1,1,L)=(1,1,0)=(1,L+1,0)=(L+1,1,0)=
c    (L+1,L+1,0)=(1,L+1,L)=(L+1,1,L)=
c    (L+1,L+1,L)
       itype=nexttype(ConcInit)
       ijatom(1,1,L)=itype
       ijatom(1,1,0)=itype
       ijatom(1,L+1,0)=itype
       ijatom(L+1,1,0)=itype
       ijatom(L+1,L+1,0)=itype
       ijatom(1,L+1,L)=itype
       ijatom(L+1,1,L)=itype
       ijatom(L+1,L+1,L)=itype
c    (1,L,1)=(1,0,1)=(L+1,0,1)=(1,0,L+1)=
c    (L+1,0,L+1)=(L+1,L,1)=(1,L,L+1)=
c    (L+1,L,L+1)
       itype=nexttype(ConcInit)
       ijatom(1,L,1)=itype
       ijatom(1,0,1)=itype
       ijatom(L+1,0,1)=itype
       ijatom(1,0,L+1)=itype
       ijatom(L+1,0,L+1)=itype
       ijatom(L+1,L,1)=itype
       ijatom(1,L,L+1)=itype
       ijatom(L+1,L,L+1)=itype
c    (L,1,1)=(0,1,1)=(L,L+1,1)=(L,1,L+1)=
c    (L,L+1,L+1)=(0,L+1,1)=(0,1,L+1)=
c    (0,L+1,L+1)
       itype=nexttype(ConcInit)
       ijatom(L,1,1)=itype
       ijatom(0,1,1)=itype
       ijatom(L,L+1,1)=itype
       ijatom(L,1,L+1)=itype
       ijatom(L,L+1,L+1)=itype
       ijatom(0,L+1,1)=itype
       ijatom(0,1,L+1)=itype
       ijatom(0,L+1,L+1)=itype
c    (1,L,L)=(1,0,0)=(L+1,L,L)=(L+1,0,0)=
c    (1,L,0)=(1,0,L)=(L+1,L,0)=(L+1,0,L)
       itype=nexttype(ConcInit)
       ijatom(1,L,L)=itype
       ijatom(1,0,0)=itype
       ijatom(L+1,L,L)=itype
       ijatom(L+1,0,0)=itype
       ijatom(1,L,0)=itype
       ijatom(1,0,L)=itype
       ijatom(L+1,L,0)=itype
       ijatom(L+1,0,L)=itype
c    (L,1,L)=(0,1,0)=(L,L+1,L)=(0,L+1,0)=
c    (L,1,0)=(0,1,L)=(L,L+1,0)=(0,L+1,L)
       itype=nexttype(ConcInit)
       ijatom(L,1,L)=itype
       ijatom(0,1,0)=itype
       ijatom(L,L+1,L)=itype
       ijatom(0,L+1,0)=itype
       ijatom(L,1,0)=itype
       ijatom(0,1,L)=itype
       ijatom(L,L+1,0)=itype
       ijatom(0,L+1,L)=itype
c    (1,L,L)=(1,0,0)=(L+1,L,L)=(L+1,0,0)=
c    (1,L,0)=(1,0,L)=(L+1,L,0)=(L+1,0,L)
       itype=nexttype(ConcInit)
       ijatom(1,L,L)=itype
       ijatom(1,0,0)=itype
       ijatom(L+1,L,L)=itype
       ijatom(L+1,0,0)=itype
       ijatom(1,L,0)=itype
       ijatom(1,0,L)=itype
       ijatom(L+1,L,0)=itype
       ijatom(L+1,0,L)=itype
c    (L,L,L)=(0,0,0)=(L,L,0)=(L,0,L)=
c    (0,L,L)=(0,0,L)=(L,0,0)=(0,L,0)
       itype=nexttype(ConcInit)
       ijatom(L,L,L)=itype
       ijatom(0,0,0)=itype
       ijatom(L,L,0)=itype
       ijatom(L,0,L)=itype
       ijatom(0,L,L)=itype
       ijatom(0,0,L)=itype
       ijatom(L,0,0)=itype
       ijatom(0,L,0)=itype
c    (1,1,k)=(L+1,1,k)=(1,L+1,k)=(L+1,L+1,k)
       do 30 k=2,L-1
         itype=nexttype(ConcInit)
         ijatom(1,1,k)=itype
         ijatom(L+1,1,k)=itype
         ijatom(1,L+1,k)=itype
         ijatom(L+1,L+1,k)=itype
30     continue
c    (1,L,k)=(L+1,L,k)=(1,0,k)=(L+1,0,k)
       do 31 k=2,L-1
         itype=nexttype(ConcInit)
         ijatom(1,L,k)=itype
         ijatom(L+1,L,k)=itype
         ijatom(1,0,k)=itype
         ijatom(L+1,0,k)=itype
31     continue
c    (L,1,k)=(L,L+1,k)=(0,1,k)=(0,L+1,k)
       do 32 k=2,L-1
         itype=nexttype(ConcInit)
         ijatom(L,1,k)=itype
         ijatom(L,L+1,k)=itype
         ijatom(0,1,k)=itype
         ijatom(0,L+1,k)=itype
32     continue
c    (L,L,k)=(L,0,k)=(0,L,k)=(0,0,k)
       do 33 k=2,L-1
         itype=nexttype(ConcInit)
         ijatom(L,L,k)=itype
         ijatom(L,0,k)=itype
         ijatom(0,L,k)=itype
         ijatom(0,0,k)=itype
33     continue
c    (1,j,1)=(1,j,L+1)=(L+1,j,1)=(L+1,j,L+1)
       do 40 j=2,L-1
         itype=nexttype(ConcInit)
         ijatom(1,j,1)=itype
         ijatom(1,j,L+1)=itype
         ijatom(L+1,j,1)=itype
         ijatom(L+1,j,L+1)=itype
40     continue
c    (1,j,L)=(1,j,0)=(L+1,j,L)=(L+1,j,0)
       do 41 j=2,L-1
         itype=nexttype(ConcInit)
         ijatom(1,j,L)=itype
         ijatom(1,j,0)=itype
         ijatom(L+1,j,L)=itype
         ijatom(L+1,j,0)=itype
41     continue
c    (L,j,1)=(L,j,L+1)=(0,j,1)=(0,j,L+1)
       do 42 j=2,L-1
         itype=nexttype(ConcInit)
         ijatom(L,j,1)=itype
         ijatom(L,j,L+1)=itype
         ijatom(0,j,1)=itype
         ijatom(0,j,L+1)=itype
42     continue
c    (L,j,L)=(0,j,L)=(L,j,0)=(0,j,0)
       do 43 j=2,L-1
         itype=nexttype(ConcInit)
         ijatom(L,j,L)=itype
         ijatom(0,j,L)=itype
         ijatom(L,j,0)=itype
         ijatom(0,j,0)=itype
43     continue
c    (i,1,1)=(i,1,L+1)=(i,L+1,1)=(i,L+1,L+1)
       do 50 i=2,L-1
         itype=nexttype(ConcInit)
         ijatom(i,1,1)=itype
         ijatom(i,1,L+1)=itype
         ijatom(i,L+1,1)=itype
         ijatom(i,L+1,L+1)=itype
50     continue
c    (i,1,L)=(i,L+1,L)=(i,1,0)=(i,L+1,0)
       do 51 i=2,L-1
         itype=nexttype(ConcInit)
         ijatom(i,1,L)=itype
         ijatom(i,L+1,L)=itype
         ijatom(i,1,0)=itype
         ijatom(i,L+1,0)=itype
51     continue
c    (i,L,1)=(i,L,L+1)=(i,0,1)=(i,0,L+1)
       do 52 i=2,L-1
         itype=nexttype(ConcInit)
         ijatom(i,L,1)=itype
         ijatom(i,L,L+1)=itype
         ijatom(i,0,1)=itype
         ijatom(i,0,L+1)=itype
52     continue
c    (i,L,L)=(i,L,0)=(i,0,L)=(i,0,0)
       do 53 i=2,L-1
         itype=nexttype(ConcInit)
         ijatom(i,L,L)=itype
         ijatom(i,L,0)=itype
         ijatom(i,0,L)=itype
         ijatom(i,0,0)=itype
53     continue
c    (i,j,1)=(i,j,L+1)
       do 300 i=2,L-1
        do 301 j=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(i,j,1)=itype
	 ijatom(i,j,L+1)=itype
301     continue
300     continue
c    (i,j,L)=(i,j,0)
       do 302 i=2,L-1
        do 303 j=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(i,j,L)=itype
	 ijatom(i,j,0)=itype
303     continue
302     continue
c    (i,1,k)=(i,L+1,k)
       do 400 i=2,L-1
        do 401 k=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(i,1,k)=itype
	 ijatom(i,L+1,k)=itype
401     continue
400     continue
c    (i,L,k)=(i,0,k)
       do 402 i=2,L-1
        do 403 k=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(i,L,k)=itype
	 ijatom(i,0,k)=itype
403     continue
402     continue
c    (1,j,k)=(L+1,j,k)
       do 500 j=2,L-1
        do 501 k=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(1,j,k)=itype
	 ijatom(L+1,j,k)=itype
501     continue
500     continue
c    (L,j,k)=(0,j,k)
       do 502 j=2,L-1
        do 503 k=2,L-1
	 itype=nexttype(ConcInit)
	 ijatom(L,j,k)=itype
	 ijatom(0,j,k)=itype
503     continue
502     continue
c    (i,j,k) has no borders
       do 600 i=2,L-1
        do 601 j=2,L-1
	 do 602 k=2,L-1
	 ijatom(i,j,k)=nexttype(ConcInit)
602     continue
601     continue
600     continue
       return
      end
