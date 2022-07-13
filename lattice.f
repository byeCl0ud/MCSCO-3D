        subroutine lattice(
     2       nnnumber, nnnum_p,
     3       idx, jdy )
         implicit none
         integer nnnumber, nnnum_p
         integer idx(nnnum_p), jdy(nnnum_p)

          write(6,*) ' Simple square lattice '
          nnnumber=4
          idx(1)=1
          jdy(1)=0
          idx(2)=0
          jdy(2)=1
          idx(3)=-1
          jdy(3)=0
          idx(4)=0
          jdy(4)=-1

       return
       end

       subroutine random_initial(seed)
          integer seed
          real*8 r_int
          common /rand/r_int

          r_int = seed

       end
