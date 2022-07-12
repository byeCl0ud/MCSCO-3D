        subroutine lattice(
     1       lt,
     2       nnnumber, nnnum_p,
     3       idx, jdy )
         implicit none
         integer lt
         integer nnnumber, nnnum_p
         integer idx(nnnum_p), jdy(nnnum_p)

         if( lt .eq. 1 )then
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
         else if( lt .eq. 2 )then
          write(6,*) ' Triangular lattice '
          nnnumber=6
          idx(1)=1
          jdy(1)=0
          idx(2)=0
          jdy(2)=1
          idx(3)=-1
          jdy(3)=1
          idx(4)=-1
          jdy(4)=0
          idx(5)=0
          jdy(5)=-1
          idx(6)=1
          jdy(6)=-1
         else
          write(6,*) ' Rombohedral lattice '
          nnnumber=2
          idx(1)=1
          jdy(1)=0
          idx(2)=-1
          jdy(2)=0
         end if

       return
       end

       subroutine random_initial(seed)
          integer seed
          real*8 r_int
          common /rand/r_int

          r_int = seed

       end
