        subroutine readinput(
     1         inputFile,
     2         nMC,
     3         L, Ti, Tf, nTs, nrp,
     4         D, Jcoop, r, q )

        character inputFile*128
        integer L, nTs, q, nrp
        real*8 D, Jcoop, r, Ti, Tf

       open(UNIT=11, FILE=inputFile,
     ^      ACCESS='SEQUENTIAL',
     ^      STATUS='OLD' )

          read(11,*) L
          read(11,*) Ti
          read(11,*) Tf
          read(11,*) nTS
          read(11,*) nrp
          read(11,*) D
          read(11,*) Jcoop
          read(11,*) r
          read(11,*) q
          read(11,*) seed
       close(11)

c        WRITE(*,*) 'D=',D,' J=',Jcoop,' r=',r

        r=log(r)
        call random_initial(seed)
        nMC=L+2
        return
        end
