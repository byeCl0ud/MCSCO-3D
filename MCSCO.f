        program MCSCO  
        
c Monte Carlo simulation of Ising-like SCO model in 3D

        implicit none 

        character inputFile*128
        character outputFile*128

        parameter( 
     1   inputFile='input.txt' ,
     2   outputFile='output.txt' )

       integer L, nTs, q, nMC, nrp
       real*8 D, Jcoop, r, Ti, Tf, ConcInit
       parameter( ConcInit=0.5 )

       WRITE(*,*) 'Welcome to SCO modeling with 3D MC'

       call readinput(
     1         inputFile,
     2         nMC,
     3         L, Ti, Tf, nTs, nrp,
     4         D, Jcoop, r, q )

      WRITE(*,*) 'Input load O.K.'

      call corepart(
     1          outputFile,
     2          nMC,
     3          L, Ti, Tf, nTs, nrp,
     4          D, Jcoop, r, q,
     5          ConcInit )

      WRITE(*,*) 'All done, no errors have  occured'

      end
