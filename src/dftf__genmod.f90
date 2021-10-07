        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct  7 11:37:04 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DFTF__genmod
          INTERFACE 
            SUBROUTINE DFTF(N,F,H)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: F(0:N-1)
              REAL(KIND=8), INTENT(IN) :: H
            END SUBROUTINE DFTF
          END INTERFACE 
        END MODULE DFTF__genmod
