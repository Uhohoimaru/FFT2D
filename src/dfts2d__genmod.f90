        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct  7 11:37:04 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DFTS2D__genmod
          INTERFACE 
            SUBROUTINE DFTS2D(NX,NY,FX,FY,Z,MFX,MFY,MZ)
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              REAL(KIND=8), INTENT(IN) :: FX(0:NX-1)
              REAL(KIND=8), INTENT(IN) :: FY(0:NY-1)
              COMPLEX(KIND=8) :: Z(0:NX-1,0:NY-1)
              REAL(KIND=8), INTENT(OUT) :: MFX(0:NX-1)
              REAL(KIND=8), INTENT(OUT) :: MFY(0:NY-1)
              COMPLEX(KIND=8), INTENT(OUT) :: MZ(0:NX-1,0:NY-1)
            END SUBROUTINE DFTS2D
          END INTERFACE 
        END MODULE DFTS2D__genmod
