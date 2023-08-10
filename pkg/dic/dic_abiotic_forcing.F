#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_ABIOTIC_FORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_ABIOTIC_FORCING(
     I                      PTR_DIC, 
     I                      PTR_ALK,
     I                      PTR_PO4,
#ifdef ALLOW_O2
     I                      PTR_O2 ,
#endif
     I                      bi, bj, iMin, iMax, jMin, jMax,
     I                      myIter, myTime, myThid )

C !DESCRIPTION:
C updates all the tracers for the effects of air-sea exchange

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C  PTR_DIC              :: dissolved inorganic carbon
C  PTR_ALK              :: alkalinity
C  PTR_O2               :: oxygen
C  bi, bj               :: current tile indices
C  myIter               :: current timestep
C  myTime               :: current time
C  myThid               :: thread number
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gDIC   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gALK   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_O2
      _RL  PTR_O2 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  gO2    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
      INTEGER bi, bj, iMin, iMax, jMin, jMax
      INTEGER myIter
      _RL myTime
      INTEGER myThid
#ifdef ALLOW_PTRACERS

C !LOCAL VARIABLES: ====================================================
C  i,j,k                  :: loop indices
C  SURA                   :: tendency of alkalinity due to freshwater
C  SURC                   :: tendency of DIC due to air-sea exchange
C                            and virtual flux
C  SURO                   :: tendency of O2 due to air-sea exchange
      _RL  SURA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURO(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER i,j,k
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('DIC_ABIOTIC_FORCING',myThid)
#endif
      IF ( useThSIce .OR. useSEAICE .OR. useCoupler ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('DIC_FIELDS_UPDATE',myThid)
#endif
        CALL DIC_FIELDS_UPDATE(
     I                    bi, bj, myTime, myIter, myThid )
      ENDIF
      
       DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           SURA(i,j)   =0. _d 0
           SURC(i,j)   =0. _d 0
           SURO(i,j)   =0. _d 0
         ENDDO
       ENDDO

C carbon air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('DIC_SURFFORCING',myThid)
#endif
       CALL DIC_SURFFORCING(
     I                    PTR_DIC, PTR_ALK, PTR_PO4,
     O                    SURC,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

C alkalinity air-sea interaction
C Only if def ALLOW_OLD_VIRTUALFLUX is anything actually done
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('ALK_SURFFORCING',myThid)
#endif
       CALL ALK_SURFFORCING(
     I                    PTR_ALK,
     O                    SURA,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

#ifdef ALLOW_O2
C oxygen air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('O2_SURFFORCING',myThid)
#endif
       CALL O2_SURFFORCING(
     I                    PTR_O2,
     O                    SURO,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#endif

       DO j=jMin,jMax
         DO i=iMin,iMax
               gALK(i,j,1)=gALK(i,j,1)+SURA(i,j)
               gDIC(i,j,1)=gDIC(i,j,1)+SURC(i,j)
#ifdef ALLOW_O2
               gO2(i,j,1) =gO2(i,j,1)+SURO(i,j)
#endif
         ENDDO
       ENDDO
       IF ( useOBCS ) THEN
         DO j=jMin,jMax
          DO i=iMin,iMax
            GDIC(i,j,1) = GDIC(i,j,1)*maskInC(i,j,bi,bj)
            GALK(i,j,1) = GALK(i,j,1)*maskInC(i,j,bi,bj)
#ifdef ALLOW_O2
            GO2(i,j,1)  =  GO2(i,j,1)*maskInC(i,j,bi,bj)
#endif
          ENDDO
         ENDDO
       ENDIF

#ifdef ALLOW_TIMEAVE
C save averages
       IF ( PTRACERS_taveFreq.GT.0. ) THEN
         DO j=jMin,jMax
          DO i=iMin,iMax
              SURave(i,j,bi,bj)    =SURave(i,j,bi,bj)+
     &                              SURC(i,j)*deltaTClock
#ifdef ALLOW_O2
              SUROave(i,j,bi,bj)   =SUROave(i,j,bi,bj)+
     &                              SURO(i,j)*deltaTClock
#endif
              pCO2ave(i,j,bi,bj)   =pCO2ave(i,j,bi,bj)+
     &                              pCO2(i,j,bi,bj)*deltaTClock
              pHave(i,j,bi,bj)     =pHave(i,j,bi,bj)+
     &                              pH(i,j,bi,bj)*deltaTClock
              fluxCO2ave(i,j,bi,bj)=fluxCO2ave(i,j,bi,bj)+
     &                           fluxCO2(i,j,bi,bj)*deltaTClock
          ENDDO
         ENDDO
C         DIC_timeAve(bi,bj) = DIC_timeAve(bi,bj)+deltaTClock
       ENDIF
#endif /* ALLOW_TIMEAVE*/

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_FILL(pCO2   ,'DICPCO2 ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(fluxCO2,'DICCFLX ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH     ,'DICPHAV ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(SURC   ,'DICTFLX ',0,1 ,2,bi,bj,myThid)
#ifdef ALLOW_O2
        CALL DIAGNOSTICS_FILL(SURO   ,'DICOFLX ',0,1 ,2,bi,bj,myThid)
#endif
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('DIC_ABIOTIC_FORCING',myThid)
#endif
#endif /* ALLOW_PTRACERS */

       RETURN
       END
