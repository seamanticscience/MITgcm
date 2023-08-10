#include "DIC_OPTIONS.h"

CStartOfInterface
      SUBROUTINE DIC_ABIOTIC_DIAGS( myTime, myIter, myThid )
C     *==========================================================*
C     | SUBROUTINE DIC_ABIOTIC_DIAGS
C     | o Do DIC abiotic diagnostic output
C     *==========================================================*
C     *==========================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"

C     == Routine arguments ==
C     myTime - Current time of simulation ( s )
C     myIter - Iteration number
C     myThid -  Number of this instance of INI_FORCING
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEndOfInterface

#ifdef ALLOW_TIMEAVE
C     == Local variables ==
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      INTEGER bi, bj
      CHARACTER*(10) suff
#ifdef ALLOW_MNC
      CHARACTER*(1) pf
#endif

C     Dump files and restart average computation if needed
      IF (DIFFERENT_MULTIPLE(PTRACERS_taveFreq,myTime,deltaTClock)) THEN

C      Normalize by integrated time
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         CALL TIMEAVE_NORMALIZE(SURave,DIC_AbioTimeAve, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(SUROave,DIC_AbioTimeAve, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(pCO2ave,DIC_AbioTimeAve, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(pHave,DIC_AbioTimeAve, 1 ,
     &                 bi,bj,myThid)
         CALL TIMEAVE_NORMALIZE(fluxCO2ave,DIC_AbioTimeAve, 1 ,
     &                 bi,bj,myThid)
        ENDDO
       ENDDO

       IF ( PTRACERS_timeave_mdsio )  THEN
         IF ( rwSuffixType.EQ.0 ) THEN
           WRITE(suff,'(I10.10)') myIter
         ELSE
           CALL RW_GET_SUFFIX( suff, myTime, myIter, myThid )
         ENDIF
         CALL WRITE_FLD_XY_RL('DIC_Surtave.',suff,SURave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_SurOtave.',suff,SUROave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_pHtave.',suff,pHave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_pCO2tave.',suff,pCO2ave,
     &        myIter,myThid)
         CALL WRITE_FLD_XY_RL('DIC_fluxCO2ave.',suff,fluxCO2ave,
     &        myIter,myThid)
       ENDIF

#ifdef ALLOW_MNC
       IF ( PTRACERS_timeave_mnc ) THEN
         IF ( writeBinaryPrec .EQ. precFloat64 ) THEN
          pf(1:1) = 'D'
         ELSE
          pf(1:1) = 'R'
         ENDIF
         CALL MNC_CW_SET_UDIM('dic_tave', -1, myThid)
         CALL MNC_CW_RL_W_S('D','dic_tave',0,0,'T',myTime,myThid)
         CALL MNC_CW_SET_UDIM('dic_tave', 0, myThid)
         CALL MNC_CW_I_W_S('I','dic_tave',0,0,'iter',myIter,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_SUR_ave',SURave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_SURO_ave',SUROave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_pH_ave',pHave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_pCO2_ave',pCO2ave,myThid)
         CALL MNC_CW_RL_W(
     &        pf,'dic_tave',0,0,'dic_fluxCO2_ave',fluxCO2ave,myThid)
       ENDIF
#endif

C      Reset averages to zero
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         CALL TIMEAVE_RESET(SURave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(SUROave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(pCO2ave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(pHave,1,bi,bj,myThid)
         CALL TIMEAVE_RESET(fluxCO2ave,1,bi,bj,myThid)
         DIC_AbioTimeAve(bi,bj) = 0. _d 0
        ENDDO
       ENDDO

      ENDIF
#endif /* ALLOW_TIMEAVE */
      RETURN
      END
