#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_ABIOTIC_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_ABIOTIC_INIT( myThid )

C !DESCRIPTION:
C  DIC diagnostic initialization

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"

c!INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_TIMEAVE

C !LOCAL VARIABLES: ====================================================
      INTEGER bi, bj
c     CHARACTER*(MAX_LEN_MBUF) suff

C set arrays to zero if first timestep
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          CALL TIMEAVE_RESET(SURave,   1,  bi, bj, myThid)
          CALL TIMEAVE_RESET(SUROave,   1,  bi, bj, myThid)
          CALL TIMEAVE_RESET(pCO2ave,   1,  bi, bj, myThid)
          CALL TIMEAVE_RESET(pHave,   1,  bi, bj, myThid)
          CALL TIMEAVE_RESET(fluxCO2ave,   1,  bi, bj, myThid)
c         print*,'QQ start timeave'
          DIC_AbioTimeAve(bi,bj) = 0. _d 0
         ENDDO
        ENDDO
#endif /* ALLOW_TIMEAVE */

      RETURN
      END
