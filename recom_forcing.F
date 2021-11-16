C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_forcing.F,v 1.49 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

CBOP
C     !ROUTINE: RECOM_FORCING
C     !INTERFACE:
      SUBROUTINE RECOM_FORCING( 
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE RECOM_FORCING
C     | o Calling interface for REcoM
C     *==================================================================
C     |
C     | RECOM_FORCING is the interface for calling a the bio-geo-chemical
C     | REcoM model (Markus Schartau)
C     | that provides source and sink terms for the passive tracers.
C     | the ptracers-package must be enabled for this interface to work
C     |
C     | RECOM_FORCING is called from GCHEM_FORCING, analogously to
C     | Stephanie Dutkiewicz's dic_forcing.
C     |
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
# ifdef ALLOW_GCHEM
#  include "GCHEM.h"
# endif
# ifdef ALLOW_RECOM
#  include "RECOM.h"
# endif /* ALLOW_RECOM */
#endif /* ALLOW_PTRACERS */
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
#endif /* ALLOW_SEAICE */
#ifdef ALLOW_THSICE
# include "THSICE_VARS.h"
#endif /* ALLOW_THSICE */
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
#endif /* ALLOW_EXF */
#ifdef ALLOW_LONGSTEP
# include "LONGSTEP.h"
#endif
C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     iMin, iMax, jMin, jMax - loop ranges
C     bi, bj - tile indices
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      INTEGER iMin, iMax, jMin, jMax
      integer bi, bj
      integer mythid 
      integer myiter
      _RL     myTime

#ifdef ALLOW_RECOM
C     !LOCAL VARIABLES :
C     == Local variables ==
C     i,j,k,iTracer - loop indices
C     kSurface      - surface index
C     krLoc         - local vertical loop boundary klowC
C     recom_dic     - dissolved inorganic carbon (intermediate variable)
C     recom_alk     - total alkalinity (intermediate variable)

C     benthicLayerOld/New - auxillary varialbes
C     angleOfInclination - of incoming solar irradiation
C     daylen        - length of the day (function of latitude and time of year)
C     radForce      - short wave radiation (needs to be replaced by a proper
C                     field !
      INTEGER i,j,k,iTracer,krLoc, kSurface
CCV      INTEGER iDIC, iALK, iFe

      _RL     recom_dic   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     recom_alk   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     dicFluxLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     hPlusLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     windSpeedLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     atmospCO2Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     saltLoc     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     thetaLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     co2_ppm     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     co2_mol     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef RECOM_CALCIFICATION
      _RL     benthicLayer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:4)
#else
      _RL     benthicLayer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:3)
#endif
C     these need to get the proper dimension ... etc.
      _RS     cosAngleOfIncidence(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     lengthOfDay(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     PARadiation(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     feDustFluxLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy) 
      _RS     surfpco2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     resuspension variables
      _RL     bottomStress
      _RL     erosionCoeff
      _RL     erosionFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CMS   implementing LithPartFlux
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
      _RS     ironLithParticleFluxLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif  /* RECOM_IRON_LITHOGENIC_PARTICLES */
C     diagnostics auxilliary fields
      INTEGER idiag, ndiags3d_used
      _RL     recom_diag_integral(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          nDiagsREcoM3D)
      _RL     recom_diag_3Dfields(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          Nr,nDiagsREcoM3D)
      _RL     recom_diag_2Dfields(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          nDiagsREcoM2D)
C    Local gptr :: work-space for time-stepping
      _RL     gPtrLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &               PTRACERS_num)
CEOP

CCV: tracer indices are now set in RECOM.h
CCV      iDIC = 2
CCV      iALK = 3
CCV      iFe  = PTRACERS_numInUse
C     surface layer index
      ksurface = 1
C     resuspension
      erosionCoeff = recom_FeErosionRate*(1.-recom_porosity)
     &     * recom_sizeFraction / recomCritShearStress

C     initialise some local fields
      DO j=jMin,jMax
       DO i=iMin,iMax
        feDustFluxLoc(i,j) = 0. _d 0
        erosionFlux  (i,j) = 0. _d 0
CMS   initialise LithParticleFlux too
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
        ironLithParticleFluxLoc(i,j) = 0. _d 0
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
       ENDDO
      ENDDO
      DO idiag = 1,nDiagsREcoM3D
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          recom_diag_3Dfields(i,j,k,idiag) = 0. _d 0
         ENDDO 
        ENDDO
       ENDDO
       DO j=jMin,jMax
        DO i=iMin,iMax
         recom_diag_integral(i,j,idiag) = 0. _d 0
        ENDDO 
       ENDDO
      ENDDO

C     solar radiation
      CALL RECOM_INSOLATION(
     O     cosAngleOfIncidence, lengthOfDay, PARadiation,
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )

#ifdef ALLOW_EXF
      IF ( recom_windFile .EQ. ' '  ) THEN
C     if no wind has been read copy wind speed of exf package (which
C     may or may not have been read from a file)
       DO j=jMin,jMax
        DO i=iMin,iMax
         windSpeed(i,j,bi,bj) = wspeed(i,j,bi,bj)
        ENDDO
       ENDDO  
      ENDIF 
#endif /* ALLOW_EXF */
C     Compute DIC Flux
      DO j=1-OLy, sNy+OLy
       DO i=1-OLx, sNx+OLx
C     initialise local arrays
C     recom_co2flux should not have to deal with dic/alk=0, so make
C     sure that recom_tiny is larger than zero.
        recom_dic(I,J) = MAX(recom_tiny,
     &       pTracer(I,J,kSurface,bi,bj,iDIC))
        recom_alk(I,J) = MAX(recom_tiny,
     &       pTracer(I,J,kSurface,bi,bj,iALK))
        dicFluxLoc(I,J)   = 0. _d 0
        co2_ppm(I,J)      = 0. _d 0
        co2_mol(I,J)      = 0. _d 0
        hPlusLoc(I,J)     = hPlus(I,J,bi,bj)
        windSpeedLoc(I,J) = windSpeed(I,J,bi,bj)
        atmospCO2Loc(I,J) = atmospCO2(I,J,bi,bj)
C        print*,'atmCO2', atmospCO2Loc
#ifdef ALLOW_LONGSTEP
        thetaLoc(I,J)     = LS_theta(I,J,kSurface,bi,bj)
C     this is a terrible hack and we should never have negative
C     salinities
        saltLoc(I,J)      = MAX(0.,LS_salt(I,J,kSurface,bi,bj) )
#else
        thetaLoc(I,J)     = theta(I,J,kSurface,bi,bj)
C     this is a terrible hack and we should never have negative
C     salinities
        saltLoc(I,J)      = MAX(0.,salt(I,J,kSurface,bi,bj) )
#endif /*  ALLOW_LONGSTEP */
CML        saltLoc(I,J)      = salt (I,J,kSurface,bi,bj)
CML        print '(A,2I4,3E10.2)','ml-dic',i,j,hPlus(I,J,bi,bj),
CML     &       pTracer(I,J,kSurface,bi,bj,iDIC),
CML     &       pTracer(I,J,kSurface,bi,bj,iALK)
       ENDDO
      ENDDO
      CALL RECOM_CO2FLUX(
     I     atmospCO2loc, recom_dic, recom_alk,
     I     thetaLoc, saltLoc, windSpeedLoc,
     O     dicFluxLoc, co2_ppm, co2_mol,
     U     hPlusLoc,
     I     iMin, iMax, jMin, jMax )
CML      DO j=jMin,jMax
CML       DO i=iMin,iMax
CML        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
CML         CALL RECOM_CO2FLUX(
CML     I        atmospCO2loc(I,J), recom_dic(I,J), recom_alk(I,J),
CML     I        thetaLoc(I,J), saltLoc(I,J), windSpeedLoc(I,J),
CML     O        dicFluxLoc(I,J), co2_ppm(I,J), co2_mol(I,J),
CML     U        hPlusLoc(I,J),
CML     I        iMin, iMax, jMin, jMax )
CMLCMLC co2_ppm_dummy is micro atmosheres indeed (compare recom_co2_flux)
CMLCML         print '(A,8E10.2)','ml-dic',
CMLCML     &        dicFluxLoc,atmospCO2(I,J,bi,bj),
CMLCML     &        recom_dic,recom_alk, 
CMLCML     &        theta(I,J,kSurface,bi,bj),salt(I,J,kSurface,bi,bj),
CMLCML     &        windSpeed(I,J,bi,bj), hPlus(I,J,bi,bj)    
CML        ENDIF
CML       ENDDO
CML      ENDDO
C     now copy back some variables
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
         surfpco2(I,J,bi,bj) = co2_ppm(I,J)   *maskC(I,J,kSurface,bi,bj)
         dicFlux (I,J,bi,bj) = dicFluxLoc(I,J)*maskC(I,J,kSurface,bi,bj)
         hPlus   (I,J,bi,bj) = hPlusLoc(I,J)  *maskC(I,J,kSurface,bi,bj)
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustDeposit(i,j,bi,bj)
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
CMS calculate ironLithogenicParticleFluxLoc
CMS feDustDeposit :: iron dust deposition at the surface, in mumol/m^2/s
CMS vgl.: RECOM.h l.150
CMS so get it in mug/m^2/s by multipliing with M = 56 g/mol
CMS fe is 3.5% of dust so reverse it by multipliing with 100/3.5
CMS divide by 1000 to get it in mg/m^2/s
	 ironLithParticleFluxLoc(i,j) = feDustDeposit(i,j,bi,bj)*56*0.1/3.5
#endif  /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDIF
       ENDDO
      ENDDO
C     Apply sea-ice mask according to the two different sea ice models
C     if necessary
#ifdef ALLOW_SEAICE
      IF ( useSeaice .AND. .NOT.useThSice ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         dicFlux(i,j,bi,bj) = dicFlux(i,j,bi,bj)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustFluxLoc(i,j)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
	 ironLithParticleFluxLoc(i,j) = ironLithParticleFluxLoc(i,j)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SEAICE */
#ifdef ALLOW_THSICE
      IF ( useThSice ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         dicFlux(i,j,bi,bj) = dicFlux(i,j,bi,bj)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustFluxLoc(i,j)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
	 ironLithParticleFluxLoc(i,j) = ironLithParticleFluxLoc(i,j)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDDO
       ENDDO
      ENDIF	
#endif /* ALLOW_THSICE */
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
C     determine loop ranges
         krLoc = klowC(i,j,bi,bj)
C     only if wet point:
C     surface forcing for individual tracers:
C     DIC
C     use surfaceForcingPtr to store forcing for the next time step to be
C     used with kpp (and pray that it is correct); IMPORTANT: convert to
C     units per second (dicFlux is in DIC/m^2/day)
         surfaceForcingPtr(i,j,bi,bj,iDIC) = 
     &        - dicFlux(i,j,bi,bj)/recom_secondsPerDay
CCV#ifdef ALLOW_RECOM_FEINPUT
C     Iron:
C     feDustFlux is iron, needs to be converted into dissolved iron, also
C     obey MITgcm sign conventions (positive = upward flux)
         surfaceForcingPtr(i,j,bi,bj,iFe) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * constantIronSolubility*feDustFluxLoc(i,j)
         pTracer(i,j,ksurface,bi,bj,iFe) =
     &        pTracer(i,j,ksurface,bi,bj,iFe)
     &        - surfaceForcingPtr(i,j,bi,bj,iFe)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
C     resuspension of sediment leads to iron flux into the ocean at 
C     the bottom, which is assumed to be proportional to an *imaginary* 
C     sediment erosion flux according to Ariathurai and Arulanandan (1978).
C     Stress units are m^2/s^2 which is a bit unusual for MITgcm, but
C     we stick with them for now. Also bottomStress has been scaled by
C     the bottomDragLinear/Quadratic already
#ifdef ALLOW_LONGSTEP 
         bottomStress = 0.25 _d 0 *
     &        ((LS_uVel(I,J,krLoc,bi,bj)+LS_uVel(I+1,J,krLoc,bi,bj))**2 
     &        +(LS_vVel(I,J,krLoc,bi,bj)+LS_vVel(I,J+1,krLoc,bi,bj))**2)
#else
         bottomStress = 0.25 _d 0 *
     &        ( (uVel(I,J,krLoc,bi,bj)+uVel(I+1,J,krLoc,bi,bj))**2 
     &        + (vVel(I,J,krLoc,bi,bj)+vVel(I,J+1,krLoc,bi,bj))**2 )
#endif /* ALLOW_LONGSTEP */
C     limit the stress and thus the erosion flux (may be useful in
C     very shallow regions with high velocities)
         bottomStress = MIN( bottomStress, recomShearStressMax )
         erosionFlux(i,j) = erosionCoeff
     &        * MAX( 0. _d 0, bottomStress-recomCritShearStress )
         pTracer(i,j,krLoc,bi,bj,iFe) =
     &        pTracer(i,j,krLoc,bi,bj,iFe)
     &        + erosionFlux(i,j)*dTtracerLev(krLoc)
     &        * recip_drF(krLoc) * _recip_hFacC(i,j,krLoc,bi,bj)
CMS Implementing RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
         surfaceForcingPtr(i,j,bi,bj,iDust) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * ironLithParticleFluxLoc(i,j)
c         surfaceForcingPtr(i,j,bi,bj,iDust) = 
c     &        - maskC(i,j,kSurface,bi,bj) * 1.0
         pTracer(i,j,ksurface,bi,bj,iDust) =
     &        pTracer(i,j,ksurface,bi,bj,iDust)
     &        - surfaceForcingPtr(i,j,bi,bj,iDust)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
CCV#endif /* ALLOW_RECOM_FEINPUT */
C     copy global fields to local arrays
         BenthicLayer(i,j,1) = benthicLayerN(i,j,bi,bj)
         benthicLayer(i,j,2) = benthicLayerC(i,j,bi,bj)
CCV#ifdef ALLOW_RECOM_SILICATE
         benthicLayer(i,j,3) = benthicLayerSi(i,j,bi,bj)
CCV#else /* not ALLOW_RECOM_SILICATE */
CCV         benthicLayer(i,j,3) = 0.
CCV#endif /* ALLOW_RECOM_SILICATE */
#ifdef RECOM_CALCIFICATION
         benthicLayer(i,j,4) = benthicLayerCaCO3(i,j,bi,bj)
#endif /* RECOM_CALCIFICATION */

C     end if wet point
        ENDIF
       ENDDO
      ENDDO
C     call the bio-geo-chemical model
      CALL RECOM_SMS(iMin,iMax,jMin,jMax,bi,bj,
     I     Nr, klowC, PTRACERS_num, nDiagsREcoM3D, nDiagsREcoM2D,
     I     dTtracerLev(kSurface), 
     I     maskC, hFacC, recip_hFacC, drF, recip_drF, recip_drC,
#ifdef ALLOW_LONGSTEP
     I     LS_theta, 
#else
     I     theta, 
#endif /* ALLOW_LONGSTEP */
     I     dicFlux, PARadiation, cosAngleOfIncidence, 
     I     pTracer,
     U     benthicLayer, 
     O     gPtrLoc, recom_diag_3Dfields, recom_diag_2Dfields )
C     copy new local tendency term to global array and perform
C     explicit (Eulerian) time step
      DO iTracer=1,PTRACERS_numInUse
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          pTracer(i,j,k,bi,bj,iTracer)=max(recom_tiny,
     &         pTracer(i,j,k,bi,bj,iTracer)
     &         +gPtrLoc(i,j,k,bi,bj,iTracer))*maskC(i,j,k,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     store benthic layer concentrations to global array
      DO j=jMin,jMax
       DO i=iMin,iMax
        benthicLayerN(i,j,bi,bj) = benthicLayer(i,j,1)
        benthicLayerC(i,j,bi,bj) = benthicLayer(i,j,2)
CCV#ifdef ALLOW_RECOM_SILICATE
        benthicLayerSi(i,j,bi,bj) = benthicLayer(i,j,3)
CCV#endif /* ALLOW_RECOM_SILICATE */
#ifdef RECOM_CALCIFICATION
        benthicLayerCaCO3(i,j,bi,bj) = benthicLayer(i,j,4)
#endif /* RECOM_CALCIFICATION */
       ENDDO
      ENDDO
#ifdef ALLOW_DIAGNOSTICS
      IF (useDiagnostics) THEN
C     vertical integral of 3D diagnostics
       DO idiag = 1,nDiagsRecom3D
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
           recom_diag_integral(i,j,idiag) = 
     &          recom_diag_integral(i,j,idiag) 
     &          + recom_diag_3Dfields(i,j,k,idiag)
     &          * drF(k) * _hFacC(I,J,k,bi,bj)
          Enddo
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DIAGNOSTICS
      IF (useDiagnostics) THEN
#ifdef RECOM_2CLASSES
       IF ( nDiagsREcoM3D .GT. 0 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,1),
     &                       'net_pps ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,1),
     &                       'NETPPVIS', 0,  1, 2, bi, bj, myThid)
       ENDIF
       IF ( nDiagsREcoM3D .GT. 1 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,2),
     &                       'net_ppd ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,2),
     &                       'NETPPVID', 0,  1, 2, bi, bj, myThid)
       ENDIF
#else
       IF ( nDiagsREcoM3D .GT. 0 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,1),
     &                       'net_pp  ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,1),
     &                       'NETPPVI ', 0,  1, 2, bi, bj, myThid)
       ENDIF
#endif /* RECOM_2CLASSES */
#ifdef RECOM_2CLASSES
       IF ( nDiagsREcoM3D .GT. 2 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,3),
     &                       'gr_pps  ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,3),
     &                       'GRPPVIS ', 0,  1, 2, bi, bj, myThid)
       ENDIF
       IF ( nDiagsREcoM3D .GT. 3 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,4),
     &                       'gr_ppd  ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,4),
     &                       'GRPPVID ', 0,  1, 2, bi, bj, myThid)
       ENDIF
#else
       IF ( nDiagsREcoM3D .GT. 1 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,2),
     &                       'gross_pp', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,2),
     &                       'GROSPPVI', 0,  1, 2, bi, bj, myThid)
       ENDIF
#endif /* RECOM_2CLASSES */
#ifdef RECOM_2CLASSES
       IF ( nDiagsREcoM3D .GT. 4 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,5),
     &                       'net_nass', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,5),
     &                       'NETNAVIS', 0,  1, 2, bi, bj, myThid)
       ENDIF
       IF ( nDiagsREcoM3D .GT. 5 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,6),
     &                       'net_nasd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,6),
     &                       'NETNAVID', 0,  1, 2, bi, bj, myThid)
       ENDIF
#else
       IF ( nDiagsREcoM3D .GT. 2 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,3),
     &                       'net_nas ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,3),
     &                       'NETNAVI ', 0,  1, 2, bi, bj, myThid)
       ENDIF
#endif /* RECOM_2CLASSES */
#ifdef RECOM_2CLASSES
       IF ( nDiagsREcoM3D .GT. 6 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,7),
     &                       'n_assims', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,7),
     &                       'GRNAVIS ', 0,  1, 2, bi, bj, myThid)
       ENDIF
       IF ( nDiagsREcoM3D .GT. 7 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,8),
     &                       'n_assimd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,8),
     &                       'GRNAVID ', 0,  1, 2, bi, bj, myThid)
       ENDIF
#else
       IF ( nDiagsREcoM3D .GT. 3 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,4),
     &                       'n_assim', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,4),
     &                       'GROSNAVI', 0,  1, 2, bi, bj, myThid)
       ENDIF
#endif /* RECOM_2CLASSES */
#ifdef RECOM_MAREMIP
       ndiags3d_used = 8
       if ( nDiagsREcom3D . GE. (ndiags3d_used+28) ) then
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+1),
     &            'reminer ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+2),
     &            'graztot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+3),
     &            'resptot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+4),
     &            'calcprod', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+5),
     &            'calcdiss', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+6),
     &            'rdoczoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+7),
     &            'rpoczoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+8),
     &            'cdiapoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+9),
     &            'cphypoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+10),
     &            'cdiadoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+11),
     &            'cphydoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+12),
     &            'cresdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+13),
     &            'cresphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+14),
     &            'creszoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+15),
     &            'grazdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+16),
     &            'grazphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+17),
     &            'par3d   ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+18),
     &            'nlimdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+19),
     &            'nlimphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+20),
     &            'llimdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+21),
     &            'llimphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+22),
     &            'felimdia', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+23),
     &            'felimphy', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+24),
     &            'silimdia', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+25),
     &            'sink_pon', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+26),
     &            'sink_poc', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+27),
     &            'sink_bsi', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+28),
     &            'sinkcalc', 0, Nr, 2, bi, bj, myThid)
       endif
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
       ndiags3d_used = 8
#ifdef RECOM_MAREMIP
     &    + 28
#endif
       if ( nDiagsREcom3D . GE. (ndiags3d_used+6) ) then
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+1),
     &            'ligsrrem', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+2),
     &            'ligsrdoc', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+3),
     &            'ligremin', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+4),
     &            'ligphot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+5),
     &            'liguptak', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+6),
     &            'lig_coag', 0, Nr, 2, bi, bj, myThid)
       endif
#endif
       IF ( nDiagsREcoM2D .GE. 4) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,1),
     &                       'SEDFN   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,2),
     &                       'SEDFC   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,3),
     &                       'SEDFSI  ', 0,  1, 2, bi, bj, myThid)
#ifdef RECOM_CALCIFICATION
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,4),
     &                       'SEDFCALC', 0,  1, 2, bi, bj, myThid)
#endif
       ENDIF
       IF ( nDiagsREcoM2D .GE. 9) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,5),
     &                       'DIAFN   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,6),
     &                       'DIAFC   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,7),
     &                       'DIAFALK ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,8),
     &                       'DIAFSI  ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,9),
     &                       'DIAFFE  ', 0,  1, 2, bi, bj, myThid)
       ENDIF

#ifdef RECOM_EXPORT_DIAGNOSTICS
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,10),
     &                       'EXPORTN ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,11),
     &                       'EXPORTC ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,12),
     &                       'EXPCALC ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,13),
     &                       'EXPORTSI', 0,  1, 2, bi, bj, myThid)
#endif

       CALL DIAGNOSTICS_FILL(PARadiation,
     &                      'PARSURF ', 0,  1, 2, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(surfaceForcingPtr(1-Olx,1-Oly,1,1,iDIC),
     &                      'CO2Flx  ', 0,  1, 1, bi, bj, myThid) 
       CALL DIAGNOSTICS_FILL(surfpco2,
     &                      'pCO2surf ', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(hPlus,
     &                      'HPlus   ', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(surfaceForcingPtr(1-Olx,1-Oly,1,1,iFe),
     &                      'FeSrfFlx', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(erosionFlux,
     &                      'FeBtmFlx', 0,  1, 2, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(benthicLayerN,
     &                      'NBENTHOS', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(benthicLayerC,
     &                      'CBENTHOS', 0,  1, 1, bi, bj, myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
       CALL DIAGNOSTICS_FILL(benthicLayerSi,
     &                      'SBENTHOS', 0,  1, 1, bi, bj, myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
#ifdef RECOM_CALCIFICATION
       CALL DIAGNOSTICS_FILL(benthicLayerCaCO3,
     &                      'CALCBENT', 0,  1, 1, bi, bj, myThid)
#endif /* RECOM_CALCIFICATION */ 
      ENDIF
      
        ndiags3d_used = 8
#ifdef RECOM_MAREMIP
     &    + 28
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
#ifdef RECOM_IRON_TWO_LIGANDS
     &    + 11
#else        
     &    + 6
#endif
#endif

      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+1),
     &                       'chlsyphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+2),
     &                       'chlsydia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+3),
     &                       'kochlphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+4),
     &                       'kochldia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+5),
     &                       'ppmaxphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+6),
     &                       'ppmaxdia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+7),
     &                       'pphotphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+8),
     &                       'pphotdia', 0, Nr, 2, bi, bj, myThid)      
      
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RECOM */

      RETURN
      END

#ifdef ALLOW_RECOM
      SUBROUTINE RECOM_FORCING_DUMMY( 
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "RECOM.h"
      INTEGER iMin, iMax, jMin, jMax
      integer bi, bj
      integer mythid 
      integer myiter
      _RL     myTime
      INTEGER i,j,k,iTracer,krLoc, kSurface
      ksurface = 1
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
C--------------------------------------
CCV The meaning of the following part is completely unclear to me. Delete?
CCV#ifdef ALLOW_RECOM_FEINPUT
C     Iron:
         iTracer = 4
C     feDustFlux is iron, needs to be converted into dissolved iron, also
C     obey MITgcm sign conventions (positive = upward flux)
         surfaceForcingPtr(i,j,bi,bj,iTracer) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * constantIronSolubility*feDustDeposit(i,j,bi,bj)
         pTracer(i,j,ksurface,bi,bj,iTracer) =
     &        pTracer(i,j,ksurface,bi,bj,iTracer)
     &        - surfaceForcingPtr(i,j,bi,bj,iTracer)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
         pTracer(i,j,ksurface,bi,bj,5) =
     &        pTracer(i,j,ksurface,bi,bj,5)
     &        - surfaceForcingPtr(i,j,bi,bj,iTracer)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
CCV#endif /* ALLOW_RECOM_FEINPUT */
C--------------------------------------
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END
#endif /* ALLOW_RECOM */
#ifdef ALLOW_RECOM
CBOP
C     !ROUTINE: BGC_DUMMY_MODEL
C     !INTERFACE:
      SUBROUTINE BGC_DUMMY_MODEL( 
     INr, kr, nTracer, rC, drF, bgc_temp, bgc_ptracer, 
     Obgc_gPtr, 
     ImyTime, myIter, myThid 
     &)

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE BGC_DUMMY_MODEL
C     | o Dummy routine for a bio-geo-chemical model
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     Nr - number of vertical points
C     kr - index of lowest wet tracer grid cell (k>kr means below bottom)
C     nTracer - number biological tracers
C     rC - depth of grid cell centers
C     drF - thickness of vertical grid cells
C     bgc_temp - local vertical array of temperature (potential!)
C     bgc_pTracer - local vertical array of biological tracers
C     bgc_gPtr - local vertical array of biological tracer tendencies
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      INTEGER Nr, kr, nTracer
      _RS     rC(1:Nr),drF(1:Nr)
      _RL     bgc_temp(1:Nr)
      _RL     bgc_pTracer(1:Nr,1:nTracer)
      _RL     bgc_gPtr(1:Nr,1:nTracer)
      integer mythid 
      integer myiter
      _RL     myTime

C     !LOCAL VARIABLES :
C     == Local variables ==
      INTEGER k, iTracer
      CHARACTER*(3) NVAR
      CHARACTER*(128) FMT
CEOP

C--   the bio-geo-chemical model goes here
C     this is just a place-holder
      print *, 'ml-bgc: myIter = ', myIter, ', myTime = ', myTime
      print *, 'ml-bgc:   k, bgc_pTracer'
      write(NVAR,'(I3)') nTracer
      FMT = '(A,I3,2F10.4,'//NVAR//'E10.2)'
      print *, 'ml-bgc:   k     depth  bgc_temp  bgc_pTracer'
      DO k = 1,Nr
       write(*,FMT) ' ml-bgc: ', k, rC(k), bgc_temp(k),
     &      (bgc_pTracer(k,iTracer),iTracer=1,nTracer)
      ENDDO 
      print *, 'ml-bgc:   k     depth  bgc_temp  bgc_gpTr'
      DO k = 1,Nr
       DO iTracer = 1,nTracer
        bgc_gPtr(k,iTracer) = float(k+iTracer)
       ENDDO
       write(*,FMT) ' ml-bgc: ', k, rC(k), bgc_temp(k),
     &     (bgc_gPtr(k,iTracer),iTracer=1,nTracer)
      ENDDO
C--

      RETURN
      END
#endif /* ALLOW_RECOM */

