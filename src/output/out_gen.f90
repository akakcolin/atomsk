module out_gen
  use comv
  use constants
  use functions
  use messages
  use files
  use sorting
  use subroutines

  implicit none
contains


  SUBROUTINE write_gen(H,P,comment,AUXNAMES,AUX,outputfile)
    !
    CHARACTER(LEN=*),INTENT(IN):: outputfile
    CHARACTER(LEN=1):: answer
    CHARACTER(LEN=2):: species

    CHARACTER(LEN=4096):: msg, temp
    CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
    CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
    LOGICAL:: isreduced
    INTEGER:: i, j, last, tmp_id

    INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
    LOGICAL:: contiguous   !are atom with same species contiguous? must they be packed?
    REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
    REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
    REAL(dp),DIMENSION(:,:),POINTER:: Ppoint  !pointer to P or Q
    REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),TARGET:: P
    REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: Q  !copy of P (if needed)
    REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),TARGET:: AUX !auxiliary properties
    REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: AUX2 !auxiliary properties (temporary)
    REAL(dp),DIMENSION(:,:),POINTER:: AUXpoint  !pointer to AUX or AUX2

    CHARACTER(LEN=2), dimension(:), allocatable  ::elements
    !
    !
    !Initialize variables

    tmp_id =0
    contiguous = .TRUE.

    !
    WRITE(msg,*) 'entering WRITE_GEN'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !
    !
100 CONTINUE
    !First make sure that atoms of same species are all contiguous
    DO i=1,SIZE(P,1)
       last=i
       DO j=i+1,SIZE(P,1)
          IF( P(j,4)==P(i,4) ) THEN
             IF(j==last+1) THEN
                last=j
             ELSE
                contiguous=.FALSE.
                GOTO 120
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
120 CONTINUE
    IF(.NOT.contiguous) THEN
       !Display a message
       nwarn=nwarn+1
       CALL ATOMSK_MSG(3703,(/msg,outputfile/),(/0.d0/))
       READ(*,*) answer
       IF(answer.NE.langyes .AND. answer.NE.langBigYes) THEN
          contiguous = .TRUE.
       ENDIF
    ENDIF
    !
    IF(.NOT.contiguous) THEN
       !If species are not contiguous then rearrange them
       !Note that the arrays P and AUX must remain untouched
       !Therefore the columns fixx, fixy, fixz of AUX are copied into Q for sorting


       ALLOCATE( Q( SIZE(P,1), SIZE(P,2) ) )
       Q(:,:) = P(:,:)

       !
       !Sort atoms (and aux. prop.) to pack atoms with the same species
       CALL PACKSORT(Q,4,newindex)
       Ppoint=>Q
       !
       IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)>0 ) THEN
          !Save which atoms are fixed in AUX2
          ALLOCATE( AUX2( SIZE(AUX,1) , 3 ) )
          AUX2(:,1) = Q(:,5)
          AUX2(:,2) = Q(:,6)
          AUX2(:,3) = Q(:,7)

          AUXpoint=>AUX2
       ENDIF
       !
       !Count the number of atoms for each atomic species
       CALL FIND_NSP(Ppoint(:,4),atypes)
       !
       CALL ATOMSK_MSG(3003,(/''/),(/atypes(:,1)/))
       !
    ELSE
       !If everything is fine, just use P in the following
       Ppoint=>P
       AUXpoint=>AUX
       !Count the number of atoms for each atomic species
       CALL FIND_NSP(Ppoint(:,4),atypes)
    ENDIF
    !
    !
    !
200 CONTINUE
    WRITE(msg,*) "Number of different species: ", SIZE(atypes(:,1))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    IF(ofu.NE.6) THEN
       OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
    ENDIF

    !
    ALLOCATE( elements( size(atypes,1) ))

    CALL FIND_IF_REDUCED(H,P,isreduced)
    !Write atom coordinates
    IF(isreduced) THEN
        WRITE(ofu,*) " ", size(P,1), ' F'
    ELSE
        WRITE(ofu,*) " ", size(P,1), ' S'
    ENDIF

    temp=""
    DO i=1,SIZE(atypes,1)
      CALL ATOMSPECIES(atypes(i,1),species)
      temp = TRIM(ADJUSTL(temp))//"  "//species
    ENDDO
    WRITE(ofu,*) TRIM(ADJUSTL(temp))

    !
    210 FORMAT(i6,2X, i6, 2X,  3(f16.8,2X))

    DO i=1,SIZE(Ppoint,1)
       call ATOMSPECIES(Ppoint(i,4), species)
       if(.not.ANY(elements==species)) then
               tmp_id = tmp_id+1
               elements(tmp_id)=species
       end if

       WRITE(ofu,210) i, tmp_id, Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)

    ENDDO


    write(ofu,*) '0.0000 0.000 0.0000'

    WRITE(ofu,201) H(1,1), H(1,2), H(1,3)
    WRITE(ofu,201) H(2,1), H(2,2), H(2,3)
    WRITE(ofu,201) H(3,1), H(3,2), H(3,3)
    201 FORMAT(3(f16.8,2X))
    GOTO 500
    !
    !
    !
500 CONTINUE
    CLOSE(ofu)
    msg = "GEN"
    temp = outputfile
    CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
    !
    !
    !
1000 CONTINUE
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
    IF(ALLOCATED(elements)) DEALLOCATE(elements)
    NULLIFY(Ppoint)
    NULLIFY(AUXpoint)
    !
    !
    !
  END SUBROUTINE write_gen


end module out_gen
