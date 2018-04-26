    call init(randGenPool, env, iseed, oldCompat=.true.)
    call randGenPool%getGenerator(env, randomFreq)
    allocate(eigenValues(3*nAtom))
    allocate(basis(3*nAtom, 3*nAtom))
    allocate(coordcm(3,nAtom))
    allocate(invpartner(nAtom))
    allocate(newhess(3*nAtom, 3*nAtom))
    allocate(dotpr(3*nAtom, 3*nAtom))
    allocate(redmas(3*nAtom))

    n3 = 3*nAtom

    !mass scaling
    iCount = 0
    do ii = 1, nAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nAtom
          do ll = 1, 3
            jCount = jCount + 1
            pDynMatrix(jCount,iCount) = pDynMatrix(jCount,iCount) &
                & / (sqrt(mass(ii)) * sqrt(mass(jj)))
          end do
        end do
      end do
    end do


    xvec = 0.0_dp
    do ii = 1, nAtom
        do jj = 1, 3
            xvec(jj) = xvec(jj) + mass(ii)*coord0(jj,ii)
        end do
    end do
    prod = sum(mass(:))
    xvec(:) = xvec(:) / prod

    do ii = 1, nAtom
        coordcm(:,ii) = coord0(:,ii) - xvec(:)
    end do

    centinv = .false.
    diffinv = 0.0_dp
    do ii = 1, nAtom
        distinv = 1.0d2
        do jj = 1, nAtom
            if( ii == jj ) then
                do kk = 1, 3
                    xvec(kk) = coordcm(kk,ii)
                end do
            else 
                do kk = 1, 3
                    xvec(kk) = coordcm(kk,ii) + coordcm(kk,jj)
                end do
            end if   
            rr = sqrt(dot_product(xvec,xvec))
            if(rr<distinv) then
                distinv = rr
                invpartner(ii) = jj
                if( ii == jj ) then
                    invpartner(ii) = 0
                end if
            end if
        end do
        diffinv = diffinv+distinv
    end do
    if(diffinv < 1.0E-2_dp) then
        centinv = .true.
    end if


    !############################
    !## construct translations ##
    !############################
    basis = 0.0_dp
    do ii=1, nAtom
      basis(1+3*(ii-1),1) = sqrt(mass(ii))
      basis(2+3*(ii-1),2) = sqrt(mass(ii))
      basis(3+3*(ii-1),3) = sqrt(mass(ii))
    end do
    !#########################
    !## construct rotations ##
    !#########################
    do ii=1, nAtom
      basis(2+3*(ii-1),4)= -sqrt(mass(ii))*coordcm(3,ii)
      basis(3+3*(ii-1),4)=  sqrt(mass(ii))*coordcm(2,ii)
      basis(3+3*(ii-1),5)= -sqrt(mass(ii))*coordcm(1,ii)
      basis(1+3*(ii-1),5)=  sqrt(mass(ii))*coordcm(3,ii)
      basis(1+3*(ii-1),6)= -sqrt(mass(ii))*coordcm(2,ii)
      basis(2+3*(ii-1),6)=  sqrt(mass(ii))*coordcm(1,ii)
    end do

    !######################################################
    !## overlaping matrix of the rotational eigenvectors ##
    !######################################################
    dotpr = 0.0_dp
    do ii = 4, 6
      do jj = 4, 6
        dotpr(ii,jj) = dot_product(basis(:,ii),basis(:,jj))
      end do
    end do

  !#################################################################
  !## determine the number of independent rotational eigenvectors ##
  !#################################################################
    det = dotpr(4,4)*dotpr(5,5)*dotpr(6,6)+dotpr(4,5)*dotpr(5,6)*dotpr(6,4)
    det = det+dotpr(4,6)*dotpr(5,4)*dotpr(6,5)-dotpr(4,5)*dotpr(5,4)*dotpr(6,6)
    det = det-dotpr(4,6)*dotpr(5,5)*dotpr(6,4)-dotpr(4,4)*dotpr(5,6)*dotpr(6,5)

    if(det<offsetlinear) then
      nindep = 5
    else 
      nindep = 6
    end if


    if(nindep==5) then
      det1 = dotpr(4,4)*dotpr(5,5)-dotpr(4,5)*dotpr(5,4)
      det2 = dotpr(4,4)*dotpr(6,6)-dotpr(4,6)*dotpr(6,4)
      if( det1 < offsetlinear ) then
        if( det2 < offsetlinear ) then
          basis(:,4) = basis(:,6)
        else 
          basis(:,5) = basis(:,6)
        end if
      end if
      basis(:,6) = 0.0_dp
    end if
    !#########################################################
    !## normalize rotational and translational eigenvectors ##
    !#########################################################
    do ii = 1,nindep
      rnorm= sqrt(dot_product(basis(:,ii),basis(:,ii)))
      basis(:,ii) = basis(:,ii)/rnorm
    end do
    do ii = 5, nindep
      do jj = 4, ii-1 
          prod = dot_product(basis(:,ii),basis(:,jj))
          call daxpy(n3,-prod,basis(1,jj),1,basis(1,ii),1)
      end do
      rnorm = sqrt(dot_product(basis(:,ii),basis(:,ii)))
      call dscal(n3,1.0_dp/rnorm,basis(1,ii),1)
    end do

    !###################################################
    !## add 3N-nindep orthonormal vibrational vectors ##
    !###################################################
    do ii = nindep+1, n3
20  continue
      call getRandom(randomFreq, basis(:,ii))
      basis(:,ii) = 100.0_dp*basis(:,ii)
      do jj=1, ii-1 
          prod = dot_product(basis(:,ii),basis(:,jj))
          call daxpy(n3,-prod,basis(1,jj),1,basis(1,ii),1)
      end do
      rnorm = sqrt(dot_product(basis(:,ii),basis(:,ii)))
      if( rnorm < 1.0E-2_dp) goto 20
      basis(:,ii) = basis(:,ii)/rnorm
    end do


    if(centinv) then
      do ii = nindep+1, n3 
        do jj = 1, n3 
          iCount = (jj+2)/3
          jCount = invpartner(iCount)
          kk   = jj-3*(iCount-1)+3*(jCount-1)
          if( jCount == 0 ) then
            newhess(jj,ii-nindep) = 0.0_dp
            dotpr(jj,ii-nindep) = basis(jj,ii)
          else if( iCount > jCount ) then
            xvalplus  = basis(jj,ii) - basis(kk,ii)
            xvalminus = basis(jj,ii) + basis(kk,ii)
            newhess(jj,ii-nindep) =  xvalplus
            newhess(kk,ii-nindep) = -xvalplus
            dotpr(jj,ii-nindep)   =  xvalminus
            dotpr(kk,ii-nindep)   =  xvalminus
          end if 
        end do
      end do

    !############################################
    !## project out translations and rotations ##
    !############################################
      do ii = 1, n3-nindep 
        do jj = 1, 3 
          prod = dot_product(dotpr(:,ii),basis(:,jj))
          call daxpy(n3,-prod,basis(1,jj),1,dotpr(1,ii),1)
        end do
        do jj = 4, nindep
          prod = dot_product(newhess(:,ii),basis(:,jj))
          call daxpy(n3,-prod,basis(1,jj),1,newhess(1,ii),1)
        end do
      end do 

!#################################################
!## construct "gerade" vibrational eigenvectors ##
!#################################################
      licz = nindep
      do ii = 1, n3-nindep
        rnorm = sqrt(dot_product(newhess(:,ii), newhess(:,ii)))
        if( rnorm > 1.0d-4) then
          call dscal(n3,1.0_dp/rnorm,newhess(1,ii),1)
          do jj = nindep+1, licz
            prod = dot_product(newhess(:,ii),basis(:,jj))
            call daxpy(n3,-prod,basis(1,jj),1,newhess(1,ii),1)
          end do
          rnorm = sqrt(dot_product(newhess(:,ii), newhess(:,ii)))
          if( rnorm>1.0d-2 ) then
            licz = licz+1
            newhess(:,ii) = newhess(:,ii)/rnorm
            basis(:,licz) = newhess(:,ii)
          end if
        end if
      end do
      ngerade = licz

!###################################################
!## construct "ungerade" vibrational eigenvectors ##
!###################################################
      do ii = 1, n3-nindep 
        rnorm =  sqrt(dot_product(dotpr(:,ii),dotpr(:,ii)))
        if(rnorm>1.0d-2) then
          call dscal(n3,1/rnorm,dotpr(1,ii),1)
          do jj = ngerade+1, licz
            prod = dot_product(dotpr(:,ii),basis(:,jj))
            call daxpy(n3,-prod,basis(1,jj),1,dotpr(1,ii),1)
          end do
          rnorm = sqrt(dot_product(dotpr(:,ii),dotpr(:,ii)))
          if(rnorm > 1.0d-2) then
            licz = licz+1
            call dscal(n3,1.0_dp/rnorm,dotpr(1,ii),1)
            basis(:,licz) = dotpr(:,ii)
          end if
        end if
      end do
    end if


!######################################################
!## transform hessian to the new ROT+TRANS+VIB basis ##
!######################################################
    call dgemm('T', 'N', n3, n3, n3, 1.0_dp, basis, n3, pDynMatrix, n3, 0.0_dp, newhess, n3)
    call dgemm('N', 'N', n3, n3, n3, 1.0_dp, newhess, n3, basis, n3, 0.0_dp, pDynMatrix, n3)
    newhess = pDynMatrix

    if (centinv) then
        newhess(nindep+1:ngerade,ngerade+1:n3) = 0.0_dp
        newhess(ngerade+1:n3,nindep+1:ngerade) = 0.0_dp
    end if

    eigenValues = 0.0_dp
    call heev(newhess(nindep+1:,nindep+1:), eigenValues(nindep+1:), 'U', 'V')


    pDynMatrix(:,1:nindep) = basis(:, 1:nindep)
    call dgemm('N', 'N', n3, n3-nindep, n3-nindep, 1.0_dp, basis(1,nindep+1), n3, newhess(nindep+1,nindep+1), n3, 0.0_dp, pDynMatrix(1,nindep+1), n3)


    do ii=1, nAtom
        massx=1.0_dp/sqrt(mass(ii))
        do jj =1, 3
            kk = jj+3*(ii-1)
            call dscal(n3,massx,pDynMatrix(kk,1),n3)
        end do
    end do
    do ii=1, n3
        redmas(ii)=1.0_dp/dot_product(pDynMatrix(:, ii), pDynMatrix(:,ii))
    end do
    do ii=1, n3
        pDynMatrix(:,ii) = pDynMatrix(:,ii)*sqrt(redmas(ii))
    end do

  ! take square root of modes (allowing for imaginary modes) and print
    eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)

!    write(*,*) 'Low frequencies:'
!    write(*,'(6f8.2)') eigenValues(1:nindep)*Hartree__cm

    write(stdOut,*)'Vibrational modes (cm-1):'
    do ii = nindep+1, 3*nAtom
      write(stdOut,'(f8.2)')eigenValues(ii)*Hartree__cm
    end do
    write(*,*)

    open(newunit=jj, file='FREQ.DAT', action="write", status="replace", form="formatted")
      write(jj,*)'Vibrational modes (cm-1):'
      do ii = nindep+1, 3*nAtom
        write(jj,'(f8.2)')eigenValues(ii)*Hartree__cm
      end do
      write(jj,*)
    close(jj)

    if (tNormalModes) then
        fdNormalModes = 201
        open(fdNormalModes, file=normalModesOut, STATUS='UNKNOWN')
        write(*,*) 'Writing normal modes in "vibrarions.molden"'


        write(fdNormalModes,*)"Entering Gaussian System"
        write(fdNormalModes,*)"                        Standard orientation:"
        write(fdNormalModes,*)"---------------------------------------------------------------------"
        write(fdNormalModes,*)"Center     Atomic     Atomic              Coordinates (Angstroms)"
        write(fdNormalModes,*)"Number     Number      Type              X           Y           Z   " 
        write(fdNormalModes,*)"---------------------------------------------------------------------"
        do ii=1, nAtom
            write(fdNormalModes,'(i5,6x,i5,13x,i1,4x,3f12.6)') ii, getAtomicNumberFromSymbol(speciesName(species(ii))),0,(0.529177208300000035*coordcm(jj,ii),jj=1,3)
        end do
        write(fdNormalModes,*)"---------------------------------------------------------------------"
        write(fdNormalModes,*)"  208 basis functions,   360 primitive gaussians,   234 cartesian basis functions"
        write(fdNormalModes,*)"   21 alpha electrons       21 beta electrons"
        write(fdNormalModes,*)
        write(fdNormalModes,*)"**********************************************************************"
        write(fdNormalModes,*) '' 
        write(fdNormalModes,*)"   and normal coordinates:"
        do ii=nindep+1, n3, 3
            kk=min(ii+2,n3)
            write(fdNormalModes,'(3(I22,1x))')(jj-nindep,jj=ii,kk)
            write(fdNormalModes,'(21x,3(a1,22x))')('A',jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Frequencies -- ",(eigenValues(jj)*Hartree__cm,jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Red. masses -- ",(redmas(jj),jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Frc consts  -- ",0.1,0.1,0.1
            write(fdNormalModes,'(a,3(f10.4,13x))')" IR Inten    -- ",(0.0_dp,jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Raman Activ -- ",(0.0_dp,jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Depolar (P) -- ",(0.0_dp,jj=ii,kk)
            write(fdNormalModes,'(a,3(f10.4,13x))')" Depolar (U) -- ",(0.0_dp,jj=ii,kk)
            write(fdNormalModes,*)"Atom AN      X      Y      Z        X      Y      Z        X      Y      Z"
            do jj=1, nAtom
                write(fdNormalModes,'(2i4,3(2x,3f7.2))')jj,getAtomicNumberFromSymbol(speciesName(species(jj))) ,((pDynMatrix(3*(jj-1)+licz,iCount),licz=1,3),iCount=ii,kk)
            end do
        end do
        write(fdNormalModes,*) '' 
        close(fdNormalModes)
    end if

    open(newunit=info, file=autotestTag, action="write", status="old", position="append")
    call writeTagged(info, tag_frequencies, eigenValues(nindep+1:))
    close(info)


    DEALLOCATE(eigenValues)
    DEALLOCATE(basis)
    DEALLOCATE(newhess)
    DEALLOCATE(coordcm)
    DEALLOCATE(invpartner)
    DEALLOCATE(dotpr)
    DEALLOCATE(redmas)
!    end if

!  if (tWriteDetailedOut) then
!    close(fdUser)
!  end if
