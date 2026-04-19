subroutine comenv_lambda(KW,M0,L,R,MENVD,LAMBDA,id,LAMBF)
    use track_support
    implicit none
    real(dp), intent(in):: M0,L,R,MENVD,LAMBDA
    integer, intent(in) :: KW
    real(dp), intent(out) :: LAMBF

    REAL(dp) :: CELAMF, RZAMSF
    EXTERNAL CELAMF, RZAMSF
    integer, intent(in), optional :: id
    real(dp):: RZAMS
    integer :: idd
    type(track), pointer :: t

    idd = 1
    if(present(id)) idd = id

    t => tarr(idd)

    if (t% is_he_track) then
        RZAMS = 10.d0**t% tr(i_logR,ZAMS_HE_EEP)
    elseif(kw>= He_MS .and. kw<=He_GB .and. use_sse_NHe)then
        RZAMS = RZAMSF(M0)
    else
        RZAMS = 10.d0**t% tr(i_logR,ZAMS_EEP)
    endif

    if (LAMBDA < -1.5d0 .and. i_binding_energy_re > 0) then
        LAMBF = - (t%pars%mass * (t%pars%mass - t%pars%core_mass)) / &
         ((t%pars%binding_energy_re / 3.8d48) * R)
        write(0,'(A,I3,A,F8.3,A,F8.3,A,F8.3,A,ES12.4,A,ES12.4,A,F10.5)') &
            '[CE lambda] KW=',KW,' M=',t%pars%mass,' Mc=',t%pars%core_mass, &
            ' R=',R,' BE_re=',t%pars%binding_energy_re, &
            ' BE=',t%pars%binding_energy,' LAMBF(be_re)=',LAMBF
    else if (LAMBDA < -0.5d0 .and. i_binding_energy > 0) then
        LAMBF = - (t%pars%mass * (t%pars%mass - t%pars%core_mass)) / &
         ((t%pars%binding_energy / 3.8d48) * R)
        write(0,'(A,I3,A,F8.3,A,F8.3,A,F8.3,A,ES12.4,A,F10.5)') &
            '[CE lambda] KW=',KW,' M=',t%pars%mass,' Mc=',t%pars%core_mass, &
            ' R=',R,' BE=',t%pars%binding_energy,' LAMBF(be)=',LAMBF
    else
        LAMBF = CELAMF(KW,M0,L,R,RZAMS,MENVD,LAMBDA)
        write(0,'(A,I3,A,F8.3,A,F8.3,A,F10.5)') &
            '[CE lambda] KW=',KW,' M=',M0,' R=',R, &
            ' RZAMS=',RZAMS,' LAMBF(analytic)=',LAMBF
    endif

    LAMBF = MIN(100.0d0, LAMBF)

    nullify(t)
end subroutine comenv_lambda