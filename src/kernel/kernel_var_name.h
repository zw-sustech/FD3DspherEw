!-------------------------------------------------------------------------------
!-- kernel var name in nc files
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-- isotropic medium
!-------------------------------------------------------------------------------
!-- lambda+mu
    character (len=*),parameter ::                     &
        VNAME_KER_PHA_LA_LM = 'phase_lambda_LM'     ,  & 
        VNAME_KER_PHA_MU_LM = 'phase_mu_LM'         ,  & 
        VNAME_KER_AMP_LA_LM = 'amplitude_lambda_LM' ,  & 
        VNAME_KER_AMP_MU_LM = 'amplitude_mu_LM'
!-- kappa+mu
    character (len=*),parameter ::                    &
        VNAME_KER_PHA_KA_KM = 'phase_kappa_KM'     ,  & 
        VNAME_KER_PHA_MU_KM = 'phase_mu_KM'        ,  & 
        VNAME_KER_AMP_KA_KM = 'amplitude_kappa_KM' ,  & 
        VNAME_KER_AMP_MU_KM = 'amplitude_mu_KM'
!-- Vp+Vs
    character (len=*),parameter ::                       &
        VNAME_KER_PHA_VP         = 'phase_Vp'         ,  &
        VNAME_KER_PHA_VS         = 'phase_Vs'         ,  &
        VNAME_KER_AMP_VP         = 'amplitude_Vp'     ,  &
        VNAME_KER_AMP_VS         = 'amplitude_Vs'

!-------------------------------------------------------------------------------
!-- VTI medium
!-------------------------------------------------------------------------------
    character (len=*),parameter ::                                &
        VNAME_KER_PHA_VTI_VPH        = 'phase_Vph'           ,  &
        VNAME_KER_PHA_VTI_VPV        = 'phase_Vpv'           ,  &
        VNAME_KER_PHA_VTI_VSH        = 'phase_Vsh'           ,  &
        VNAME_KER_PHA_VTI_VSV        = 'phase_Vsv'           ,  &
        VNAME_KER_PHA_VTI_ETA        = 'phase_eta'           ,  &
        VNAME_KER_AMP_VTI_VPH    = 'amplitude_Vsh'       ,  &
        VNAME_KER_AMP_VTI_VPV    = 'amplitude_Vph'       ,  &
        VNAME_KER_AMP_VTI_VSH    = 'amplitude_Vpv'       ,  &
        VNAME_KER_AMP_VTI_VSV    = 'amplitude_Vsv'       ,  &
        VNAME_KER_AMP_VTI_ETA    = 'amplitude_eta'       ,  &
        VNAME_KER_PHA_VTI_A          = 'phase_A'             ,  &
        VNAME_KER_PHA_VTI_C          = 'phase_C'             ,  &
        VNAME_KER_PHA_VTI_N          = 'phase_N'             ,  &
        VNAME_KER_PHA_VTI_L          = 'phase_L'             ,  &
        VNAME_KER_PHA_VTI_F          = 'phase_F'             ,  &
        VNAME_KER_AMP_VTI_A      = 'amplitude_A'         ,  &
        VNAME_KER_AMP_VTI_C      = 'amplitude_C'         ,  &
        VNAME_KER_AMP_VTI_N      = 'amplitude_N'         ,  &
        VNAME_KER_AMP_VTI_L      = 'amplitude_L'         ,  &
        VNAME_KER_AMP_VTI_F      = 'amplitude_F'

!-------------------------------------------------------------------------------
!-- Generial ANISO medium
!-------------------------------------------------------------------------------
    character (len=*),parameter ::                                &
        VNAME_KER_PHA_C11        = 'phase_C11'               ,  &
        VNAME_KER_PHA_C12        = 'phase_C12'               ,  &
        VNAME_KER_PHA_C13        = 'phase_C13'               ,  &
        VNAME_KER_PHA_C14        = 'phase_C14'               ,  &
        VNAME_KER_PHA_C15        = 'phase_C15'               ,  &
        VNAME_KER_PHA_C16        = 'phase_C16'               ,  &
        VNAME_KER_PHA_C22        = 'phase_C22'               ,  &
        VNAME_KER_PHA_C23        = 'phase_C23'               ,  &
        VNAME_KER_PHA_C24        = 'phase_C24'               ,  &
        VNAME_KER_PHA_C25        = 'phase_C25'               ,  &
        VNAME_KER_PHA_C26        = 'phase_C26'               ,  &
        VNAME_KER_PHA_C33        = 'phase_C33'               ,  &
        VNAME_KER_PHA_C34        = 'phase_C34'               ,  &
        VNAME_KER_PHA_C35        = 'phase_C35'               ,  &
        VNAME_KER_PHA_C36        = 'phase_C36'               ,  &
        VNAME_KER_PHA_C44        = 'phase_C44'               ,  &
        VNAME_KER_PHA_C45        = 'phase_C45'               ,  &
        VNAME_KER_PHA_C46        = 'phase_C46'               ,  &
        VNAME_KER_PHA_C55        = 'phase_C55'               ,  &
        VNAME_KER_PHA_C56        = 'phase_C56'               ,  &
        VNAME_KER_PHA_C66        = 'phase_C66'               ,  &
        VNAME_KER_AMP_C11        = 'amplitude_C11'               ,  &
        VNAME_KER_AMP_C12        = 'amplitude_C12'               ,  &
        VNAME_KER_AMP_C13        = 'amplitude_C13'               ,  &
        VNAME_KER_AMP_C14        = 'amplitude_C14'               ,  &
        VNAME_KER_AMP_C15        = 'amplitude_C15'               ,  &
        VNAME_KER_AMP_C16        = 'amplitude_C16'               ,  &
        VNAME_KER_AMP_C22        = 'amplitude_C22'               ,  &
        VNAME_KER_AMP_C23        = 'amplitude_C23'               ,  &
        VNAME_KER_AMP_C24        = 'amplitude_C24'               ,  &
        VNAME_KER_AMP_C25        = 'amplitude_C25'               ,  &
        VNAME_KER_AMP_C26        = 'amplitude_C26'               ,  &
        VNAME_KER_AMP_C33        = 'amplitude_C33'               ,  &
        VNAME_KER_AMP_C34        = 'amplitude_C34'               ,  &
        VNAME_KER_AMP_C35        = 'amplitude_C35'               ,  &
        VNAME_KER_AMP_C36        = 'amplitude_C36'               ,  &
        VNAME_KER_AMP_C44        = 'amplitude_C44'               ,  &
        VNAME_KER_AMP_C45        = 'amplitude_C45'               ,  &
        VNAME_KER_AMP_C46        = 'amplitude_C46'               ,  &
        VNAME_KER_AMP_C55        = 'amplitude_C55'               ,  &
        VNAME_KER_AMP_C56        = 'amplitude_C56'               ,  &
        VNAME_KER_AMP_C66        = 'amplitude_C66'

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
