!             Program block renorm-main.f    
!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO 
!             E-MAIL:   mhjensen@fys.uio.no
!             LANGUAGE: F95 
!
!             Program to renormalize the nucleon-nucleon interaction
!             using 
!             1)  a g-matrix technique via matrix inversion. It is energy dependent. 
!                 option in renorm.ini  :   g-matrix
!             2)  a no-core technique via a similarity transformation.
!                 Sharp cutoff in oscillator space
!                 option in renorm.ini  :   no-core
!             3)  a momentum space renormalization via similarity transformation. 
!                 Sharp cutoff in momentum space
!                 option in renorm.ini  :   vlowk
!             4)  Flow renormalization group (rg) equations in momentum space  (k)
!                 Soft cutoff in momentum space
!                 option in renorm.ini  :   v-krg
!             5)  Flow renormalization group (rg) equations in oscillator space  (n)
!                 option in renorm.ini  :   v-nrg
!                 Soft cutoff in oscillator space
!             6)  You can also produce just the bare interaction, using the vbare option
!
!             Options 2 and 5 can also be used to derive a renormalized Coulomb interaction
!             for spherical quantum dots in three dimensions
!
!             Run simply as vrenorm.exe, it picks up the input file renorm.ini
!     Main program starts here
!
PROGRAM  v_effective
  USE ang_mom_functions
  USE constants
  USE inifile
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE
 ! Oldfashioned  common blocks   for Bonn-type potentials 
  INTEGER ::  kread,kwrite,kpunch,kda
  COMMON /crdwrt/ kread,kwrite,kpunch,kda(9)
  CHARACTER (LEN=120) :: infilename, outputfile, int_Mscheme
  CHARACTER (LEN=120) :: int_Jscheme
  CHARACTER (LEN=120) :: Jscheme_spdata, Mscheme_spdata
  LOGICAL :: fail

  !     Read in and organize all relevant data, the file renorm.ini contains all data
  infilename   = 'renorm.ini'
  call ini_open(infilename, 5, fail, .false.)
  IF (fail) STOP 'Error opening parameter file, probably wrong file, use renorm.ini as name'
  !     open the output files
  outputfile = ini_read_string('output_run')
  OPEN(unit=6,file=outputfile)
  kread = 14; kwrite = 6
  int_Mscheme = ini_read_string('int_Mscheme')
  OPEN(UNIT=14,file=int_Mscheme)
  !     file for the renormalized interaction, unformatted output
  int_Jscheme = Ini_Read_String('int_Jscheme')
! in case of binary  OPEN(UNIT=7,FILE=renorminteraction_file,form='binary')
  OPEN(UNIT=7,FILE=int_Jscheme)
  !     file for single-particle data and calculation information
  Jscheme_spdata = Ini_Read_String('Jscheme_spdata')
  Mscheme_spdata = Ini_Read_String('Mscheme_spdata')
  OPEN(UNIT=8,FILE=Jscheme_spdata)
  OPEN(UNIT=9,FILE=Mscheme_spdata)
  type_of_renormv = ini_read_string('type_of_renormv')
  coulomb_included = ini_read_string('coulomb_included')
  physical_system = ini_read_string('physical_system')
  !  Setup data and start calculation
  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     CALL setup_nuclearsp_data
     CALL setup_nuclearveff_data
  CASE('atomic_physics')
     CALL setup_atomicsp_data
     CALL setup_atomicveff_data
  END SELECT
  !     Factorials for 3j, 6j and 9j symbols and for the moshinsky transf coeffs
  CALL commons_to_angmom  
  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     !     setup of the effective interactions
     SELECT CASE (type_of_renormv) 
     CASE('g-matrix')
        !    Setup max number of mesh points needed
        CALL setup_ho_cutoff
        CALL setup_g_matrix
     CASE('no-core')
        !    Setup max number of mesh points needed
        CALL setup_ho_cutoff
        CALL setup_nocore
     CASE('vlowk')
        !  for vlowk this is fixed by input variables
        n_rel = n_k1+n_k2
        CALL setup_vlowk
     CASE('vbare')
        !  for vlowk this is fixed by input variables
        n_rel = n_k1+n_k2
        CALL setup_vlowk
     CASE('v-krg')
        !  mesh points for solving differential equation
        n_rel = n_k1+n_k2
        CALL setup_vkrg
     CASE('v-nrg')
        !  here we use the no-core infrastructure in an oscillator basis
        CALL setup_ho_cutoff
        CALL setup_nocore
     END SELECT
  CASE('atomic_physics')
     CALL setup_nocore
  END SELECT
! Prepare for transformation to m-scheme, first sp basis
  CALL mscheme_orbits
! read in two-body matrix elements in J-scheme
  CALL read_gmatrix
! Then convert two-body matrix elements to m-scheme
  CALL convert_to_mscheme

END PROGRAM v_effective
!
!     This subroutine reads in variables needed to specify the
!     effective interaction 
!     Read in energy variables for the g-matrix, effective interaction
!     and restriction on excitations in oscillator energy
!
SUBROUTINE setup_nuclearveff_data
  USE constants
  USE inifile
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE 
  INTEGER :: i!, number_twobody_elements,number_pp_elements, number_pn_elements, number_nn_elements 
  REAL(DP) :: interval
  !  Extract variables
  type_of_pot = Ini_Read_String('type_of_pot')
  csb = Ini_Read_String('csb_choice')
  cib = Ini_Read_String('cib_choice')
  SELECT CASE (type_of_renormv) 
  CASE('g-matrix')
     n_startenergy_g = Ini_Read_Int('n_startenergy_g')
     IF ( n_startenergy_g < 0) THEN
        WRITE(6,*) 'Number of starting energies less than zero'
        STOP 'startingenergygwrong'
     ENDIF
     !  always odd number of starting energies
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) )
     e_start_g(1) = Ini_Read_Double('first_startingenergy')
     e_start_g(n_startenergy_g) = Ini_Read_Double('last_startingenergy')
     interval = (e_start_g(n_startenergy_g)-e_start_g(1))/n_startenergy_g
     DO i = 2, n_startenergy_g-1
        e_start_g(i)=e_start_g(i-1)+interval
     ENDDO 
  CASE('no-core')
     n_startenergy_g = 0
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) ) 
     e_start_g = 0.0_dp
  CASE('v-nrg')
     n_startenergy_g = 0
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) ) 
     e_start_g = 0.0_dp
  CASE('vlowk')
     n_startenergy_g = 0
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) ) 
     e_start_g = 0.0_dp
     n_k1 = Ini_Read_Int('n_k1')
     n_k2 = Ini_Read_Int('n_k2')
     n_k = n_k1+n_k2
     k_cutoff = Ini_Read_Double('k_cutoff')
     k_max = Ini_Read_Double('k_max')
        IF ((type_of_pot == 'OPEP').OR.(type_of_pot == 'Tensorinteraction')   &
             .OR.(type_of_pot == 'LSinteraction') ) THEN
           n_k2 = 0
           n_k1 = n_rel
           n_rel = n_k1
           k_cutoff = k_max
        ENDIF
  CASE('vbare')
     n_startenergy_g = 0
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) ) 
     e_start_g = 0.0_dp
     n_k1 = Ini_Read_Int('n_k1')
     n_k2 = Ini_Read_Int('n_k2')
     n_k = n_k1+n_k2
     k_cutoff = Ini_Read_Double('k_cutoff')
     k_max = Ini_Read_Double('k_max')
  CASE('v-krg')
     n_startenergy_g = 0
     IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
     ALLOCATE ( e_start_g (n_startenergy_g) ) 
     e_start_g = 0.0_dp
     n_k1 = Ini_Read_Int('n_k1')
     n_k2 = Ini_Read_Int('n_k2')
     n_k = n_k1+n_k2
     k_cutoff = Ini_Read_Double('k_cutoff')
     k_max = Ini_Read_Double('k_max')
  END SELECT

  !  find total number of two-body matrix elements
  CALL find_number_twobody!(number_twobody_elements, number_pp_elements, &
!       number_pn_elements, number_nn_elements )  

END SUBROUTINE setup_nuclearveff_data
!
!     Reads in and allocates sp data for the nuclear case
!     Read comments below in order to properly understand the input
! 
SUBROUTINE setup_nuclearsp_data
  USE single_particle_orbits  
  USE inifile
  USE constants
  IMPLICIT NONE
  INTEGER :: i, n1, l1, j1, t1, shell1, nnn, nn1, ll1, lll, j, l12
  CHARACTER(LEN= 10) space, model
  INTEGER :: model_space, states, max_space
  INTEGER, ALLOCATABLE, DIMENSION(:) ::  n_big,l_big,j_big,tzp_big, shell_big
  CHARACTER(LEN= 10), ALLOCATABLE, DIMENSION(:) :: space_big, model_big

  !   Extract values for various variables
  hbar_omega = Ini_Read_Double('hbar_omega')
  oscl=hbarc/SQRT(p_massave*hbar_omega)
  pauli_operator = Ini_Read_String('pauli_operator')
  max_space = Ini_Read_Int('max_space')
  lab_lmax = Ini_Read_Int('lab_lmax')
  jmin = Ini_Read_Int('jmin')
  jmax = Ini_Read_Int('jmax')
  mass_nucleus = Ini_Read_Int('mass_nucleus')
  lab_nmax = Ini_Read_Int('lab_nmax')
  WRITE(8,'(68H   ----> Oscillator parameters, Model space and single-particle data)')
  WRITE(8,'(65HMass number A of chosen nucleus (important for CoM corrections): ,I10)') mass_nucleus 
  WRITE(8,'(30HOscillator length and energy: ,E12.6,2X,E12.6)')  oscl, hbar_omega
  ! Here we read the modelspace max l and n and the 2n+l for the last shell
  square_calculation = .FALSE.

  SELECT CASE (type_of_renormv) 
  CASE('g-matrix')
     SELECT CASE (pauli_operator)
     CASE ('square')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nmax  = model_space
        nlmax = nlmax_model
        square_calculation = .TRUE.
     CASE ('triangular')
        model_space = 2*lab_nmax+lab_lmax
        nlmax_model = model_space
        nlmax = 2*lab_nmax+lab_lmax
        nmax = lab_nmax
     CASE ('wings')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nlmax = 2*lab_nmax   !  we make a triangular truncation here.
        nmax = lab_nmax
     END SELECT
  CASE('no-core')
     model_space = 2*lab_nmax+lab_lmax
     nlmax_model = model_space
     nlmax = max_space
     nmax = max_space/2
  CASE('v-nrg')
     model_space = 2*lab_nmax+lab_lmax
     nlmax_model = model_space
     nlmax = max_space
     nmax = max_space/2
  CASE('vlowk')
     SELECT CASE (pauli_operator)
     CASE ('square')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nmax  = model_space
        nlmax = nlmax_model
        square_calculation = .TRUE.
     CASE ('triangular')
        model_space = 2*lab_nmax+lab_lmax
        nlmax_model = model_space
        nlmax = 2*lab_nmax+lab_lmax
        nmax = lab_nmax
     CASE ('wings')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nlmax = 2*lab_nmax   !  we make a triangular truncation here.
        nmax = lab_nmax
     END SELECT
  CASE('vbare')
     SELECT CASE (pauli_operator)
     CASE ('square')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nmax  = model_space
        nlmax = nlmax_model
        square_calculation = .TRUE.
     CASE ('triangular')
        model_space = 2*lab_nmax+lab_lmax
        nlmax_model = model_space
        nlmax = 2*lab_nmax+lab_lmax
        nmax = lab_nmax
     CASE ('wings')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nlmax = 2*lab_nmax   !  we make a triangular truncation here.
        nmax = lab_nmax
     CASE ('talent')
        model_space = 2*lab_nmax+lab_lmax
        nlmax_model = lab_lmax
        nlmax = model_space 
        nmax = lab_nmax
     END SELECT
  CASE('v-krg')
     SELECT CASE (pauli_operator)
     CASE ('square')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nmax  = model_space
        nlmax = nlmax_model
        square_calculation = .TRUE.
     CASE ('triangular')
        model_space = 2*lab_nmax+lab_lmax
        nlmax_model = model_space
        nlmax = 2*lab_nmax+lab_lmax
        nmax = lab_nmax
     CASE ('wings')
        model_space = lab_lmax 
        nlmax_model = 2*model_space
        nlmax = 2*lab_nmax   !  we make a triangular truncation here.
        nmax = lab_nmax
     END SELECT
  END SELECT
  !  Relative lmax
  IF ( lmax < jmax+1) lmax = jmax+1
  WRITE(8,*) 'Min and max value of partial wave ang. mom', jmin, jmax
  WRITE(8,*) 'Max value of relative orb mom or cm orb mom,  l or L= ', lmax
  WRITE(8,*) 'Max value of relative n:', nmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for large space:', nlmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for model space:', nlmax_model
  CALL count_single_particle_orbits(states)
  neutron_data%total_orbits = states
  proton_data%total_orbits = states
  ALLOCATE( n_big(2*states), l_big(2*states), j_big(2*states), tzp_big(2*states), shell_big(2*states))
  ALLOCATE( space_big(2*states), model_big(2*states))
  !     Setup all possible orbit information
  all_orbit%total_orbits=neutron_data%total_orbits+ proton_data%total_orbits
  !     Allocate space for all single-particle data
  CALL allocate_sp_array(neutron_data,neutron_data%total_orbits) 
  CALL allocate_sp_array(proton_data,proton_data%total_orbits) 
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  !     Read neutron single-particle data
  CALL setup_single_particle_orbits(model_space, states, &
       n_big,l_big,j_big,tzp_big, shell_big, space_big, model_big)
  DO i=1, neutron_data%total_orbits 
     neutron_data%nn(i)=n_big(2*i)
     neutron_data%ll(i)=l_big(2*i); neutron_data%jj(i)=j_big(2*i)
     neutron_data%itzp(i)=tzp_big(2*i)
     neutron_data%nshell(i)=shell_big(2*i)
     neutron_data%orbit_status(i)=space_big(2*i) 
     neutron_data%model_space(i)=model_big(2*i)
     neutron_data%e(i)=hbar_omega*(shell_big(2*i)+1.5)
     neutron_data%evalence(i)= 0.0_dp
     !     Neutrons are in the internal structure always even numbers 
     all_orbit%nn(i*2)=neutron_data%nn(i)
     all_orbit%ll(i*2)=neutron_data%ll(i)
     all_orbit%jj(i*2)=neutron_data%jj(i)
     all_orbit%nshell(i*2)=neutron_data%nshell(i)
     all_orbit%e(i*2)=neutron_data%e(i)
     all_orbit%evalence(i*2) = neutron_data%evalence(i)
     all_orbit%itzp(i*2)=neutron_data%itzp(i)
     all_orbit%orbit_status(i*2)=neutron_data%orbit_status(i)
     all_orbit%model_space(i*2)=neutron_data%model_space(i)
  ENDDO
  DO i=1,  proton_data%total_orbits 
     proton_data%nn(i)=n_big(2*i-1)
     proton_data%ll(i)=l_big(2*i-1)
     proton_data%jj(i)=j_big(2*i-1)
     proton_data%itzp(i)=tzp_big(2*i-1)
     proton_data%nshell(i)=shell_big(2*i-1)
     proton_data%orbit_status(i)=space_big(2*i-1) 
     proton_data%model_space(i)=model_big(2*i-1)
     proton_data%evalence(i)  = 0.0_dp
     proton_data%e(i)=hbar_omega*(shell_big(2*i-1)+1.5)
     !     protons are in the internal structure always odd numbers
     all_orbit%nn(i*2-1)=proton_data%nn(i)
     all_orbit%ll(i*2-1)=proton_data%ll(i)
     all_orbit%jj(i*2-1)=proton_data%jj(i)
     all_orbit%nshell(i*2-1)=proton_data%nshell(i)
     all_orbit%evalence(i*2-1) = proton_data%evalence(i)
     all_orbit%e(i*2-1)=proton_data%e(i)
     all_orbit%itzp(i*2-1)=proton_data%itzp(i)
     all_orbit%orbit_status(i*2-1)=proton_data%orbit_status(i)
     all_orbit%model_space(i*2-1)=proton_data%model_space(i)
  ENDDO
  WRITE(8,*) 'Total number of single-particle orbits', all_orbit%total_orbits 
!  WRITE(8,'(72HLegend:         n   l  2j  tz    2n+l  energy evalence  particle or hole)')
  WRITE(8,'(101HLegend:         n     l     2j   tz    2n+l  HO-energy     evalence     particle/hole  inside/outside)') 
  j_lab_max = 0
  DO i=1, all_orbit%total_orbits
     WRITE(8,'(7HNumber:,6(I4,2X),2X,E12.6,2X,E12.6,2X,A10,2X,A10)') i, all_orbit%nn(i), all_orbit%ll(i), &
          all_orbit%jj(i), all_orbit%itzp(i), all_orbit%nshell(i), all_orbit%e(i), all_orbit%evalence(i), &
          all_orbit%orbit_status(i), all_orbit%model_space(i)
     IF ( all_orbit%jj(i) > j_lab_max) j_lab_max = all_orbit%jj(i)
  ENDDO
  ! Here we fix itzmin, itzmax and j_lab_min
  itzmin = -1; itzmax=1; j_lab_min=0
  DEALLOCATE( space_big, model_big); DEALLOCATE( n_big, l_big, j_big, tzp_big, shell_big)

END SUBROUTINE setup_nuclearsp_data
!
!     This subroutine reads in variables needed to specify the
!     effective interaction for atomic systems, viz 3dim spherical electrons
!     For atomic systems, only nocore and vnrg are the options for renormalizing the
!     Coulomb
!
SUBROUTINE setup_atomicveff_data
  USE constants
  USE inifile
  USE wave_functions
  USE stored_bare_interaction
  IMPLICIT NONE 
  INTEGER :: i!, number_twobody_elements,number_pp_elements, number_pn_elements, number_nn_elements 
  CHARACTER (LEN=120) :: renorminteraction_file
  REAL(DP) :: interval
  !  Extract variables
  renorminteraction_file = Ini_Read_String('renorminteraction_file')
  OPEN(UNIT=7,FILE=renorminteraction_file)
  n_startenergy_g = 0
  IF ( MOD(n_startenergy_g,2) == 0) n_startenergy_g = n_startenergy_g+1
  ALLOCATE ( e_start_g (n_startenergy_g) ) 
  e_start_g = 0.0_dp
  !  find total number of two-body matrix elements
  CALL find_number_twobody!(number_twobody_elements, number_pp_elements, &
!       number_pn_elements, number_nn_elements )  
  number_nn_elements = number_twobody_elements
  WRITE(7,'(34H   ----> Interaction part         )')
  WRITE(7,'(30HElectron-electron interaction:,A20)')
  WRITE(7,'(20HType of calculation:,A20)') type_of_renormv
  WRITE(7,'(38HNumber and value of starting energies:,I4)')  n_startenergy_g  
  WRITE(7,'(10(1X,E12.6) )') (e_start_g(i), i=1,n_startenergy_g )
  WRITE(7,'(38HTotal number of twobody matx elements:,4I15 )') number_twobody_elements, number_pp_elements, &
       number_pn_elements, number_nn_elements 
  WRITE(7,'(91HMatrix elements with the following legend, NOTE no hbar_omega/A for Hcom, p_ip_j and r_ir_j)')
  WRITE(7,'(103H  Tz Par  2J   a   b   c   d          <ab|V|cd>        <ab|Hcom|cd>     <ab|r_ir_j|cd>   <ab|p_ip_j|cd>)') 

END SUBROUTINE setup_atomicveff_data
!
!     Reads in and allocates sp data for the nuclear case
!     Read comments below in order to properly understand the input
! 
SUBROUTINE setup_atomicsp_data
  USE single_particle_orbits  
  USE inifile
  USE constants
  IMPLICIT NONE
  INTEGER :: i, n1, l1, j1, t1, shell1, nnn, nn1, ll1, lll, j, l12
  CHARACTER(LEN= 10) space, model
  INTEGER :: model_space, states, max_space
  INTEGER, ALLOCATABLE, DIMENSION(:) ::  n_big,l_big,j_big,tzp_big, shell_big
  CHARACTER(LEN= 10), ALLOCATABLE, DIMENSION(:) :: space_big, model_big
  CHARACTER (LEN=120) :: spdata_file

  !   Extract values for various variables
  spdata_file = Ini_Read_String('spdata_file')
  OPEN(UNIT=8,FILE=spdata_file)
  !   we use atomic units for electrons!
  hbar_omega = Ini_Read_Double('hbar_omega')
  atomic_strength =  Ini_Read_Double('lambda_value')
  oscl  = 1.0/SQRT(hbar_omega)
  pauli_operator = Ini_Read_String('pauli_operator')
  max_space = Ini_Read_Int('max_space')
  lab_lmax = Ini_Read_Int('lab_lmax')
  jmin = Ini_Read_Int('jmin')
  jmax = Ini_Read_Int('jmax')
  number_electrons = Ini_Read_Int('number_electrons')
  lab_nmax = Ini_Read_Int('lab_nmax')
  WRITE(8,'(68H   ----> Oscillator parameters, Model space and single-particle data)')
  WRITE(8,'(65HNumber of electrons (important for CoM corrections):             ,I10)') number_electrons
  WRITE(8,'(30HOscillator length and energy: ,E12.6,2X,E12.6)')  oscl, hbar_omega
  ! Here we read the modelspace max l and n and the 2n+l for the last shell
  square_calculation = .FALSE.
  model_space = 2*lab_nmax+lab_lmax
  nlmax_model = model_space
  nlmax = max_space
  nmax = max_space/2
  !  Relative lmax
  IF ( lmax < jmax+1) lmax = jmax+1
  WRITE(8,*) 'Min and max value of partial wave ang. mom', jmin, jmax
  WRITE(8,*) 'Max value of relative orb mom or cm orb mom,  l or L= ', lmax
  WRITE(8,*) 'Max value of relative n:', nmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for large space:', nlmax
  WRITE(8,*) 'Max value of 2*n + l+ cm 2*N +L for model space:', nlmax_model
  CALL count_single_particle_orbits(states)
  all_orbit%total_orbits = states
  ALLOCATE( n_big(2*states), l_big(2*states), j_big(2*states), tzp_big(2*states), shell_big(2*states))
  ALLOCATE( space_big(2*states), model_big(2*states))
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  CALL setup_single_particle_orbits(model_space, states, &
       n_big,l_big,j_big,tzp_big, shell_big, space_big, model_big)
  DO i=1, all_orbit%total_orbits 
     all_orbit%nn(i)=n_big(i)
     all_orbit%ll(i)=l_big(i)
     all_orbit%jj(i)=j_big(i)
     all_orbit%itzp(i)=tzp_big(i)
     all_orbit%nshell(i)=shell_big(i)
     all_orbit%orbit_status(i)=space_big(i) 
     all_orbit%model_space(i)=model_big(i)
     all_orbit%e(i) = (shell_big(i)+1.5)
     all_orbit%evalence(i)= 0.0_dp
  ENDDO
  WRITE(8,*) 'Total number of single-particle orbits', all_orbit%total_orbits 
  !  WRITE(8,'(72HLegend:         n   l  2j  tz    2n+l  energy evalence  particle or hole)')
  WRITE(8,'(101HLegend:         n     l     2j   tz    2n+l  HO-energy     evalence     particle/hole  inside/outside)') 
  j_lab_max = 0
  DO i=1, all_orbit%total_orbits
     WRITE(8,'(7HNumber:,6(I4,2X),2X,E12.6,2X,E12.6,2X,A10,2X,A10)') i, all_orbit%nn(i), all_orbit%ll(i), &
          all_orbit%jj(i), all_orbit%itzp(i), all_orbit%nshell(i), all_orbit%e(i), all_orbit%evalence(i), &
          all_orbit%orbit_status(i), all_orbit%model_space(i)
     IF ( all_orbit%jj(i) > j_lab_max) j_lab_max = all_orbit%jj(i)
  ENDDO
  ! Here we fix itzmin, itzmax and j_lab_min, itzmin = itzmax = 1 for electrons only (artificial)
  itzmin = 1; itzmax=1; j_lab_min=0
  DEALLOCATE( space_big, model_big); DEALLOCATE( n_big, l_big, j_big, tzp_big, shell_big)

END SUBROUTINE setup_atomicsp_data
!
!    setup of single-particle basis
!
SUBROUTINE setup_single_particle_orbits(model_space, norb, nn, ll, jj, itzp, nshell, space, model)
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: norb, model_space
  INTEGER, DIMENSION(norb) :: na,la,ja
  INTEGER, DIMENSION(2*norb), INTENT(INOUT) ::  nn,ll,jj,itzp, nshell
  CHARACTER(LEN= 10), DIMENSION(2*norb), INTENT(INOUT)::   space, model
  INTEGER :: icount, j_max, j_min, j, n, l

  icount = 0
  DO n = 0, lab_nmax
     DO l = 0, lab_lmax
        IF ( n+n+l > nlmax) CYCLE
        IF ( (square_calculation) .AND. ( n+n+l > nlmax_model/2) ) CYCLE
        j_min = l+l-1; j_max = l+l+1
        IF (j_min < 0) j_min = 1
        DO j = j_max, j_min, -2
           icount = icount + 1
           ja(icount)=j
           la(icount)=l
           na(icount)=n
        ENDDO
     ENDDO
  ENDDO
  IF ( icount /= norb) THEN
     WRITE(6,*) 'Error in allocation of single-particle states'
     STOP 'setup_single_particle_orbits'
  ENDIF
  ! now sort the sp orbits with increasing 2n+l and largest j value first
  CALL sort_spbasis(icount,na,la,ja)

  SELECT CASE (physical_system) 
  CASE('nuclear_physics')
     DO j = 1, icount
        jj(j*2)=ja(j)                     ! neutrons are 
        ll(j*2)=la(j)                     ! always even numbers
        nn(j*2)=na(j)                     ! protons are odd
        nshell(j*2)=2*na(j)+la(j)
        space(j*2)='particle'
        model(j*2)='inside'
        jj(j*2-1)=ja(j)                    
        ll(j*2-1)=la(j)
        nn(j*2-1)=na(j)
        nshell(j*2-1)=2*na(j)+la(j)
        space(j*2-1)='particle'
        model(j*2-1)='inside'
        itzp(j*2-1)=-1                   ! protons, nucl. phys. def.
        itzp(j*2)=1                      ! Neutrons, nucl. phys. def.
     ENDDO
  CASE('atomic_physics')
     DO j = 1, icount
        jj(j)=ja(j)                  
        ll(j)=la(j)                  
        nn(j)=na(j)                  
        nshell(j)=2*na(j)+la(j)
        space(j)='particle'
        model(j)='inside'
        itzp(j)=1                    
     ENDDO
  END SELECT

END SUBROUTINE setup_single_particle_orbits
!
!   Find maximum number of single-particle states constrained by n, l and j
!
SUBROUTINE count_single_particle_orbits(icount)
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: icount
  INTEGER :: n,l, j, j_max, j_min

  icount = 0
  DO n = 0, lab_nmax
     DO l = 0, lab_lmax
        IF ( n+n+l > nlmax) CYCLE
        IF ( (square_calculation) .AND. ( n+n+l > nlmax_model/2) ) CYCLE
        j_min = l+l-1; j_max = l+l+1
        IF (j_min < 0) j_min = 1
        DO j = j_max, j_min, -2
           icount = icount + 1
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE count_single_particle_orbits
!
!           Sort the sp basis according to shells and increasing j-value
!
SUBROUTINE sort_spbasis(count,na,la,ja)
  USE constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: count 
  INTEGER, DIMENSION(count), INTENT(INOUT) :: na,la,ja
  INTEGER :: i, j, t1, t2, t3

  !  sort according to shell
  DO i=1, count
     DO j=i, count
        IF( na(i)+na(i)+la(i) > na(j)+na(j)+la(j) ) THEN 
           t1 = na(i); t2 = la(i); t3 = ja(i)
           na(i) = na(j);   la(i) = la(j);     ja(i) =  ja(j)
           na(j) = t1;   la(j) =  t2;      ja(j) = t3
        ENDIF
     ENDDO
  ENDDO
  !  then sort according to ascending j-value
  DO i=1, count
     DO j=i, count
        IF( na(i)+na(i)+la(i) ==  na(j)+na(j)+la(j) ) THEN 
           IF ( ja(i) < ja(j) ) THEN
              t1 = na(i); t2 = la(i); t3 = ja(i)
              na(i) = na(j);   la(i) = la(j);     ja(i) =  ja(j)
              na(j) = t1;   la(j) =  t2;      ja(j) = t3
           ENDIF
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE sort_spbasis
!
!                  SUBROUTINE to fix the cutoff in mesh points
!                  for the harmonic oscillator wave functions
!                  which depend on the oscillator parameter.    
!
SUBROUTINE setup_ho_cutoff
  USE single_particle_orbits
  USE constants
  USE partial_waves
  USE wave_functions  
  IMPLICIT NONE
  INTEGER :: h,nh, lh, iq, number_of_iterations, int_points
  REAL(DP) :: sigma, norm, oscl_r
  REAL(DP) :: qmin, qmax, sum_hf, sum_norm(0:lmax)
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: q_points, weight_q
  REAL(DP) :: z_rel, factor, xp, ph, cx(0:200), contrib

  sigma = 1.0_dp; norm = 1.0_dp; int_points = 10
  number_of_iterations = 0
  nh = nmax; oscl_r=oscl*SQRT(2.)
  !  for momentum space if we set z= 0.5*(oscl_r*k)**2 > 60-70, then
  !  EXP(-60) < 10E-27. Making it twice as large ensures that we account
  !  for extensions due to the Laguerre polynoms which depend on z**(2n).
  cutoff = 2*SQRT(60.0_dp)/oscl
  qmin = 0.0D0; qmax = cutoff
  DO WHILE( (number_of_iterations < 20) .AND. (ABS(sigma) > 1E-4) )         
     ALLOCATE ( q_points(int_points), weight_q(int_points) )
     CALL gauss_legendre(qmin,qmax,q_points,weight_q,int_points)
     sum_norm  = 0.0D0
     DO h=0, lmax
        lh = h
        sum_hf = 0.
        factor = 0.5D0*((nh+lh+2)*LOG(2.D0)+fac(nh)-dfac(2*nh+2*lh+1)-0.5D0*LOG(pi))
        factor = EXP(factor)
        DO iq=1, int_points
           z_rel= q_points(iq)*q_points(iq)*oscl_r*oscl_r
           CALL laguerre_general( nh, lh+0.5D0, z_rel, cx )
           xp = EXP(-z_rel*0.5)*((q_points(iq)*oscl_r)**lh)*cx(nh)
           contrib = xp*factor*(oscl_r**(1.5D0)) 
           sum_hf=sum_hf+ weight_q(iq)*(contrib*q_points(iq))**2
        ENDDO
        sum_norm(h) = sum_hf
     ENDDO
     sigma = 0.D0
     DO h = 0, lmax
        sigma = sigma +ABS( sum_norm(h)-norm)
     ENDDO
     sigma = sigma/(lmax+1)
     number_of_iterations = number_of_iterations +1 
     !     WRITE(6,*) 'Sigma for this iteration', sigma, number_of_iterations
     DEALLOCATE ( weight_q, q_points)
     IF (ABS(sigma) > 1E-4) THEN 
        int_points = int_points+10
     ENDIF
  ENDDO
  n_rel = int_points
  n_cm_mom = int_points
  WRITE(6,'(41H New HO cutoff and # integration points: ,E12.6,2X,I5)') &
       cutoff, n_rel

END SUBROUTINE setup_ho_cutoff
!
!           Find total number of two-body mtx-elements
!           mtx elements in the rel-CoM frame
!
SUBROUTINE find_number_twobody !(number_twobody_elements, number_pp_elements, number_pn_elements, number_nn_elements )
  USE single_particle_orbits
  USE configurations
  USE constants
  USE stored_bare_interaction
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmatrix_configs 
  INTEGER ::  p_parity, ang_mom, isospin_z, i, j, ia, ib, ic, id, pq_confs
!  INTEGER, INTENT(INOUT) :: number_twobody_elements, number_pp_elements, number_pn_elements, number_nn_elements 

  number_twobody_elements = 0; number_pp_elements = 0
  number_pn_elements = 0; number_nn_elements = 0
  !     loop over isospin projection
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative         
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=j_lab_min,j_lab_max
           !     find all possible configurations, large space and model space 

           IF (type_of_renormv=='no-core')  THEN 
              CALL number_nocore_confs(ang_mom,p_parity,isospin_z,gmatrix_configs)
           ELSEIF (type_of_renormv=='v-nrg')  THEN 
              CALL number_nocore_confs(ang_mom,p_parity,isospin_z,gmatrix_configs)
           ELSE
              CALL number_gmatrix_confs(ang_mom,p_parity,isospin_z,gmatrix_configs)
           ENDIF
           pq_confs=gmatrix_configs%number_confs
           ALLOCATE(gmatrix_configs%config_ab(pq_confs+pq_confs) )
           IF (type_of_renormv=='no-core')  THEN 
              CALL setup_nocore_configurations &
                   (ang_mom,p_parity,isospin_z,gmatrix_configs)   
           ELSEIF (type_of_renormv=='v-nrg')  THEN 
              CALL setup_nocore_configurations &
                   (ang_mom,p_parity,isospin_z,gmatrix_configs)   
           ELSE
              CALL setup_gmatrix_configurations &
                   (ang_mom,p_parity,isospin_z,gmatrix_configs)   
           ENDIF
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
           DO i=1,gmatrix_configs%number_confs
              ia= gmatrix_configs%config_ab(i*2-1)
              ib= gmatrix_configs%config_ab(i*2)
              DO j=i,gmatrix_configs%number_confs
                 ic= gmatrix_configs%config_ab(j+j-1)
                 id= gmatrix_configs%config_ab(j+j)
                 IF ( isospin_z == -1) number_pp_elements = number_pp_elements +1
                 IF ( isospin_z == 0) number_pn_elements = number_pn_elements +1
                 IF ( isospin_z == 1) number_nn_elements = number_nn_elements +1
                 number_twobody_elements = number_twobody_elements+1
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE find_number_twobody

SUBROUTINE mscheme_orbits
  USE single_particle_orbits
  USE wave_functions
  USE constants
  IMPLICIT NONE
  CHARACTER(LEN= 10) ::   space, model
  INTEGER :: i, j, l, m, n, shell, isoz, count
  INTEGER :: a, b, ma, la, ja, na, ta, mb, lb, jb, nb, tb, ia, ib
  REAL(DP) :: energy, e_ho, e_hodiag

  ! First we find the number of orbits in mscheme based on the j-scheme number of
  ! orbits and the relevant quantum numbers
  count = 0
  DO i = 1, all_orbit%total_orbits
     !  note that the angular momentum is 2*j, thus no factor of 2 for 2j+1
     count = count + all_orbit%jj(i)+1
  ENDDO
  mscheme_basis%total_orbits = count
  !     Allocate space in heap for all single-particle data
  CALL allocate_sp_array(mscheme_basis,mscheme_basis%total_orbits) 
  WRITE(6,*) 'Quantum numbers n, l, j, shell, m, isopin, energy, particle/hole, status'
  count = 0
  DO i=1, all_orbit%total_orbits 
     n = all_orbit%nn(i)
     l=  all_orbit%ll(i)
     j= all_orbit%jj(i)
     shell = all_orbit%nshell(i)
     energy = all_orbit%e(i)
     isoz= all_orbit%itzp(i)   
     space = all_orbit%orbit_status(i)
     model = all_orbit%model_space(i)
     DO m = -j,j,2
        count = count + 1
        IF ( count > mscheme_basis%total_orbits+1) THEN
           WRITE(6,*) 'Serious counting error, more states than possible'
           STOP
        ENDIF
        mscheme_basis%mvalue(count) = m
        mscheme_basis%nn(count) = n 
        mscheme_basis%ll(count)=  l
        mscheme_basis%jj(count)=  j
        mscheme_basis%nshell(count)= shell
        mscheme_basis%e(count)= energy
        mscheme_basis%itzp(count)= isoz
        mscheme_basis%orbit_status(count) = space 
        mscheme_basis%model_space(count)=  model
     ENDDO
  ENDDO
! Then print the quantum numbers
  WRITE(9,'(45HLegend:         n     l     2j   2mj   2tz   )') 
  DO i=1,   mscheme_basis%total_orbits
     WRITE(9,1001) i, mscheme_basis%nn(i), mscheme_basis%ll(i), mscheme_basis%jj(i), &
          mscheme_basis%mvalue(i), mscheme_basis%itzp(i)
  ENDDO
1001 FORMAT('Orbit number:', I3,2X,  5(1X,I3))

END SUBROUTINE mscheme_orbits
!
!       This function sets up all possible two-body configs in m-scheme and
!       performs the transformation to m-scheme, the basis is given by T_z, 
!       parity and two-body M, which varies from -J to J
!
SUBROUTINE convert_to_mscheme
  USE single_particle_orbits
  USE configurations
  USE constants
  USE ang_mom_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: mscheme_configs
  INTEGER ::  p_parity, m, isospin_z, n_confs, ket, bra, i, jab_min, &
       jab_max, jcd_min, jcd_max
  INTEGER :: na, nb, nc, nd, la, lb, lc, ld, ja, jb, jc, jd, jtmin, jtmax, jtot
  INTEGER :: ma, mb, mc, md, itza, itzb, itzc, itzd, a, b, c, d
  INTEGER :: ia, ib, ic, id
  REAL(DP), ALLOCATABLE :: final_v(:,:)
  REAL(DP) :: angmom_fact, sum
  REAL(DP), DIMENSION(n_startenergy_g) :: ans
  !                          
  !     loop over isospin projection
  !
  WRITE(14,*) ' Matrix elements < ab | V | cd> in mscheme'   
  WRITE(14,*) 'Legend: a  b  c  d and  < ab | V | cd>'   
  DO isospin_z=itzmin,itzmax 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO m=-j_lab_max,j_lab_max
           !     find all possible configurations 
           CALL number_mscheme_confs &
                (m,p_parity,isospin_z,mscheme_configs)
           IF (mscheme_configs%number_confs <= 0 ) CYCLE
           n_confs=mscheme_configs%number_confs
           ALLOCATE(mscheme_configs%config_ab(n_confs+n_confs) )
           CALL setup_mscheme_configurations &
                (m,p_parity,isospin_z,mscheme_configs)
           ALLOCATE(final_v(n_confs, n_confs)); final_v = 0.
           !          here we set up the matrix elements in m-scheme
           DO bra = 1, n_confs
              !          set up the relevant quantum numbers for |ab>
              ia=mscheme_configs%config_ab(bra*2-1)
              ib=mscheme_configs%config_ab(bra*2)
              ! orbit a
              na = mscheme_basis%nn(ia); la= mscheme_basis%ll(ia)
              ja=mscheme_basis%jj(ia); ma=mscheme_basis%mvalue(ia)
              itza=mscheme_basis%itzp(ia)
              nb = mscheme_basis%nn(ib); lb=mscheme_basis%ll(ib)
              itzb=mscheme_basis%itzp(ib); mb=mscheme_basis%mvalue(ib)           
              jb=mscheme_basis%jj(ib)
              ! get identifier in jj-scheme
              b = 0; a = 0
              DO i = 1, all_orbit%total_orbits
                 IF ( (all_orbit%nn(i)==na).AND.(all_orbit%ll(i)==la).AND.  & 
                      (all_orbit%jj(i)==ja).AND.(all_orbit%itzp(i)==itza  )) a = i
                 IF ( (all_orbit%nn(i)==nb).AND.(all_orbit%ll(i)==lb).AND.  & 
                      (all_orbit%jj(i)==jb).AND.(all_orbit%itzp(i)==itzb  )) b = i
              ENDDO
              IF ( b == 0) CYCLE;     IF ( a == 0) CYCLE
              jab_min = ABS(ja-jb); jab_max = (ja+jb)
              DO ket=1,mscheme_configs%number_confs
                 ic=mscheme_configs%config_ab(ket+ket-1)
                 id=mscheme_configs%config_ab(ket+ket)
                 nc = mscheme_basis%nn(ic); lc= mscheme_basis%ll(ic)
                 jc=mscheme_basis%jj(ic); mc=mscheme_basis%mvalue(ic)
                 itzc=mscheme_basis%itzp(ic)
                 nd = mscheme_basis%nn(id); ld=mscheme_basis%ll(id)
                 itzd=mscheme_basis%itzp(id); md=mscheme_basis%mvalue(id)           
                 jd=mscheme_basis%jj(id)
                 ! get identifier in jj-scheme
                 c=0; d = 0
                 DO i = 1, all_orbit%total_orbits
                    IF ( (all_orbit%nn(i)==nd).AND.(all_orbit%ll(i)==ld).AND.  & 
                         (all_orbit%jj(i)==jd).AND.(all_orbit%itzp(i)==itzd  )) d = i
                    IF ( (all_orbit%nn(i)==nc).AND.(all_orbit%ll(i)==lc).AND.  & 
                         (all_orbit%jj(i)==jc).AND.(all_orbit%itzp(i)==itzc  )) c = i
                 ENDDO
                 IF ( d == 0) CYCLE; IF ( c == 0) CYCLE
                 jcd_min = ABS(jc-jd); jcd_max = (jc+jd)
                 jtmin = MAX( jcd_min, jab_min); jtmax = MIN(jab_max,jcd_max)
                 IF ( jtmin > jtmax) CYCLE
                 !   Now call the function which sets up the interaction 
                 !   in jj-coupling and fetches the interaction
                 sum = 0.
                 DO jtot = jtmin, jtmax, 2
                  CALL pphhmtx(a,b,c,d,jtot/2,ans)
                   angmom_fact = tjs(ja,jb,jtot,ma,mb,-2*m)* &
                         tjs(jc,jd,jtot,mc,md,-2*m)* &
                         (jtot+1)*(-1)**( (ja -jb + jc - jd) / 2 )
                    sum = sum + angmom_fact*ans(n_startenergy_g/2+1)
                 ENDDO
                 final_v(bra,ket) = sum
              ENDDO
           ENDDO
           !     print it
           CALL print_mscheme(final_v,mscheme_configs)
           !     free space
           DEALLOCATE(mscheme_configs%config_ab); DEALLOCATE(final_v)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE convert_to_mscheme
!
!                    Print out m-scheme mtx elements
!
SUBROUTINE print_mscheme(gna,mscheme_configs)
  USE configurations
  USE constants
  USE single_particle_orbits
  USE wave_functions
  IMPLICIT NONE
  TYPE (configuration_descriptor), INTENT(IN)  :: mscheme_configs
  INTEGER :: i, j, ia, ib, ic,id
  DOUBLE PRECISION, DIMENSION(mscheme_configs%number_confs, &
       mscheme_configs%number_confs), INTENT(IN)  :: gna

  DO i=1,mscheme_configs%number_confs
     ia=mscheme_configs%config_ab(i*2-1)
     ib=mscheme_configs%config_ab(i*2)
     DO j=1,mscheme_configs%number_confs
        ic=mscheme_configs%config_ab(j+j-1)
        id=mscheme_configs%config_ab(j+j)
        WRITE(14,1001) ia, ib, ic, id, gna(j,i)
     ENDDO
  ENDDO
1001 FORMAT(1X,I2,1X,I2,1X,I2,1x,I2,5X,(E12.6,1X))

END SUBROUTINE print_mscheme
!
!     This function returns the matrix for V
!
SUBROUTINE pphhmtx(ja,jb,jc,jd,jt,gmtpn)
  USE stored_bare_interaction
  USE single_particle_orbits
  USE constants
  IMPLICIT NONE
  LOGICAL TRIAG
  INTEGER, INTENT(IN) :: ja,jb,jc,jd,jt
  INTEGER :: iph
  REAL(DP), INTENT(INOUT), DIMENSION(n_startenergy_g) :: gmtpn
  REAL(DP) :: dij, delta

  gmtpn=0.
  IF(2*all_orbit%nn(ja)+all_orbit%ll(ja) +2*all_orbit%nn(jb)+all_orbit%ll(jb) > nlmax ) RETURN
  IF(2*all_orbit%nn(jc)+all_orbit%ll(jc) +2*all_orbit%nn(jd)+all_orbit%ll(jd) > nlmax ) RETURN
  IF((all_orbit%itzp(ja)+all_orbit%itzp(jb)) /= &
       (all_orbit%itzp(jc)+all_orbit%itzp(jd))) RETURN
  IF((-1)**(all_orbit%ll(ja)+all_orbit%ll(jb)) /=  &
       (-1)**(all_orbit%ll(jc)+all_orbit%ll(jd))) RETURN
  IF((ja == jb).AND.(MOD(jt,2)/=0)) RETURN       
  IF((jc == jd).AND.(MOD(jt,2)/=0)) RETURN       
  IF(triag(all_orbit%jj(ja),all_orbit%jj(jb),2*jt)) RETURN
  IF(triag(all_orbit%jj(jc),all_orbit%jj(jd),2*jt)) RETURN
  CALL mtx_elements(ja,jb,jc,jd,jt,gmtpn)

END SUBROUTINE pphhmtx
