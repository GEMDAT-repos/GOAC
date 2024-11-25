!__author__ = "Konstantin Koester"
!__copyright__ = "Copyright 2024, GOAC"
!__credits__ = ["Konstantin Koester", "Tobias Binninger", "Payam Kaghazch"]
!__license__ = "MIT"
!__version__ = "0.1.0"
!__maintainer__ = ""
!__email__ = "p.kaghazchi@fz-juelich.de"
!__status__ = "Development"

!Routine for Greedy heuristic
subroutine greedy(species_nums, species_occs, shared_sites, const_e, a_e, b_e, n, tol, & 
                  energies, solutions, max_site_num, num_species)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n
        !f2py intent(in) max_site_num, num_species
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        logical:: best_solutions(n, n, num_species, max_site_num), &
                  current_solution(num_species, max_site_num)
        integer:: i, j, l, kill
        logical:: occupied, valid_solutions, solution_unique
        real(kind=prec):: energy, current_energy, best_energies(n,n)
        logical, intent(out):: solutions(n, num_species, max_site_num)
        real(kind=prec), intent(out):: energies(n) 
        !f2py intent(out) energies, solutions

        solutions = .FALSE.
        kill = 0

        do while(.not. valid_solutions(num_species, max_site_num, n, species_nums, &
                                      species_occs, shared_sites, solutions) .and. kill .lt. 10d7)
                kill = kill +1
                best_energies = HUGE(1.0D+0)
                best_solutions = .FALSE.
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, current_solution, current_energy)
                do l=1, n 
                        do i=1, num_species
                                !Check if species occ already fullfilled
                                if(COUNT(solutions(l,i,:)) .lt. species_occs(i)) then
                                        do j=1, species_nums(i)
                                                !Check if current site is occupied by spcies or shared species
                                                if(.not. solutions(l,i,j) .and. &
                                                   .not. occupied(i, j, num_species, max_site_num, &
                                                                  solutions(l,:,:), shared_sites)) then
                                                        current_solution = solutions(l,:,:)
                                                        current_solution(i,j) = .TRUE.
                                                        current_energy = energy(num_species, max_site_num, species_nums, &
                                                                                current_solution, const_e, a_e, b_e)
                                                        if(current_energy .lt. MAXVAL(best_energies(l,:),1)) then
                                                                if(tol .gt. 0.d0 .and. & 
                                                                   MINVAL(ABS(best_energies(l,:)-current_energy),1) &
                                                                   .lt. tol) then
                                                                        cycle
                                                                end if
                                                                best_solutions(l,MAXLOC(best_energies(l,:),1),:,:) &
                                                                              = current_solution
                                                                best_energies(l,MAXLOC(best_energies(l,:),1)) &
                                                                              = current_energy
                                                        end if
                                                end if
                                        end do
                                end if
                        end do
                end do
                !$OMP END PARALLEL DO
                solutions = .FALSE.
                energies = HUGE(1.0D+0)
                do i=1, n
                        do j=1, n
                                if(best_energies(i,j) .lt. MAXVAL(energies,1)) then
                                        if(tol .gt. 0.d0 .and. &
                                           MINVAL(ABS(energies-best_energies(i,j)),1) &
                                           .lt. tol) then
                                                        cycle
                                        end if
                                        if(solution_unique(num_species, max_site_num, n, species_nums, &
                                                           best_solutions(i,j,:,:), solutions)) then
                                                solutions(MAXLOC(energies,1),:,:) = best_solutions(i,j,:,:)
                                                energies(MAXLOC(energies,1)) = best_energies(i,j)
                                        end if
                                end if
                        end do
                end do
        end do
        if(kill .ge. 10d7) then
                write(*,*) "WARNING: Subroutine killed because of too many while iterations"
        end if
        return
end subroutine greedy

subroutine random_samples(species_nums, species_occs, shared_sites, const_e, a_e, b_e, n, &
                          energies, solutions, max_site_num, num_species)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n
        !f2py intent(in) max_site_num, num_species
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species)
        real(kind=prec), intent(in):: const_e
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        integer:: i, l, kill, random_site, itarray(8), inumin
        integer,dimension(:),allocatable::iseed
        logical:: valid_solutions, occupied
        real(kind=prec):: energy, random
        logical, intent(out):: solutions(n, num_species, max_site_num)
        real(kind=prec), intent(out):: energies(n)
        !f2py intent(out) energies, solutions

        !Prepare random numbers
        call random_seed(size=inumin)
        allocate(iseed(inumin))
        call date_and_time(values=itarray)
        iseed=itarray(8)+itarray(7)*itarray(6)
        call random_seed(put=iseed)
        do i=1, 10
                call random_number(random)
        end do
        
        solutions = .FALSE.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, kill, random, random_site)
        do l=1, n
                do i=1, num_species
                        kill = 0
                        do while(COUNT(solutions(l,i,:)) .lt. species_occs(i) .and. kill .lt. 10d7)
                                kill = kill + 1
                                call random_number(random)
                                random_site = int(random*species_nums(i)+1d0)
                                if(.not. solutions(l,i,random_site) .and. & 
                                        .not. occupied(i, random_site, num_species, &
                                        max_site_num, solutions(l,:,:), shared_sites)) then
                                        solutions(l,i,random_site) = .TRUE.
                                end if
                        end do
                        if(kill .ge. 10d7) then
                                write(*,*) "Warning: Aborted while because of too many iterations!"
                        end if
                end do
                energies(l) = energy(num_species, max_site_num, species_nums, &
                                     solutions(l,:,:), const_e, a_e, b_e)
                
                if(.not. valid_solutions(num_species, max_site_num, 1, species_nums, &
                         species_occs, shared_sites, solutions(l,:,:))) then
                        write(*,*) "Warning: Generated non valid solutions!"
                end if
        end do
        !$OMP END PARALLEL DO
        deallocate(iseed)
        return
end subroutine random_samples

subroutine local_minimizer(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, &
                           stop_time, stop_no_improve_steps, &
                           max_site_num, num_species, n_solutions, lowest_energies, lowest_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n_solutions, n, stop_time, stop_no_improve_steps
        !f2py intent(in) max_site_num, num_species, n_solutions
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species), &
                              solutions(n_solutions, num_species, max_site_num)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        logical, allocatable:: opt_solutions(:,:,:,:)
        real(kind=prec), allocatable:: opt_energies(:,:)
        real(kind=prec), intent(out):: lowest_energies(n_solutions, n)
        logical, intent(out):: lowest_solutions(n_solutions, n, num_species, max_site_num)
        !f2py intent(out) lowest_energies, lowest_solutions
        logical:: energy_changed, solution_unique
        integer:: i,j,k,kill, steps_no_improve
        real(kind=prec):: energy, start, now, start_time, prev_min_energy

        call cpu_time(start)
        start_time = OMP_GET_WTIME()
        allocate(opt_solutions(n_solutions, n, num_species, max_site_num))
        allocate(opt_energies(n_solutions, n))
        lowest_energies = HUGE(1.0D+0)
        lowest_solutions = .FALSE.
        call branch_n_bound(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, &
                            max_site_num, num_species, n_solutions, opt_energies, opt_solutions)
        do i=1, n_solutions
                lowest_solutions(i,1,:,:) = solutions(i,:,:) 
                lowest_energies(i,1) = energy(num_species, max_site_num, species_nums, &
                                              solutions(i,:,:), const_e, a_e, b_e)
                do j=1, n
                        if(opt_energies(i,j) .lt. MAXVAL(lowest_energies(i,:),1)) then
                               lowest_solutions(i,MAXLOC(lowest_energies(i,:),1),:,:) = opt_solutions(i,j,:,:)
                               lowest_energies(i,MAXLOC(lowest_energies(i,:),1)) = opt_energies(i,j)
                        end if
                end do 
        end do
        deallocate(opt_solutions)
        deallocate(opt_energies)

        allocate(opt_solutions(n, n, num_species, max_site_num))
        allocate(opt_energies(n, n))
        steps_no_improve = 0
        prev_min_energy = HUGE(1.0D+0)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(energy_changed, kill, opt_energies, opt_solutions, j, k, now)
        do i=1, n_solutions
                if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                        cycle
                elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                        cycle
                end if
                energy_changed = .TRUE.
                kill = 0
                do while(energy_changed .and. kill .lt. 10d7)
                        kill = kill + 1
                        energy_changed = .FALSE.
                        call branch_n_bound(species_nums, species_occs, shared_sites, lowest_solutions(i,:,:,:), &
                                            const_e, a_e, b_e, n, tol, max_site_num, num_species, n, &
                                            opt_energies, opt_solutions)
                        do j=1, n
                                do k=1, n
                                        if(opt_energies(j,k) .lt. MAXVAL(lowest_energies(i,:),1)) then
                                                if(tol .gt. 0.d0 .and. &
                                                   MINVAL(ABS(lowest_energies(i,:)-opt_energies(j,k)),1) &
                                                   .lt. tol) then
                                                        cycle
                                                end if
                                                if(solution_unique(num_species, max_site_num, n, species_nums, &
                                                                   opt_solutions(j,k,:,:), lowest_solutions(i,:,:,:))) then
                                                        lowest_solutions(i,MAXLOC(lowest_energies(i,:),1),:,:) = & 
                                                                         opt_solutions(j,k,:,:)
                                                        lowest_energies(i,MAXLOC(lowest_energies(i,:),1)) = opt_energies(j,k)
                                                        energy_changed = .TRUE.
                                                end if
                                        end if
                                end do
                        end do
                end do
                if(kill .ge. 10d7) then
                        write(*,*) "WARNING: Subroutine killed because of too many while iteations"
                end if
                if((tol .gt. 0.d0 .and. ABS(MINVAL(lowest_energies)-prev_min_energy) .lt. tol) .or. & 
                   (tol .lt. 0.d0 .and. MINVAL(lowest_energies) .ge. prev_min_energy)) then
                        steps_no_improve = steps_no_improve+1
                else
                        steps_no_improve = 0
                        prev_min_energy = MINVAL(lowest_energies)
                end if
                call cpu_time(now)
                write(*,*) "CPU-TIME:", (now-start), "GLOBAL-MINIMUM:", MINVAL(lowest_energies)
        end do
        !$OMP END PARALLEL DO
        deallocate(opt_solutions)
        deallocate(opt_energies)
        return
end subroutine local_minimizer

subroutine branch_n_bound(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, & 
                          max_site_num, num_species, n_solutions, opt_energies, opt_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n_solutions, n
        !f2py intent(in) max_site_num, num_species, n_solutions
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species), &
                              solutions(n_solutions, num_species, max_site_num)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        logical, intent(out):: opt_solutions(n_solutions, n, num_species, max_site_num)
        real(kind=prec), intent(out):: opt_energies(n_solutions, n)
        !f2py intent(out) opt_energies, opt_solutions
        integer:: i,j,k,l,m
        logical:: current_solution(num_species, max_site_num), valid_solutions
        real(kind=prec):: current_energy, energy

        opt_energies = HUGE(1.0D+0)
        opt_solutions = .FALSE.
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k, m, current_solution, current_energy)
        do l=1, n_solutions
                do i=1, num_species
                        do j=1, species_nums(i)
                                do k=j+1, species_nums(i)
                                        if(solutions(l,i,k) .neqv. solutions(l,i,j)) then
                                                current_solution = solutions(l,:,:)
                                                current_solution(i,j) = solutions(l,i,k)
                                                current_solution(i,k) = solutions(l,i,j)
                                                do m=1, num_species
                                                        if(i .ne. m .and. shared_sites(i,m)) then
                                                                if(current_solution(i,j) .and. current_solution(m,j)) then
                                                                        current_solution(m,j) = solutions(l,m,k)
                                                                        current_solution(m,k) = solutions(l,m,j)
                                                                        exit
                                                                elseif (current_solution(i,k) .and. current_solution(m,k)) then
                                                                        current_solution(m,j) = solutions(l,m,k)
                                                                        current_solution(m,k) = solutions(l,m,j)
                                                                        exit
                                                                end if
                                                        end if
                                                end do
                                                if(.not. valid_solutions(num_species, max_site_num, 1, species_nums, &
                                                                         species_occs, shared_sites, current_solution)) then
                                                        write(*,*) "ERROR: Invald Solution was considered!"
                                                end if
                                                current_energy = energy(num_species, max_site_num, species_nums, &
                                                                        current_solution, const_e, a_e, b_e)
                                                if(current_energy .lt. MAXVAL(opt_energies(l,:),1)) then
                                                        if(tol .gt. 0.d0 .and. &
                                                           MINVAL(ABS(opt_energies(l,:)-current_energy),1) &
                                                           .lt. tol) then
                                                                cycle
                                                        end if
                                                        opt_solutions(l,MAXLOC(opt_energies(l,:),1),:,:) &
                                                                = current_solution
                                                        opt_energies(l,MAXLOC(opt_energies(l,:),1)) &
                                                                = current_energy
                                                end if
                                        end if
                                end do
                        end do
                end do
        end do
        !$OMP END PARALLEL DO
        return
end subroutine branch_n_bound

!Check if current solution is unique in best solutions
function solution_unique(num_species, max_site_num, n, species_nums, current_solution, best_solutions)
        implicit none
        integer, intent(in):: n, num_species, max_site_num 
        integer, intent(in):: species_nums(num_species)
        logical, intent(in):: current_solution(num_species, max_site_num), best_solutions(n, num_species, max_site_num)
        logical:: solution_unique
        integer:: i,j,k

        do k=1, n
                solution_unique = .FALSE.
                do i=1, num_species
                        do j=1, species_nums(i)
                                if(current_solution(i,j) .neqv. best_solutions(k,i,j)) then
                                        solution_unique = .TRUE.
                                        exit
                                end if
                        end do
                        if(solution_unique) then
                                exit
                        end if 
                end do
                if(.not. solution_unique) then
                        return
                end if
        end do
        return
end function solution_unique

!Check if all solutions are valid
function valid_solutions(num_species, max_site_num, n, species_nums, species_occs, shared_sites, solutions)
        implicit none
        integer, intent(in):: num_species, max_site_num, n
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species)
        logical, intent(in):: solutions(n, num_species, max_site_num)
        logical:: valid_solutions, occupied
        integer:: i, j, l  

        valid_solutions = .TRUE.

        do l=1, n
                do i=1, num_species
                        if(COUNT(solutions(l, i,:)) .ne. species_occs(i)) then
                                valid_solutions = .FALSE.
                                return
                        end if
                        do j=1, species_nums(i)
                                if(solutions(l,i,j) .and. occupied(i, j, num_species, &
                                                                   max_site_num, &
                                                                   solutions(l,:,:), shared_sites)) then
                                        valid_solutions = .FALSE.
                                        return
                                end if
                        end do
                end do
        end do
        return
end function valid_solutions

!Calculate energy for a given solution
function energy(num_species, max_site_num, species_nums, solution, const_e, a_e, b_e)
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: num_species, max_site_num
        integer, intent(in):: species_nums(num_species)
        logical, intent(in):: solution(num_species, max_site_num)
        real(kind=prec), intent(in):: const_e
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        real(kind=prec):: energy
        integer:: i,j,k,l, start

        energy = const_e
        do i=1, num_species
                do j=1, species_nums(i)
                        if(solution(i,j)) then
                                energy = energy + a_e(i,j)
                                do k=i, num_species
                                        if(k .eq. i) then
                                                start = j+1
                                        else
                                                start = 1
                                        end if
                                        do l=start, species_nums(k)
                                                if(solution(k,l)) then
                                                        energy = energy + b_e(i,j,k,l)
                                                end if
                                        end do
                                end do
                        end if
                end do
        end do
        return
end function energy

!Check for occupations of a given site+position by shared species
function occupied(i_in, j_in, num_species, max_site_num, solution, shared_sites)
        implicit none
        integer, intent(in):: i_in, j_in, num_species, max_site_num
        logical, intent(in):: shared_sites(num_species, num_species)
        logical, intent(in):: solution(num_species, max_site_num)
        logical:: occupied
        integer:: i

        do i=1, num_species
                if(i .ne. i_in .and. shared_sites(i_in,i)) then
                        if(solution(i,j_in)) then
                                occupied = .TRUE.
                                return
                        end if
                end if
        end do

        occupied = .FALSE.
        return
end function occupied                

!Replica Excahnge Monte Carlo run subroutine
subroutine remc(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, t_steps, t_repeat, kT, &
                stop_time, stop_no_improve_steps, &
                max_site_num, num_species, n_solutions, n_kT, final_energies, final_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n_solutions, n, t_steps, n_kT, t_repeat, &
                              stop_time, stop_no_improve_steps
        !f2py intent(in) max_site_num, num_species, n_solutions, n_kT
        !f2py sim_an = 1d0, sim_an_steps = 1000000000, write_steps = 100000
        real(kind=prec), intent(in):: kT(n_kT)
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species), &
                              solutions(n_solutions, num_species, max_site_num)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        real(kind=prec):: best_energies(n_solutions, n)
        logical:: best_solutions(n_solutions, n, num_species, max_site_num)
        real(kind=prec), intent(out):: final_energies(n_solutions, n)
        logical, intent(out):: final_solutions(n_solutions, n, num_species, max_site_num)
        !f2py intent(out) final_energies, final_solutions
        logical:: current_solutions(n_solutions, num_species, max_site_num), &
                  solutions_in(n_kT, n_solutions, num_species, max_site_num), &
                  flip, buffer_solution(num_species, max_site_num), &
                  final_solutions_tmp(n_kT, n_solutions, n, num_species, max_site_num)
        integer:: k, t, i, s, itarray(8), inumin, steps_no_improve
        integer,dimension(:),allocatable::iseed
        real(kind=prec):: random, energy, current_energy, compare_energy, delta, prev_min_energy, &
                          final_energies_tmp(n_kT, n_solutions, n), start, now, start_time

        !Prepare random numbers
        call random_seed(size=inumin)
        allocate(iseed(inumin))
        call date_and_time(values=itarray)
        iseed=itarray(8)+itarray(7)*itarray(6)
        call random_seed(put=iseed)
        do k=1, 10
                call random_number(random)
        end do

        call cpu_time(start)
        final_solutions = .FALSE.
        final_energies = HUGE(1.0d+0)
        final_solutions_tmp = .FALSE.
        final_energies_tmp = HUGE(1.0d+0)
        do k=1, n_kT
                solutions_in(k,:,:,:) = solutions(:,:,:)
        end do
        start_time = OMP_GET_WTIME()
        steps_no_improve = 0
        prev_min_energy = HUGE(1.0D+0)
        do t=1, t_repeat
                if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                        cycle
                elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                        cycle
                end if
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(s, i, now, best_energies, &
                !$OMP best_solutions, current_solutions, current_energy)
                do k=1, n_kT
                        if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                                cycle
                        elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                                cycle
                        end if
                        call monte_carlo(species_nums, species_occs, shared_sites, solutions_in(k,:,:,:), const_e, a_e, &
                                         b_e, n, tol, t_steps, kT(k), 1.d0, 10000000, t_steps+1, -1, -1, max_site_num, &
                                         num_species, n_solutions, best_energies, best_solutions, current_solutions)
                        solutions_in(k,:,:,:) = current_solutions(:,:,:)
                        do s=1, n_solutions
                                do i=1, n
                                        current_energy = best_energies(s,i)
                                        if(current_energy .lt. MAXVAL(final_energies_tmp(k,s,:),1)) then
                                                if(tol .gt. 0.d0 .and. &
                                                   MINVAL(ABS(final_energies_tmp(k,s,:)-current_energy),1) &
                                                   .lt. tol) then
                                                        cycle
                                                end if
                                                final_solutions_tmp(k,s,MAXLOC(final_energies_tmp(k,s,:),1),:,:) &
                                                              = best_solutions(s,i,:,:)
                                                final_energies_tmp(k,s,MAXLOC(final_energies_tmp(k,s,:),1)) &
                                                              = current_energy
                                        end if
                                end do
                                call cpu_time(now)
                                write(*,*) "Thread:", s, "Step:", t*t_steps, "kT:", kT(k), "Current-Energy:", &
                                           current_energy, "Lowest-Energy:", MINVAL(final_energies_tmp(k,s,:),1), & 
                                           "CPU-TIME:", (now-start)
                        end do
                end do
                !$OMP END PARALLEL DO
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k, i, current_energy, compare_energy, flip, delta, random, &
                !$OMP buffer_solution)
                do s=1, n_solutions
                        if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                                cycle
                        elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                                cycle
                        end if
                        do k=1, n_kT
                                do i=1, n
                                        current_energy = final_energies_tmp(k,s,i)
                                        if(current_energy .lt. MAXVAL(final_energies(s,:),1)) then
                                                if(tol .gt. 0.d0 .and. &
                                                   MINVAL(ABS(final_energies(s,:)-current_energy),1) &
                                                   .lt. tol) then
                                                        cycle
                                                end if
                                                final_solutions(s,MAXLOC(final_energies(s,:),1),:,:) &
                                                              = final_solutions_tmp(k,s,i,:,:)
                                                final_energies(s,MAXLOC(final_energies(s,:),1)) &
                                                              = current_energy
                                        end if
                                end do
                                current_energy = energy(num_species, max_site_num, species_nums, &
                                                        solutions_in(k,s,:,:), const_e, a_e, b_e)
                                do i=k+1, n_kT
                                        compare_energy = energy(num_species, max_site_num, species_nums, &
                                                                solutions_in(i,s,:,:), const_e, a_e, b_e)
                                        flip = .FALSE.
                                        delta = (1.d0/kT(i) - 1.d0/kT(k))*(current_energy-compare_energy)
                                        call random_number(random)
                                        if(delta .le. 0d0) then
                                                flip = .TRUE.
                                        elseif(random .lt. dexp(-delta)) then
                                                flip = .TRUE.
                                        end if
                                        if(flip) then
                                                write(*,*) "Exchange: Thread:", s, "kT1:", kT(k), "kT2:", kT(i)
                                                buffer_solution = solutions_in(k,s,:,:)
                                                solutions_in(k,s,:,:) = solutions_in(i,s,:,:)
                                                solutions_in(i,s,:,:) = buffer_solution
                                        end if
                                end do
                        end do
                end do
                !$OMP END PARALLEL DO
                if((tol .gt. 0.d0 .and. ABS(MINVAL(final_energies)-prev_min_energy) .lt. tol) .or. &
                   (tol .lt. 0.d0 .and. MINVAL(final_energies) .ge. prev_min_energy)) then
                        steps_no_improve = steps_no_improve+1
                else
                        steps_no_improve = 0
                        prev_min_energy = MINVAL(final_energies)
                end if
                call cpu_time(now)
                write(*,*) "CPU-TIME:", (now-start), "GLOBAL-MINIMUM:", MINVAL(final_energies(:,:))
        end do
        return
end subroutine remc

!Monte Carlo run subroutine
subroutine monte_carlo(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, &
                       t_steps, kT, sim_an, sim_an_steps, write_steps, stop_time, stop_no_improve_steps, &
                       max_site_num, num_species, n_solutions, best_energies, best_solutions, current_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        integer, intent(in):: max_site_num, num_species, n_solutions, n, t_steps, sim_an_steps, write_steps, &
                              stop_time, stop_no_improve_steps
        real(kind=prec), intent(in):: kT, sim_an
        !f2py intent(in) max_site_num, num_species, n_solutions
        !f2py sim_an = 1d0, sim_an_steps = 1000000000, write_steps = 100000
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species), &
                              solutions(n_solutions, num_species, max_site_num)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        real(kind=prec), intent(out):: best_energies(n_solutions, n)
        logical, intent(out):: best_solutions(n_solutions, n, num_species, max_site_num), &
                               current_solutions(n_solutions, num_species, max_site_num)
        !f2py intent(out) best_energies, best_solutions, current_solutions
        integer:: s, t, m, random_site, random_pos1, random_pos2, itarray(8), inumin, steps_no_improve
        integer,dimension(:),allocatable::iseed
        real(kind=prec):: energy, current_energy, new_energy, random, current_kT, start, now, &
                          start_time, prev_min_energy
        logical:: current_solution(num_species, max_site_num), new_solution(num_species, max_site_num), &
                  flip, valid_solutions

        !Prepare random numbers
        call random_seed(size=inumin)
        allocate(iseed(inumin))
        call date_and_time(values=itarray)
        iseed=itarray(8)+itarray(7)*itarray(6)
        call random_seed(put=iseed)
        do s=1, 10
                call random_number(random)
        end do


        call cpu_time(start)
        start_time = OMP_GET_WTIME()
        best_solutions = .FALSE.
        best_energies = HUGE(1.0d+0)
        steps_no_improve = 0
        prev_min_energy = HUGE(1.0D+0)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(current_solution, new_solution, current_energy, new_energy, current_kT, &
        !$OMP t, m, random, random_site, random_pos1, random_pos2, flip, now)
        do s=1, n_solutions
                best_solutions(s,1,:,:) = solutions(s,:,:)
                best_energies(s,1) = energy(num_species, max_site_num, species_nums, &
                                            solutions(s,:,:), const_e, a_e, b_e)
                current_solution = best_solutions(s,1,:,:)
                current_energy = best_energies(s,1)
                current_kT = kT
                if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                        cycle
                elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                        cycle
                end if
                do t=1, t_steps
                        if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                                cycle
                        elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                                cycle
                        end if
                        call random_number(random)
                        random_site = int(random*num_species+1d0)
                        call random_number(random)
                        random_pos1 = int(random*species_nums(random_site)+1.d0)
                        call random_number(random)
                        random_pos2 = int(random*species_nums(random_site)+1.d0)
                        if(current_solution(random_site,random_pos1) .neqv. current_solution(random_site,random_pos2)) then
                                !Flip can happen
                                new_solution = current_solution
                                new_solution(random_site,random_pos1) = current_solution(random_site,random_pos2)
                                new_solution(random_site,random_pos2) = current_solution(random_site,random_pos1)
                                do m=1, num_species
                                        !make double-flip if shared site
                                        if(random_site .ne. m .and. shared_sites(random_site,m)) then
                                                if(new_solution(random_site,random_pos1) .and. &
                                                   new_solution(m,random_pos1)) then
                                                        new_solution(m,random_pos1) = current_solution(m,random_pos2)
                                                        new_solution(m,random_pos2) = current_solution(m,random_pos1)
                                                        exit
                                                elseif(new_solution(random_site,random_pos2) .and. &
                                                       new_solution(m,random_pos2)) then
                                                        new_solution(m,random_pos1) = current_solution(m,random_pos2)
                                                        new_solution(m,random_pos2) = current_solution(m,random_pos1)
                                                        exit
                                                end if
                                        end if
                                end do
                                if(.not. valid_solutions(num_species, max_site_num, 1, species_nums, &
                                                         species_occs, shared_sites, new_solution)) then
                                        write(*,*) "ERROR: Invald Solution was considered!"
                                end if
                                new_energy = energy(num_species, max_site_num, species_nums, &
                                                    new_solution, const_e, a_e, b_e)
                                !Save solution if it is in n best soultions found
                                if(new_energy .lt. MAXVAL(best_energies(s,:),1)) then
                                if(.not.(tol .gt. 0.d0 .and. &
                                   MINVAL(ABS(best_energies(s,:)-new_energy),1) &
                                   .lt. tol)) then
                                        best_solutions(s,MAXLOC(best_energies(s,:),1),:,:) &
                                                     = new_solution
                                        best_energies(s,MAXLOC(best_energies(s,:),1)) &
                                                     = new_energy
                                end if
                                end if
                                !Decide if flip is accepted
                                flip = .FALSE.
                                if(new_energy .lt. current_energy) then
                                        flip = .TRUE.
                                else
                                        call random_number(random)
                                        if(random .lt. EXP(-(new_energy-current_energy)/current_kT)) then
                                                flip = .TRUE.
                                        end if
                                end if
                                if(flip) then
                                        current_solution = new_solution
                                        current_energy = new_energy
                                end if
                        end if
                        if(MODULO(t, write_steps) .eq. 0) then
                                call cpu_time(now)
                                write(*,*) "Thread:", s, "Step:", t, "kT:", current_kT, "Current-Energy:", &
                                           current_energy, "Lowest-Energy:", MINVAL(best_energies(s,:),1), &
                                           "CPU-TIME:", (now-start)
                                write(*,*) "CPU-TIME:", (now-start), "GLOBAL-MINIMUM:", MINVAL(best_energies)
                        end if
                        !Perform simulated annealing if sim_an < 1
                        if(MODULO(t, sim_an_steps) .eq. 0) then
                                current_kT = current_kT*sim_an
                        end if
                        if((tol .gt. 0.d0 .and. ABS(MINVAL(best_energies)-prev_min_energy) .lt. tol) .or. & 
                           (tol .lt. 0.d0 .and. MINVAL(best_energies) .ge. prev_min_energy)) then
                                steps_no_improve = steps_no_improve+1
                        else
                                steps_no_improve = 0
                                prev_min_energy = MINVAL(best_energies)
                        end if
                end do
                current_solutions(s,:,:) = current_solution
        end do
        !$OMP END PARALLEL DO   
        return
end subroutine monte_carlo

!Generation Exchange Genetic Algorithm
subroutine rega(species_nums, species_occs, shared_sites, const_e, a_e, b_e, n, tol, num_ga, ga_steps, &
                generation_size, pool_size, elite_size, mutation_rate, stop_time, stop_no_improve_steps, &
                write_steps, max_site_num, num_species, best_energies, best_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        real(kind=prec):: mutation_rate
        integer, intent(in):: max_site_num, num_species, n, ga_steps, write_steps, num_ga, &
                              generation_size, pool_size, elite_size, stop_time, stop_no_improve_steps
        !f2py intent(in) max_site_num, num_species
        !f2py write_steps = 100000
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        real(kind=prec), intent(out):: best_energies(n)
        logical, intent(out):: best_solutions(n, num_species, max_site_num)
        !f2py intent(out) best_energies, best_solutions
        integer:: i, j, t, itarray(8), inumin, steps_no_improve
        integer,dimension(:),allocatable::iseed
        real(kind=prec):: random, start, now, &
                          generation_energies(n, generation_size), &
                          start_time, prev_min_energy
        logical:: generation(n, generation_size, num_species, max_site_num)
        logical, allocatable:: buffer(:,:,:)

        !Prepare random numbers
        call random_seed(size=inumin)
        allocate(iseed(inumin))
        call date_and_time(values=itarray)
        iseed=itarray(8)+itarray(7)*itarray(6)
        call random_seed(put=iseed)
        do i=1, 10
                call random_number(random)
        end do

        allocate(buffer(int(generation_size/2), num_species, max_site_num))

        call cpu_time(start)
        start_time = OMP_GET_WTIME()

        do i=1, n
                !Build start generation from input solutions and random solutions
                generation(i,:,:,:) = .FALSE.
                generation_energies(i,:) = HUGE(1.0d+0)
                call random_samples(species_nums, species_occs, shared_sites, const_e, a_e, b_e, generation_size, &
                                    generation_energies(i,:), generation(i,:,:,:), max_site_num, num_species)
        end do

        best_energies = HUGE(1.0d+0)
        best_solutions = .FALSE.
        steps_no_improve = 0
        prev_min_energy = HUGE(1.0D+0)
        do t=1, num_ga
                if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                        cycle
                elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                        cycle
                end if

                !$OMP PARALLEL DO DEFAULT(SHARED)
                do i=1, n
                        call ga(species_nums, species_occs, shared_sites, generation(i,:,:,:), const_e, a_e, b_e, &
                                generation_size, tol, ga_steps, generation_size, pool_size, elite_size, mutation_rate, &
                                stop_time, -1, ga_steps+1, max_site_num, num_species, generation_size, &
                                generation_energies(i,:), generation(i,:,:,:))
                 end do
                 !$OMP END PARALLEL DO

                 !Generation Cross-Over
                 do i=1, n
                        do j=1, n
                                generation(j,int(generation_size/n)*(i-1)+1:int(generation_size/n)*i,:,:) = &
                                generation(i,int(generation_size/n)*(i-1)+1:int(generation_size/n)*i,:,:)
                        end do
                        do j=1, generation_size
                                if(generation_energies(i,j) .lt. MAXVAL(best_energies,1)) then
                                        if(tol .gt. 0.d0 .and. &
                                           MINVAL(ABS(best_energies-generation_energies(i,j)),1) &
                                           .lt. tol) then
                                                cycle
                                        end if
                                        best_solutions(MAXLOC(best_energies, 1),:,:) = generation(i,j,:,:)
                                        best_energies(MAXLOC(best_energies, 1)) = generation_energies(i,j)
                                end if
                        end do
                 end do

                if(MODULO(t, write_steps) .eq. 0) then
                        call cpu_time(now)
                        write(*,*) "STEP:", t, "GLOBAL MINIMUM:", MINVAL(best_energies), "CPU-TIME:", (now-start)
                end if
                if((tol .gt. 0.d0 .and. ABS(MINVAL(best_energies)-prev_min_energy) .lt. tol) .or. & 
                   (tol .lt. 0.d0 .and. MINVAL(best_energies) .ge. prev_min_energy)) then
                        steps_no_improve = steps_no_improve+1
                else
                        steps_no_improve = 0
                        prev_min_energy = MINVAL(best_energies)
                end if
        end do

        deallocate(buffer)
end subroutine rega


!Genetic Algorithm run subroutine
subroutine ga(species_nums, species_occs, shared_sites, solutions, const_e, a_e, b_e, n, tol, ga_steps, &
              generation_size, pool_size, elite_size, mutation_rate, stop_time, stop_no_improve_steps, & 
              write_steps, max_site_num, num_species, n_solutions, best_energies, best_solutions)
        !$ use omp_lib
        implicit none
        integer, parameter:: prec=8
        real(kind=prec):: mutation_rate
        integer, intent(in):: max_site_num, num_species, n_solutions, n, ga_steps, write_steps, &
                              generation_size, pool_size, elite_size, stop_time, stop_no_improve_steps
        !f2py intent(in) max_site_num, num_species, n_solutions
        !f2py write_steps = 100000
        integer, intent(in):: species_nums(num_species), species_occs(num_species)
        logical, intent(in):: shared_sites(num_species, num_species), &
                              solutions(n_solutions, num_species, max_site_num)
        real(kind=prec), intent(in):: const_e, tol
        real(kind=prec), intent(in):: a_e(num_species, max_site_num), &
                                      b_e(num_species, max_site_num, num_species, max_site_num)
        real(kind=prec), intent(out):: best_energies(n)
        logical, intent(out):: best_solutions(n, num_species, max_site_num)
        !f2py intent(out) best_energies, best_solutions
        integer:: i, j, k, l, p, t, itarray(8), inumin, parent_id1, parent_id2, &
                  diff, diff_start, steps_no_improve, num_swaps, swaps(max_site_num)
        integer,dimension(:),allocatable::iseed
        real(kind=prec):: energy, random, current_energy, random2, start, now, &
                          generation_energies(generation_size), pool_energies(pool_size), &
                          generation_energies_new(generation_size), cum_sum(generation_size), &
                          start_time, prev_min_energy, accept
        logical:: valid_solutions, buffer, acc_counted, &
                  generation(generation_size, num_species, max_site_num), pool(pool_size, num_species, max_site_num), &
                  generation_new(generation_size, num_species, max_site_num), solution_unique

        !Prepare random numbers
        call random_seed(size=inumin)
        allocate(iseed(inumin))
        call date_and_time(values=itarray)
        iseed=itarray(8)+itarray(7)*itarray(6)
        call random_seed(put=iseed)
        do i=1, 10
                call random_number(random)
        end do

        call cpu_time(start)
        start_time = OMP_GET_WTIME()

        !Build start generation from input solutions and random solutions
        generation = .FALSE.
        generation_energies = HUGE(1.0d+0)
        call random_samples(species_nums, species_occs, shared_sites, const_e, a_e, b_e, generation_size, &
                            generation_energies_new, generation_new, max_site_num, num_species)
        do i=1, generation_size
                if(i .le. n_solutions) then
                        generation(i,:,:) = solutions(i,:,:)
                        current_energy = energy(num_species, max_site_num, species_nums, &
                                                solutions(i,:,:), const_e, a_e, b_e)
                        generation_energies(i) = current_energy
                else
                        generation(i,:,:) = generation_new(i,:,:)
                        generation_energies(i) = generation_energies_new(i)
                end if
        end do

        !GA loop
        best_energies = HUGE(1.0d+0)
        best_solutions = .FALSE.
        steps_no_improve = 0
        prev_min_energy = HUGE(1.0D+0)
        do t=1, ga_steps+1
                if(stop_time .gt. 0 .and. OMP_GET_WTIME()-start_time .gt. stop_time) then
                        cycle
                elseif(stop_no_improve_steps .gt. 0 .and. steps_no_improve .gt. stop_no_improve_steps) then
                        cycle
                end if
                pool_energies = HUGE(1.0d+0)
                pool = .FALSE.
                current_energy = SUM(generation_energies(:))
                do i=1, generation_size
                        cum_sum(i) = SUM(generation_energies(1:i))/current_energy
                        if(generation_energies(i) .lt. MAXVAL(best_energies, 1)) then
                                if(tol .gt. 0.d0 .and. &
                                   MINVAL(ABS(best_energies-generation_energies(i)),1) &
                                   .lt. tol) then
                                        cycle
                                end if
                                best_solutions(MAXLOC(best_energies, 1),:,:) = generation(i,:,:)
                                best_energies(MAXLOC(best_energies, 1)) = generation_energies(i)
                        end if
                end do
                if(t .gt. ga_steps) then
                        exit
                end if
                !Pool creation
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k, l, parent_id1, parent_id2, random, &
                !$OMP random2, buffer, num_swaps, swaps, diff, diff_start, accept, acc_counted)
                do p=1, pool_size
                        !Parent selection by weighted roulete wheel
                        call random_number(random)
                        call random_number(random2)
                        parent_id1 = -1
                        parent_id2 = -1
                        do i=1, generation_size
                                if(random .lt. cum_sum(i) .and. parent_id1 .lt. 0) then
                                        parent_id1 = i
                                elseif(random2 .lt. cum_sum(i) .and. parent_id2 .lt. 0) then
                                        parent_id2 = i
                                end if
                                if(parent_id1 .gt. 0 .and. parent_id2 .gt. 0) then
                                        exit
                                end if
                        end do
                        if(parent_id1 .eq. parent_id2) then
                                write(*,*) "Parent IDs equal"
                        end if
                        
                        pool(p,:,:) = generation(parent_id1,:,:)
                        !Offspring generation
                        !diff = 0
                        !do i=1, num_species
                        !        do j=1, species_nums(i)
                        !                if(pool(p,i,j) .and. .not. generation(parent_id2,i,j)) then
                        !                        diff = diff + 1
                        !                end if
                        !        end do
                        !end do
                        do i=1, num_species
                                do j=1, species_nums(i)
                                        if(pool(p,i,j) .and. .not. generation(parent_id2,i,j)) then
                                                call random_number(random)
                                                if(random .lt. 1.d0/2.d0) then
                                                        accept = 0
                                                        num_swaps = 0
                                                        swaps = 0
                                                        do k=1, species_nums(i)
                                                                if(.not. pool(p,i,k) .and. generation(parent_id2,i,k)) then
                                                                        num_swaps = num_swaps+1
                                                                        swaps(num_swaps) = k
                                                                        acc_counted = .FALSE.
                                                                        do l=1, num_species
                                                                                if(i .ne. l .and. shared_sites(i,l)) then
                                                                                        if(pool(p,l,k) .and. &
                                                                                           generation(parent_id2,l,j)) then
                                                                                                accept = accept + 1.d0/2.d0
                                                                                                acc_counted = .TRUE.
                                                                                                exit
                                                                                        end if
                                                                                end if
                                                                        end do
                                                                        if(.not. acc_counted) then
                                                                                accept = accept + 1
                                                                        end if
                                                                end if
                                                        end do
                                                        call random_number(random)
                                                        if(num_swaps .lt. 1 .or. random .ge. accept/float(num_swaps)) then
                                                                cycle
                                                        end if
                                                        call random_number(random)
                                                        k = swaps(int(random*num_swaps+1.d0))
                                                        buffer = pool(p,i,j)
                                                        pool(p,i,j) = pool(p,i,k)
                                                        pool(p,i,k) = buffer
                                                        do l=1, num_species
                                                                !make double-flip if shared site
                                                                if(i .ne. l .and. shared_sites(l,i)) then
                                                                        if(pool(p,i,j) .and. pool(p,l,j)) then
                                                                                buffer = pool(p,l,j)
                                                                                pool(p,l,j) = pool(p,l,k)
                                                                                pool(p,l,k) = buffer
                                                                                exit
                                                                        elseif(pool(p,i,k) .and. pool(p,l,k)) then
                                                                                buffer = pool(p,l,j)
                                                                                pool(p,l,j) = pool(p,l,k)
                                                                                pool(p,l,k) = buffer
                                                                                exit
                                                                        end if
                                                                end if
                                                        end do
                                                end if
                                        end if
                                end do
                        end do

                        !diff_start = diff
                        !do i=1, num_swaps
                        !        call random_number(random)
                        !        if(random .lt. 0.125d0) then
                        !                buffer = pool(p,swaps(i,1), swaps(i,2))
                        !                pool(p,swaps(i,1), swaps(i,2)) = &
                        !                pool(p,swaps(i,1), swaps(i,3))
                        !                pool(p,swaps(i,1), swaps(i,3)) = buffer
                        !                do l=1, num_species
                        !                        !make double-flip if shared site
                        !                        if(swaps(i,1) .ne. l .and. shared_sites(swaps(i,1),l)) then
                        !                                if(pool(p,swaps(i,1), swaps(i,2)) .and. &
                        !                                   pool(p,l,swaps(i,2))) then
                        !                                        buffer = pool(p,l,swaps(i,2))
                        !                                        pool(p,l,swaps(i,2)) = pool(p,l,swaps(i,3))
                        !                                        pool(p,l,swaps(i,3)) = buffer
                        !                                        exit
                        !                                elseif(pool(p,swaps(i,1),swaps(i,3)) .and. &
                        !                                       pool(p,l,swaps(i,3))) then
                        !                                        buffer = pool(p,l,swaps(i,2))
                        !                                        pool(p,l,swaps(i,2)) = pool(p,l,swaps(i,3))
                        !                                        pool(p,l,swaps(i,3)) = buffer
                        !                                        exit
                        !                                end if
                        !                        end if
                        !                end do
                        !        end if
                        !end do
                       
                        !diff = 0
                        !do i=1, num_species
                        !        do j=1, species_nums(i)
                        !                if(pool(p,i,j) .and. .not. generation(parent_id2,i,j)) then
                        !                        diff = diff + 1
                        !                end if
                        !        end do
                        !end do
                        !write(*,*) float(diff)/float(diff_start)*100d0

                        !Mutations
                        do i=1, num_species
                                do j=1, species_nums(i)-1
                                        do k=j+1, species_nums(i)
                                                if(pool(p,i,k) .neqv. pool(p,i,j)) then
                                                        call random_number(random)
                                                        if(random .lt. mutation_rate) then
                                                                buffer = pool(p,i,j)
                                                                pool(p,i,j) = pool(p,i,k)
                                                                pool(p,i,k) = buffer
                                                                do l=1, num_species
                                                                        !make double-flip if shared site
                                                                        if(l .ne. i .and. shared_sites(l,i)) then
                                                                                if(pool(p,i,j) .and. &
                                                                                   pool(p,l,j)) then
                                                                                        buffer = pool(p,l,j)
                                                                                        pool(p,l,j) = pool(p,l,k)
                                                                                        pool(p,l,k) = buffer
                                                                                        exit
                                                                                elseif(pool(p,i,k) .and. &
                                                                                       pool(p,l,k)) then
                                                                                        buffer = pool(p,l,j)
                                                                                        pool(p,l,j) = pool(p,l,k)
                                                                                        pool(p,l,k) = buffer
                                                                                        exit
                                                                                end if
                                                                        end if
                                                                end do
                                                        end if
                                                end if
                                        end do
                                end do
                        end do

                        if(.not. valid_solutions(num_species, max_site_num, 1, species_nums, &
                                                 species_occs, shared_sites, pool(p,:,:))) then
                                write(*,*) "ERROR: Invald Solution was considered!"
                                write(*,*) pool(p,:,:)
                                call EXIT(1)
                        end if 
                        pool_energies(p) = energy(num_species, max_site_num, species_nums, &
                                                  pool(p,:,:), const_e, a_e, b_e)
                end do
                !$OMP END PARALLEL DO

                !Genration creation
                generation_new = .FALSE.
                generation_energies_new = HUGE(1.0d+0)
                !Elite
                do i=1, generation_size
                        if(generation_energies(i) .lt. MAXVAL(generation_energies_new(:elite_size),1)) then
                                if(solution_unique(num_species, max_site_num, generation_size, species_nums, &
                                                   generation(i,:,:), generation_new)) then
                                        generation_new(MAXLOC(generation_energies_new(:elite_size),1),:,:) &
                                                        = generation(i,:,:)
                                        generation_energies_new(MAXLOC(generation_energies_new(:elite_size),1)) &
                                                        = generation_energies(i)
                                end if
                        end if
                end do
                !Rest of pool
                do p=1, pool_size
                        if(pool_energies(p) .lt. MAXVAL(generation_energies_new(elite_size+1:),1)) then
                                if(solution_unique(num_species, max_site_num, generation_size, species_nums, &
                                                   pool(p,:,:), generation_new)) then
                                        generation_new(MAXLOC(generation_energies_new(elite_size+1:),1)&
                                                       +elite_size,:,:) = pool(p,:,:)
                                        generation_energies_new(MAXLOC(generation_energies_new(elite_size+1:),1)&
                                                                +elite_size) = pool_energies(p)
                                end if
                        end if
                end do

                if(MAXVAL(generation_energies_new) .ge. HUGE(1.0d+0)) then
                        call random_samples(species_nums, species_occs, shared_sites, const_e, a_e, b_e, generation_size, &
                                            generation_energies, generation, max_site_num, num_species)
                        do i=1, generation_size
                                if(generation_energies_new(i) .ge. HUGE(1.0d+0)) then
                                        generation_new(i,:,:) = generation(i,:,:)
                                        generation_energies_new(i) = generation_energies(i)
                                end if
                        end do
                end if

                generation = generation_new
                generation_energies = generation_energies_new
                if(MODULO(t, write_steps) .eq. 0) then
                        call cpu_time(now)
                        write(*,*) "Generation:", t, "Avg. Fitness:", SUM(generation_energies(:),1)/generation_size, &
                                   "Min:", MINVAL(generation_energies(:),1), "Max:", MAXVAL(generation_energies(:),1), &
                                   "CPU-TIME:", (now-start)
                end if
                if((tol .gt. 0.d0 .and. ABS(MINVAL(best_energies)-prev_min_energy) .lt. tol) .or. & 
                   (tol .lt. 0.d0 .and. MINVAL(best_energies) .ge. prev_min_energy)) then
                        steps_no_improve = steps_no_improve+1
                else
                        steps_no_improve = 0
                        prev_min_energy = MINVAL(best_energies)
                end if
        end do
        return
end subroutine ga

