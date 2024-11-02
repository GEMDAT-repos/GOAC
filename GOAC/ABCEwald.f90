module Ewald
    implicit none
    private     
    public:: getABC
                        
contains

subroutine getABC(n, r, q, lat, alpha, cutoff_real, cutoff_fourier, acc, species_site_map, &
                  num_species_site, max_num_species, const_sites, num_species, max_site_num, c, a, b, alpha_out)
        implicit none
        integer, intent(in):: n, max_num_species, num_species_site(n), num_species, max_site_num, &
                              species_site_map(n,max_num_species,2)
        !f2py intent(in) n, max_num_species, num_species_sites, num_species, max_site_num, &
        !f2py            num_species_map, site_species_map
        logical, intent(in):: const_sites(n)
        !f2py intent(in) const_sites
        real(kind=8), intent(in):: r(n,3), q(n,max_num_species), lat(3,3), acc
        !f2py intent(in) r, q, lat, acc
        real(kind=8), intent(inout):: alpha, cutoff_real, cutoff_fourier
        !f2py intent(inout) alpha, cutoff_real, cutoff_fourier
        integer:: i, ii, j, jj
        real(kind=8):: em(n, max_num_species, n, max_num_species)
        real(kind=8), intent(out):: c, a(num_species, max_site_num), &
                                    b(num_species, max_site_num, num_species, max_site_num), alpha_out
        !f2py intent(out) c, a, b, alpha_out

        c = 0.d0
        a = 0.d0
        b = 0.d0

        call matrix(n, r, q, lat, alpha, cutoff_real, cutoff_fourier, acc, max_num_species, num_species_site, em, alpha_out)
        
        do i=1, n
        do j=1, n
                if(i .ge. j .and. const_sites(i) .and. const_sites(j)) then
                        c = c + em(i,1,j,1)*2.d0
                        if(i .eq. j) c = c - em(i,1,j,1)
                        cycle
                end if
                do ii=1, num_species_site(i)
                do jj=1, num_species_site(j)
                        if((.not.(const_sites(i)) .and. const_sites(j)) .or. (.not.(const_sites(i)) &
                            .and. i .eq. j .and. ii .eq. jj)) then
                                a(species_site_map(i,ii,1), species_site_map(i,ii,2)) = &
                                a(species_site_map(i,ii,1), species_site_map(i,ii,2)) + em(i,ii,j,jj)*2.d0
                                if(i .eq. j .and. ii .eq. jj) then
                                        a(species_site_map(i,ii,1), species_site_map(i,ii,2)) = &
                                        a(species_site_map(i,ii,1), species_site_map(i,ii,2)) - em(i,ii,j,jj)
                                end if
                        elseif((.not.(const_sites(i)) .and. .not.(const_sites(j))) .and. i .ne. j) then
                                b(species_site_map(i,ii,1), species_site_map(i,ii,2), &
                                  species_site_map(j,jj,1), species_site_map(j,jj,2)) = em(i,ii,j,jj)*2.d0!&
                        end if
                end do
                end do
        end do
        end do

        return
end subroutine getABC       

subroutine matrix(n, r, q, lat, alpha, cutoff_real, cutoff_fourier, acc, max_num_species, num_species_site, em, alpha_out)
        !$ use omp_lib
        implicit none
        integer, intent(in):: n, max_num_species, num_species_site(n)
        real(kind=8), intent(in):: r(n,3), q(n,max_num_species), lat(3,3), acc
        real(kind=8), intent(inout):: alpha, cutoff_real, cutoff_fourier
        integer:: i, ii, j, jj, x, y, z, ys, zs, g_counter, k, ds_real(3), ds_fourier(3)
        real(kind=8):: lens(3), rj(3), dist, val_f, val_r, v, invdet, f_lat(3,3), g(3)
        real(kind=8), allocatable:: gs(:,:), g2s(:)
        logical:: first
        real(kind=8):: e_realm(n, max_num_species, n, max_num_species), &
                       e_fourierm(n, max_num_species, n, max_num_species), &
                       e_selfm(n, max_num_species, n, max_num_species)
        real(kind=8), parameter:: ang2bohr=0.52917721067121d0, hartree2eV=27.21138624598130d0, pi=4.d0*ATAN(1.d0)
        real(kind=8), intent(out):: em(n, max_num_species, n, max_num_species), alpha_out

        em = 0.d0
        e_realm = 0.d0
        e_fourierm = 0.d0
        e_selfm = 0.d0
       
        !Get Lattice parameter length and volume 
        do i=1,3
                lens(i) = norm2(lat(i,:))
        end do
        call cross3d(lat(1,:),lat(2,:), g)
        v = DOT_PRODUCT(g ,lat(3,:))

        !Calculate Fourier-Space lattice (inversed transposed lattice)
        invdet = 1/(lat(1,1)*lat(2,2)*lat(3,3) - lat(1,1)*lat(2,3)*lat(3,2)&
                 - lat(1,2)*lat(2,1)*lat(3,3) + lat(1,2)*lat(2,3)*lat(3,1)&
                 + lat(1,3)*lat(2,1)*lat(3,2) - lat(1,3)*lat(2,2)*lat(3,1))
        f_lat(1,1) = +invdet * (lat(2,2)*lat(3,3) - lat(2,3)*lat(3,2))
        f_lat(1,2) = -invdet * (lat(2,1)*lat(3,3) - lat(2,3)*lat(3,1))
        f_lat(1,3) = +invdet * (lat(2,1)*lat(3,2) - lat(2,2)*lat(3,1))
        f_lat(2,1) = -invdet * (lat(1,2)*lat(3,3) - lat(1,3)*lat(3,2))
        f_lat(2,2) = +invdet * (lat(1,1)*lat(3,3) - lat(1,3)*lat(3,1))
        f_lat(2,3) = -invdet * (lat(1,1)*lat(3,2) - lat(1,2)*lat(3,1))
        f_lat(3,1) = +invdet * (lat(1,2)*lat(2,3) - lat(1,3)*lat(2,2))
        f_lat(3,2) = -invdet * (lat(1,1)*lat(2,3) - lat(1,3)*lat(2,1))
        f_lat(3,3) = +invdet * (lat(1,1)*lat(2,2) - lat(1,2)*lat(2,1))
        f_lat = f_lat*2.d0*pi

        if(alpha .le. 0) alpha = (0.005d0*pi**3.d0*n/(v**2.d0))**(1.d0/6.d0)
        alpha_out = alpha
        if(cutoff_real .le. 0) cutoff_real = sqrt(abs(log(acc)))/alpha
        if(cutoff_fourier .le. 0) cutoff_fourier = 2.d0*alpha*sqrt(abs(log(acc)))

        call calc_ds(lat, cutoff_real, ds_real)
        call calc_ds(f_lat, cutoff_fourier, ds_fourier)

        !Calculate g vectors
        g_counter = 1
        allocate(gs((ds_fourier(1)+1)*(ds_fourier(2)+1)*(ds_fourier(3)+1)*4, 3))
        allocate(g2s((ds_fourier(1)+1)*(ds_fourier(2)+1)*(ds_fourier(3)+1)*4))
        gs = 0.d0
        g2s = 0.d0
        do x=0, ds_fourier(1)
                if(x .eq. 0) then
                        ys = 0
                else
                        ys = -ds_fourier(2)
                end if
                do y=ys, ds_fourier(2)
                        if(x .eq. 0 .and. y .eq. 0) then
                                zs = 1
                        else
                                zs = -ds_fourier(3)
                        end if
                        do z=zs, ds_fourier(3)
                                g = f_lat(1,:)*x + f_lat(2,:)*y + f_lat(3,:)*z
                                dist = norm2(g)
                                if(dist .le. cutoff_fourier) then
                                        gs(g_counter,:) = g
                                        g2s(g_counter) = dist**2.d0
                                        g_counter = g_counter + 1
                                end if
                        end do
                end do
        end do

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii, j, jj, k, val_f, val_r, rj, x, y, z, dist)
        do i=1, n
                first = .TRUE.
                do j=i, n
                        !Fourier Space
                        val_f = 0.d0
                        rj = r(i,:)-r(j,:)
                        do k=1, g_counter-1
                                val_f = val_f + exp(-g2s(k)/(4.d0*alpha**2.0))/g2s(k)*cos(sum(gs(k,:)*rj))
                        end do
                        !Real Space
                        val_r = 0.d0
                        do x=-ds_real(1), ds_real(1)
                        do y=-ds_real(2), ds_real(2)
                        do z=-ds_real(3), ds_real(3)
                                if(x .eq. 0 .and. y .eq. 0 .and. z .eq. 0 .and. j .eq. i) cycle
                                rj = r(j,:) + lat(1,:)*x + lat(2,:)*y + lat(3,:)*z
                                dist = norm2(r(i,:)-rj)
                                if(dist .le. cutoff_real) then
                                        val_r = val_r + erfc(alpha*dist)/dist
                                end if 
                        end do
                        end do
                        end do
                        do ii=1, num_species_site(i)
                                !Self Energy
                                if(first) e_selfm(i,ii,i,ii) = -q(i,ii)**2.d0*alpha/sqrt(pi)*ang2bohr*hartree2eV
                                do jj=1, num_species_site(j)
                                        e_fourierm(i,ii,j,jj) = val_f*q(i,ii)*q(j,jj)
                                        e_fourierm(j,jj,i,ii) = e_fourierm(i,ii,j,jj)
                                        e_realm(i,ii,j,jj) = val_r*q(i,ii)*q(j,jj)
                                        e_realm(j,jj,i,ii) = e_realm(i,ii,j,jj)
                                end do
                        end do
                        first = .FALSE.
                end do
        end do
        !$OMP END PARALLEL DO
        e_realm = e_realm*ang2bohr*hartree2eV*0.5d0
        e_fourierm = e_fourierm/v*ang2bohr*hartree2eV*pi*4.d0
        em = e_realm + e_fourierm + e_selfm

        deallocate(gs)
        deallocate(g2s)
        return

        contains
                subroutine calc_ds(lat, cutoff, ds)
                        implicit none
                        real(kind=8), intent(in):: lat(3,3), cutoff
                        real(kind=8):: ab(3), ac(3), bc(3)
                        integer, intent(out):: ds(3)

                        call cross3d(lat(1,:), lat(2,:), ab)
                        ab = ab/norm2(ab)
                        call cross3d(lat(1,:), lat(3,:), ac) 
                        ac = ac/norm2(ac)
                        call cross3d(lat(2,:), lat(3,:), bc)
                        bc = bc/norm2(bc)

                        ds(1) = INT(cutoff/ABS(DOT_PRODUCT(lat(1,:), bc))+1.d0)
                        ds(2) = INT(cutoff/ABS(DOT_PRODUCT(lat(2,:), ac))+1.d0)
                        ds(3) = INT(cutoff/ABS(DOT_PRODUCT(lat(3,:), ab))+1.d0)

                end subroutine calc_ds
   
end subroutine matrix

subroutine cross3d(a,b, cross)
        implicit none
        real(kind=8), intent(in):: a(3), b(3)
        real(kind=8), intent(out):: cross(3)

        cross = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
        return

end subroutine cross3d

end module Ewald
