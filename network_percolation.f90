	PROGRAM NETWORK_PERCOLATION
	
!*******************************USE MODULE************************************************

		use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer, C_CHAR
	    use xtc

!*******************************ASSIGN VARIABLES******************************************

		implicit none
		
		real(8), allocatable			::	com(:,:,:)										!Center of mass (molecule A)
		real(8), allocatable			::	car(:,:,:,:)									!Atom Cartesian coordinates	
		real(8), allocatable			::	comb(:,:,:)										!Center of mass (molecule B)
		real(8)							::	l(3)											!length of box, rho for g1, rho for g2
		real(8), allocatable			::	massa(:,:)										!mass of atom
		real(8), allocatable			::	massm(:)										!mass of molecule
		real(8)							::	rcutoff											!radial cutoff
		real(8)							::	rij(3), r, rcomp								!displacement between molecules
		real(8), allocatable			::	rmat(:,:)										!radius matrix
		real(8)							::	mean_cluster, prop_cluster,num, den				!numerator, denominator
		double precision, PARAMETER 	:: 	pi=3.1415926535d+0					
		integer, allocatable			::	nm(:)											!number of molecules
		integer, allocatable			::	na(:)											!number of atoms
		integer, allocatable			::	dum(:)
		integer							::	ii, jj, kk, maxvalue
		integer							::	tstart
		integer							::	pbc, cube
		integer							::	maxmolsp, imolsp								!maximum # of molecular species, current	
		integer							::	omol											!arbitrary molecule number
		integer							::	x,y,z											!dimensions x,y,z
		integer							::	tnm, tmaxstep									!total # of molecules, time steps
		integer 						::	iatom, itatom, xyz, itime						!current atom, current atom, dimensions, current time
		integer							::	imol, imola, imolb								!molecule number
		integer, allocatable			:: 	cluster(:), cluster_size(:)
		integer							::	cluster_count
		integer							::	min_cluster, max_cluster				 
		character*40					::	outfile				

 		character*40					:: infile								
   		real, allocatable				:: pos(:,:)
    	integer 						:: NATOMS, STEP, STAT
    	real 							:: box(3,3), prec, time, box_trans(3,3)
    	type(C_PTR) 					:: xd_c
    	type(xdrfile), pointer 			:: xd
    	logical 						:: ex
    	
!****************************READ SHELL SCRIPT********************************************

		read(5,*)	infile, outfile, maxmolsp	
		read(5,*)	tstart, rcutoff, pbc
		infile = trim(infile)//C_NULL_CHAR
		
!****************************OPEN XTC FILE THROUGH MODULE*********************************		
		
    	inquire(file=trim(infile),exist=ex)

    	STAT = read_xtc_natoms(infile,NATOMS)
    	allocate(pos(3,NATOMS))

    	xd_c = xdrfile_open(infile,"r")
    	call c_f_pointer(xd_c,xd)		
    	
!****************************ALLOCATE, READ, CALCULATE VARIABLES**************************

		allocate(na(maxmolsp))									
		allocate(nm(maxmolsp))									
		allocate(massa(maxmolsp,na(maxmolsp)))									
		allocate(massm(maxmolsp))	
	
		massm = 0								
		do imolsp = 1, maxmolsp								
		  read(5,*) nm(imolsp), na(imolsp)		
		end do  		
			
		read(5,*) tnm

		do imolsp = 1, maxmolsp		
		  do itatom = 1,na(imolsp)						
		    read(5,'(F7.3)') massa(imolsp,itatom)	
		  end do
		end do  		
		
		write(*,*) "Done Reading Mass"
	
		cube = pbc**3
		itime = 0
		mean_cluster = 0
		max_cluster = 0
		min_cluster = tnm
		
		allocate(com(3,maxmolsp,nm(maxmolsp)))
		allocate(comb(3,maxmolsp,cube*nm(maxmolsp)))
		allocate(car(3,maxmolsp,nm(maxmolsp),na(maxmolsp)))
		allocate(rmat(cube*nm(1),cube*nm(1)))
		allocate(cluster(cube*nm(1)))
		allocate(dum(cube*nm(1)))		
		
!******************************OPEN OUTPUT FILE*******************************************
		
		open(unit=11,file=outfile)								
		
!**************************READ ATOM POSITIONS********************************************	
	
		write(*,*) "Reading xtc"
		write(*,*) "..."
		write(*,*) "..."
		write(*,*) "..."

    	STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
        box = transpose(box_trans)
		do while ( STAT == 0 )
			itatom=0
			if (time < tstart) goto 91
			if(time == tstart) then
				l(1) = 10.0d0*box(1,1)
				l(2) = 10.0d0*box(2,2)
				l(3) = 10.0d0*box(3,3)
			end if
						
!******************************ASSEMBLE SEPARATED MOLECULES TOGETHER (PBC)****************
			do imolsp = 1, maxmolsp
				do imol = 1, nm(imolsp)
					do iatom = 1, na(imolsp)
						itatom = itatom + 1
						car(:,imolsp,imol,iatom)=pos(:,itatom)
						do xyz = 1,3
							if ((car(xyz,imolsp,imol,iatom)-car(xyz,imolsp,imol,1)) > (l(xyz)/20)) then
								car(xyz,imolsp,imol,iatom) = car(xyz,imolsp,imol,iatom) - (l(xyz)/10)
							else if ((car(xyz,imolsp,imol,iatom)-car(xyz,imolsp,imol,1)) < (-l(xyz)/20)) then
								car(xyz,imolsp,imol,iatom) = car(xyz,imolsp,imol,iatom) + (l(xyz)/10)
							end if
						end do
						
!******************************CALCULATE COM FOR EACH MOLECULE****************************

					    com(:,imolsp,imol) = com(:,imolsp,imol) + 10.0d0*car(:,imolsp,imol,iatom)*massa(imolsp,iatom)
						if ((time.eq.tstart) .and. (imol.eq.1)) then	
							massm(imolsp) = massm(imolsp) + massa(imolsp,iatom)
						end if				
					end do							
			    	com(:,imolsp,imol) = com(:,imolsp,imol) / massm(imolsp)	
				end do
			end do
		
!**************************PBC: CREATE 8 MORE BOXES***************************************

			do imolsp=1,maxmolsp
				omol=0
				do imol=1,nm(imolsp)
					do xyz=1,3	
				    	if(com(xyz,imolsp,imol) >= l(xyz)) then
							com(xyz,imolsp,imol) = com(xyz,imolsp,imol) - l(xyz)
						else if(com(xyz,imolsp,imol) < 0) then
							com(xyz,imolsp,imol) = com(xyz,imolsp,imol) + l(xyz)
						end if
					end do
					do x=0,(pbc-1)
						do y=0,(pbc-1)
							do z=0,(pbc-1)
								omol=omol+1
								comb(1,imolsp,omol)=com(1,imolsp,imol)+ x*l(1)
								comb(2,imolsp,omol)=com(2,imolsp,imol)+ y*l(2)
								comb(3,imolsp,omol)=com(3,imolsp,imol)+ z*l(3)
							end do
						end do
					end do
				end do
			end do			
	
!##############################CALCULATE MATRIX OF RADII##################################	
							
			do imola = 1, (cube*nm(1))				
				do imolb = 1, (cube*nm(1))
					rij(:) = comb(:,1,imola) - comb(:,1,imolb)
					r = 0
					do xyz=1,3
						rcomp = rij(xyz)**2
						r = r + rcomp
					end do					
					rmat(imola,imolb) = sqrt(r)
				end do
			end do	

			do imola = 1, (cube*nm(1))				
				do imolb = 1, (cube*nm(1))
					if (rmat(imola,imolb) < rcutoff) then
						rmat(imola,imolb) = 1
					else
						rmat(imola,imolb) = 0
					end if
				end do
			end do

!#########################GROUP BY CLUSTERS###############################################
		
			dum = 0
			cluster_count = 1
			cluster = 0
			
			do
				imola = minloc(cluster,DIM=1)
				cluster(imola) = cluster_count
				dum(1) = imola
				kk = 0
				ii = 1
				jj = 1
				
				do 
					do imolb = 1, cube*nm(1)
						kk = 0	
						if (cluster(imolb) > 0) goto 92					
						if (rmat(imola,imolb) > 0) then
							ii = ii + 1
							kk = 1
							dum(ii) = imolb
							cluster(imolb) = cluster_count
						end if
			92			continue
						if (kk == 0) kk = -1
					end do
					if (kk == -1 .and. ii == jj) goto 93
					jj = jj + 1
					imola = dum(jj)
				end do
				dum = 0
			93	continue
				if (minval(cluster) > 0) goto 94
				cluster_count = cluster_count + 1
			end do
94			continue

!##################################COUNT CLUSTER SIZES####################################
		
			maxvalue = maxval(cluster)
			allocate(cluster_size(maxvalue))
			cluster_size = 0 
			do imola = 1, cube*nm(1)
				cluster_size(cluster(imola)) = cluster_size(cluster(imola)) + 1
			end do

			jj=0
			num = 0
			den = 0
			do ii = 1, maxvalue
				num = num + 1/real(maxvalue)*real(cluster_size(ii))**2
				den = den + 1/real(maxvalue)*real(cluster_size(ii))
				jj = jj + real(cluster_size(ii))
			end do
			mean_cluster = mean_cluster + num/den
			if (minval(cluster_size) <= min_cluster) min_cluster = minval(cluster_size,DIM=1)
			if (maxval(cluster_size) >= max_cluster) max_cluster = maxval(cluster_size,DIM=1)
			deallocate(cluster_size)						
			itime=itime+1
			
			
			
!			if (itime==10) goto 101
			
			
			
			if (mod(itime,100) == 0) write(*,*) "Check", itime
91			CONTINUE
    		STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
		end do
		STAT = xdrfile_close(xd)
		write(*,*) "Finished reading"
				
!************************NORMALIZATION****************************************************

		write(*,*) "Calculating"
		write(*,*) "..."
		write(*,*) "..."


!				101		continue
			
		mean_cluster = mean_cluster/itime
		prop_cluster = mean_cluster/(nm(1)*cube)
			
		write(11,'(a25,f10.8)') "Proportion in Cluster: ", prop_cluster
		write(11,'(a25,f10.3)') "Mean Cluster Size: ", mean_cluster
		write(11,'(a25,i10)') "Smallest Cluster Size: ", min_cluster
		write(11,'(a25,i10)') "Largest Cluster Size: ", max_cluster

		write(*,*) "Done!"
			
		close(unit=11)												
			
		deallocate(na)									
		deallocate(nm)									
		deallocate(massa)									
		deallocate(massm)	
		deallocate(com)
		deallocate(comb)
		deallocate(car)
		deallocate(rmat)
		deallocate(cluster)
		deallocate(dum)	
			
		END PROGRAM NETWORK_PERCOLATION	
			
			
			
			
			
			
			
			
		
		
		
		
		