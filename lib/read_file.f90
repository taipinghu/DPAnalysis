subroutine read_nframe_pdb(pdbname, nframe, natom, cell, coord) 
    implicit none
    integer, intent (in) :: nframe, natom
    character*240, intent(in) :: pdbname
	real(kind=8), dimension(nframe, 6), intent(out) :: cell
    character*2, dimension(nframe, natom) :: symbol
	real(kind=8), dimension(nframe, natom, 3), intent(out) :: coord

    integer :: iframe, itemp, jtemp, iatom
    character*240 :: title, buffer, resname
    character(len=4) tmpstr

    open(10, file=TRIM(pdbname), status="old", position="rewind")
    do iframe=1, nframe
        read(10, *) 
        read(10, *) 
        read(10, *)
        read(10, "(a)") buffer
        !write(*, *) buffer
        if(buffer(1:5) == 'CRYST') then
            read(buffer(8:), *) cell(iframe,:)
        end if
        read(10, *)  ! read MODEL line
        do iatom=1, natom
            read(10, "(a)") buffer
            if(buffer(1:4)=='ATOM' .or. buffer(1:6)=='HETATM') then
                read(buffer, "(a6, i5, 2x, a2, a5, i6, 4x, 3f8.3)") tmpstr, itemp, symbol(iframe, iatom), resname, jtemp, coord(iframe, iatom, 1:3)
                !read(buffer, "(12x, a4, 14x, 3f8.3, 22x, a3)") tmpstr, symbol(iframe, iatom), coord(iframe, iatom, 1:3), element
                !write(*,*) symbol(iframe, iatom)
            end if
        end do
        read(10, *)  ! read TER line
        read(10, *)  ! read ENDMDL line 
    end do

    close(10)
end subroutine read_nframe_pdb
