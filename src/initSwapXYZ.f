c swap XYZ0 to XYZ1
c      make XYZ1 equal to XYZ0
	subroutine initSwap01(n,xyz0,xyz1)
c
	implicit none
	integer n
	real xyz0(*),xyz1(*)
        integer i
c
        do i=1,3*n
        xyz1(i) = xyz0(i)
        end do
c
        return
	end
