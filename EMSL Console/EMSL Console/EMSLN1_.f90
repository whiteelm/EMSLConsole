    program emsln1_
! ���� ������������ emsln1(ip,aa,aw,s,dL,dC)
! EMSLN - Embedded Microstrip Coupled Lines, N-strip 
! ���� - ������������ ��������� ��������������� �����
    implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
    dimension aa(6), aw(9), s(8)		! Nmax=9, ������� n<=4
    dimension bC(81), dC(81), dCC(81), dLL(81), dL(81), em(9), Um(81), dZ0(81)
    data aw/1.,1.,1.,1.,1.,1.,1.,1.,1./ ! ������ ����� �������������� ����� (��)
    data s /1.,1.,1.,1.,1.,1.,1.,1./	! ������ ����������� ������� (��)
	e2=1.	! ����������� ������ e2 < e1. ����������� n=(1..4) 
	n=2; aw(1)=2; aw(2)=2; s(1)=1; e2=1; e1=10; t=.035; h1=2; h2=3
!   e1=10; h1=1; h2=.5; t=.05! ������� (h2/h1)=(0.14..0.22..0.622)..8; e1>e2 
! �������� ������:
	aa(1) = n 	! ���������� ��������������� ����� �������� 
    aa(2) = h1 	! ������� ��������������� �������� (��)
	aa(3) = h2 	! ������� ��������������� �������� (��)
	aa(4) = t 	! ������� ��������������� ����� (��)
    aa(5) = e1 	! ����.���������. ����. �������� h1 � �������� h2
	aa(6) = e2 	! ����.���������. ������� �����.�����, e2=1
    print*,'Input EMSLN1:'; 
	write(*,'('' w'',i1,'' ='',f8.3)')(i, aw(i),i=1,n); 
	write(*,'('' s'',i1,'' ='',f8.3)')(i, s(i),i=1,n-1);
	write(*,1) ,h1,h2,t,e1,e2
1   format(' h1 =',f5.2,'   h2 =',f5.2,'  t=',f5.3/1x,'e1 = ',f5.2,10x,'  e2 =',f5.2/)
	ip=-2  
	call Main(ip,aa,aw,s,dL,dC)  ! ��������� � ������������ 
	print*,'Capacitance matrix [C] (pF/m)';		call DPRINT(dC, n)
    print*,'Inductance matrix [L] (uH/m)';		call DPRINT(dL, n)
    dCC = dC ! ��� ���������� �� ������ nroot(n,dC,bC,em,Um)
    dLL = dL
	call dminv(dLL,n,ad)
	call nroot(n,dCC,11.127*dLL,em,Um)
	print*,'Modal voltage matrix [Um] (V)'; 	call DPRINT(Um, n)
    print*,'Modal dielectric permitivities [em]';
	do i =1,n; print '(f8.2)',em(i); enddo; 
!	if(n==2) print*, 'm = sqrt(em1/em2) =', sqrt(em(1)/em(2));  print*
	if(n==2) print'(a20,f8.4)','m = sqrt(em1/em2) =', sqrt(em(1)/em(2));  print*
	call impedance(n,dC,Um,em,dZ0)
    print*,'Impedance matrix [Z0] (Ohm)'; 	call DPRINT(dZ0, n)
	end