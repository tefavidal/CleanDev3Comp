all:SimGoldbeter
SimGoldbeter:	
	gfortran -o program main-Source.f anfang.f rs-Source2D.f \
	out.f ODE-Merson.f ic-Source.f Development.f StartingTime.f flow.f \

