#The fortran compiler
FORTRAN_COMPI = ifort#ifort  

#The compilation flags
FORTRAN_OPT_FLAGS =   -i-dynamic -mcmodel=medium  -openmp -unroll4

#Target file: Dependencies files
LES_autorun:	Basic_module.o\
                continuity_module.o\
		x_momentum_module.o\
		y_momentum_module.o\
		z_momentum_module.o\
		energy_module.o\
		other_flow_data_module.o\
		filtering_module.o\
		arithmetic.o\
		flowthrough.o\
		NS_cartesian_impingingjet_main.o\
		initialize.o\
		hixon_scheme.o\
                turkel_scheme.o\
		vorticity_calculation.o\
                temperature_calculation.o\
		mean_calculation.o\
		mean_square_calculation.o\
		spl_calculation.o\
                kirchoff.o\
                lighthill_source.o\
		write.o\
		write_mean.o\
		write_mean_xyprofile.o\
		write_xyprofile.o\
		write_xzprofile.o\
		write_yzprofile.o\
		write_centerline.o\
		read_data.o\
                read_mean.o\
		filter_x.o\
		filter_y.o\
		filter_z.o\
		Boundary1.o\
		Differentiation1.o\
		Boundary_corner1.o\
		Boundary_edge1.o\
		Boundary2.o\
 		Differentiation2.o\
		Boundary_corner2.o\
		Boundary_edge2.o\
		NS_cartesian_RK2_STEP1.o\
		NS_cartesian_RK2_STEP2.o\
		NS_cartesian_RK2.o\
		bufferlayer_filter.o

#Command line
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} \
	-o LES_autorun\
		Basic_module.o\
                continuity_module.o\
		x_momentum_module.o\
		y_momentum_module.o\
		z_momentum_module.o\
		energy_module.o\
		other_flow_data_module.o\
		filtering_module.o\
		arithmetic.o\
		flowthrough.o\
		NS_cartesian_impingingjet_main.o\
		initialize.o\
		hixon_scheme.o\
		turkel_scheme.o\
		vorticity_calculation.o\
                temperature_calculation.o\
		mean_calculation.o\
		mean_square_calculation.o\
		spl_calculation.o\
                kirchoff.o\
                lighthill_source.o\
		write.o\
		write_mean.o\
		write_mean_xyprofile.o\
		write_xyprofile.o\
		write_xzprofile.o\
		write_yzprofile.o\
		write_centerline.o\
		read_data.o\
                read_mean.o\
		filter_x.o\
		filter_y.o\
		filter_z.o\
		Boundary1.o\
		Differentiation1.o\
		Boundary_corner1.o\
		Boundary_edge1.o\
		Boundary2.o\
 		Differentiation2.o\
		Boundary_corner2.o\
		Boundary_edge2.o\
		NS_cartesian_RK2_STEP1.o\
		NS_cartesian_RK2_STEP2.o\
		NS_cartesian_RK2.o\
		bufferlayer_filter.o

Basic_module.o: Basic_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Basic_module.o Basic_module.f90

continuity_module.o: continuity_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o continuity_module.o continuity_module.f90\

x_momentum_module.o: x_momentum_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o x_momentum_module.o x_momentum_module.f90

y_momentum_module.o: y_momentum_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o y_momentum_module.o y_momentum_module.f90

z_momentum_module.o: z_momentum_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o z_momentum_module.o z_momentum_module.f90

energy_module.o: energy_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o energy_module.o energy_module.f90

other_flow_data_module.o: other_flow_data_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o other_flow_data_module.o other_flow_data_module.f90

filtering_module.o: filtering_module.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o filtering_module.o filtering_module.f90

arithmetic.o: arithmetic.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o arithmetic.o arithmetic.f90

flowthrough.o: flowthrough.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o flowthrough.o flowthrough.f90

NS_cartesian_impingingjet_main.o:  NS_cartesian_impingingjet_main.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o  NS_cartesian_impingingjet_main.o  NS_cartesian_impingingjet_main.f90

initialize.o: initialize.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o initialize.o initialize.f90

hixon_scheme.o: hixon_scheme.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o hixon_scheme.o hixon_scheme.f90

turkel_scheme.o: turkel_scheme.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o turkel_scheme.o turkel_scheme.f90

vorticity_calculation.o: vorticity_calculation.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o vorticity_calculation.o vorticity_calculation.f90

temperature_calculation.o: temperature_calculation.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o temperature_calculation.o temperature_calculation.f90

mean_calculation.o: mean_calculation.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o mean_calculation.o  mean_calculation.f90

mean_square_calculation.o: mean_square_calculation.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o mean_square_calculation.o mean_square_calculation.f90

spl_calculation.o: spl_calculation.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o spl_calculation.o spl_calculation.f90

kirchoff.o: kirchoff.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o kirchoff.o kirchoff.f90

lighthill_source.o: lighthill_source.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o lighthill_source.o lighthill_source.f90

write.o: write.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write.o write.f90

write_mean.o: write_mean.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_mean.o write_mean.f90

write_mean_xyprofile.o: write_mean_xyprofile.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_mean_xyprofile.o write_mean_xyprofile.f90

write_xyprofile.o: write_xyprofile.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_xyprofile.o write_xyprofile.f90

write_xzprofile.o: write_xzprofile.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_xzprofile.o write_xzprofile.f90

write_yzprofile.o: write_yzprofile.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_yzprofile.o write_yzprofile.f90

write_centerline.o: write_centerline.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o write_centerline.o write_centerline.f90

read_data.o: read_data.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o read_data.o read_data.f90

read_mean.o: read_mean.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o read_mean.o read_mean.f90

filter_x.o: filter_x.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o filter_x.o filter_x.f90

filter_y.o: filter_y.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o filter_y.o filter_y.f90

filter_z.o: filter_z.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o filter_z.o filter_z.f90

Boundary1.o: Boundary1.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary1.o Boundary1.f90

Differentiation1.o: Differentiation1.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Differentiation1.o Differentiation1.f90

Boundary_corner1.o: Boundary_corner1.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary_corner1.o Boundary_corner1.f90

Boundary_edge1.o: Boundary_edge1.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary_edge1.o Boundary_edge1.f90

Boundary2.o: Boundary2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary2.o Boundary2.f90

Differentiation2.o: Differentiation2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Differentiation2.o Differentiation2.f90

Boundary_corner2.o: Boundary_corner2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary_corner2.o Boundary_corner2.f90

Boundary_edge2.o: Boundary_edge2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o Boundary_edge2.o Boundary_edge2.f90

NS_cartesian_RK2_STEP1.o: NS_cartesian_RK2_STEP1.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o NS_cartesian_RK2_STEP1.o NS_cartesian_RK2_STEP1.f90

NS_cartesian_RK2_STEP2.o: NS_cartesian_RK2_STEP2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o NS_cartesian_RK2_STEP2.o NS_cartesian_RK2_STEP2.f90

NS_cartesian_RK2.o: NS_cartesian_RK2.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o NS_cartesian_RK2.o NS_cartesian_RK2.f90

bufferlayer_filter.o: bufferlayer_filter.f90
	${FORTRAN_COMPI} ${FORTRAN_OPT_FLAGS} -c -o bufferlayer_filter.o bufferlayer_filter.f90


#Cleaning the directory
clean:
	rm *.o *.mod

			
#rm -f *.mod
