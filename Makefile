CC = icc
FFLAGS = -O3
all: lammps_to_cube

lammps_to_cube: 
	$(CC) $(FFLAGS) lammps_to_cube.c -o lammps_to_cube 

clean:
	rm -f lammps_to_cube
