#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP CHyM.
#
#    ICTP CHyM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP CHyM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include $(CHYM_ROOT)/makeinc
if CHYM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif

vpath %.o $(CHYM_ROOT)/share:.

## FCFLAGS += -I. `nf-config --fflags`
FCFLAGS += -I. $(NC_INCLUDES)

bin_PROGRAMS = preproc

preproc_SOURCES = preproc.F90

FCLD = $(MPIFC)

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

preproc.o: preproc.F90 mod_internal.o mod_runparams.o mod_param.o mod_mssg.o mod_crtstatic.o mod_ncio.o mod_varandtypes.o mod_mpimess.o
	 $(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
