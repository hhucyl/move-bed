########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
# Copyright (C) 2009 Sergio Galindo                                    #
# Copyright (C) 2013 William Oquendo                                   #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

PROJECT                (MECHSYSSANDBOX)
CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
INCLUDE(${MECHSYSSANDBOX_SOURCE_DIR}/UseMECHSYS.cmake)

SET(PROGS
    test_circle_move
    test_circle_spin
    test_circle2
    test_collision
    test_periodic
    test_periodic1
    test_periodic2
    test_mvbed
    test_mvbed_compact
    test_mvbedrw
    test_particle_reflect
    test_cd
)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)
