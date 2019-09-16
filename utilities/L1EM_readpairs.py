import cPickle as pickle

"""
Report the total numbr of read pairs passed to L1EM

Copyright (C) 2019 Wilson McKerrow

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

total = 0
for line in open('G_of_R_list.txt'):
	G_of_R = pickle.load(open(line.strip()))
	if G_of_R != None:
		total += G_of_R.shape[0]

print total
