################################################################################
#
#   Function to read the molecular integrals (in mulliken notation)

#
################################################################################
import numpy as np
import time


def read_moint(file_name):
	input_file = open(file_name)  # open the file
	file_content = input_file.readlines()  # read
	input_file.close()  # close
	last_line = file_content[len(file_content) - 1]  # findout  the last line
	last_line = last_line.split()
	o = int(last_line[0])  # the first and the second character of last line must be dimenion of NO and NV
	v = int(last_line[1])
	moint = np.zeros([o, v, o, v])  # alocate the integral
	for line in file_content:
		V_line = line.rstrip()
		V_line = V_line.split()
		i = int(V_line[0]) - 1  # remember python counting starts from zero, so substract one
		a = int(V_line[1]) - 1
		j = int(V_line[2]) - 1
		b = int(V_line[3]) - 1
		moint[i, a, j, b] = float(V_line[4])
	return moint


################################################################################
#
#   Function to read the fock matrix
#
#
################################################################################

def read_fock(file_name):
	input_file = open(file_name)  # open file
	file_content = input_file.readlines()  # read lines
	input_file.close()  # close file
	last_line = file_content[len(file_content) - 1]  # findout  the last line
	last_line = last_line.split()
	nbasis = int(last_line[0])  # the first of last line must be dimenion of nbasis
	fock = np.zeros([nbasis, nbasis])  # allocate Fock
	for line in file_content:
		V_line = line.rstrip()
		V_line = V_line.split()
		i = int(V_line[0]) - 1  # remember python counting starts from zero, so substract one
		j = int(V_line[1]) - 1
		fock[i, j] = float(V_line[2])
	return fock


################################################################################
#
#   Main Program
#
################################################################################
# read Viajb
moint = read_moint('moint')
# read Fock
fock = read_fock('fmatrix')

o = moint.shape[0]  # find dimension of occupied
v = moint.shape[1]  # find dimension of virtual
nbasis = fock.shape[0]  # find dimension of nbasis
emp2 = 0.0
emp2ij = np.zeros((o, o))

#   Converting fock matrix into an efficient format

fock = np.diag(fock)  # Take Fock to be diagonal
fo = fock[:5]
fv = fock[5:]
fia = fo[:, None] - fv[None, :]
fijab = fia[:, None, :, None] + fia[None, :, None, :]

#   calculate MP2 energy

time1 = time.time()  # timer start
t2 = moint.transpose(0, 2, 1, 3).conj() / fijab  # first make the amplitude
emp2 = 2 * np.einsum('ijab,iajb', t2, moint)  # now multipy with direct termp
emp2 -= np.einsum('ijab,ibja', t2, moint)  # Now multiply with exchange term
time2 = time.time()  # timer end

print('MP2 energy using Python', emp2)
print('Time taken by python  {0:.6f} seconds.'.format(time2 - time1))
