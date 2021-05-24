import numpy as np
#This function reads the 1 electron Hamiltonian from the disk
#
#On Entry
#file_name--> Name of the file to be read from
#basis--> Dimension of the basis set, has to be explicitly provided, in the present case it is 7
#
#On Exit
#A-->Numpy array with the 1 electron Hamiltonian elements(nbasis,nbasis)
##########################################################
def file_read_1e(file_name,nbasis):

 #open the file
 input_file=open(file_name)   #open the file
 #read the file using readline to file_content
 file_content=input_file.readlines()# read the content
 #close the file
 input_file.close() 
 A=np.zeros([nbasis,nbasis])

 for line in file_content:
    V_line=line.rstrip()
    V_line=V_line.split() 
    i=int(V_line[0])-1
    j=int(V_line[1])-1
    A[i][j]=float(V_line[2])
    A[j][i]=float(V_line[2])
 return A

kinetic_energy = file_read_1e('t.dat',7)
nucear_attraction_intergral = file_read_1e('v.dat',7)

H = kinetic_energy + nucear_attraction_intergral
print("H=",H)
print("k=" ,kinetic_energy)
print("v = ",nucear_attraction_intergral)
