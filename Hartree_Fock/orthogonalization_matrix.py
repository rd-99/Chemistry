import numpy as np
from numpy import linalg as LA

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

overlap = file_read_1e("s.dat",7)           # overlap matrix(7,7)
a,eigen_vector_ls = LA.eig(overlap)         # calculating eigen values and eigen-vectors

lambda_s = np.zeros((7,7))                  # writing eigen vectors as a diagonal matrix
#inv = np.linalg.inv(np.sqrt(lambda_s))
for i in range(7) :
    lambda_s[i][i] = a[i]

q = eigen_vector_ls.transpose()
p = (np.matmul(eigen_vector_ls,np.linalg.inv(np.sqrt(lambda_s))))

diagonal_S = np.matmul((np.matmul(eigen_vector_ls,lambda_s)),q)  #diagonalzing S


answer = np.matmul(p,q)   # S^(-1/2) = L_s* ((lambda_s)^(-1/2)) *(L_s)~

print(answer)
