import numpy as np

def read2e(file,nbasis):

 #read the file
 #define the 4 d array
 #open the file
 input_file=open(file)   #open the file
 #read the file using readline to file_content
 file_content=input_file.readlines()# read the content
 #close the file
 input_file.close() 
 twoe=np.zeros([nbasis,nbasis,nbasis,nbasis])
 for line in file_content:
    V_line=line.rstrip()
    V_line=V_line.split() 
    i=int(V_line[0])-1  #array index starts from 0
    j=int(V_line[1])-1
    k=int(V_line[2])-1
    l=int(V_line[3])-1    
    
    twoe[i][j][k][l]=float(V_line[4]) #### 8-fold symmetry ###
    twoe[i][j][l][k]=float(V_line[4])
    twoe[j][i][k][l]=float(V_line[4])
    twoe[j][i][l][k]=float(V_line[4])
    twoe[k][l][i][j]=float(V_line[4])
    twoe[l][k][i][j]=float(V_line[4])
    twoe[k][l][j][i]=float(V_line[4])
    twoe[l][k][j][i]=float(V_line[4])
 return twoe

 EMP2=EMP2+2*vijab[i][j][a][b]*vijab[i][j][a][b]/(fock[i][i]+fock[j][j]-fock[no+a][no+a]-fock[no+b][no+b])


#print the output
a=read2e("eri.dat",7)
for i in range(7):     
    for j in range(7):
        for k in range(7):
            for l in range(7):
                print(i+1,j+1,k+1,l+1,a[i][j][k][l])

