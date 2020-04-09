import time
import random
import math
import matplotlib.pyplot as plt
import numpy
start=time.time()

def matrix_man1(nlin, ncol):
    matrizA = []
    for i in range (nlin):
        linha=[]
        for j in range (ncol):
            linha.append(1)
        matrizA.append(linha)
   
   
    return matrizA

def matrix_man2(lin, col):
    matrizA = []
    for i in range (lin):
        linha=[]
        for j in range (col):
            if (i+j)%2==0:
                linha.append(-1)
            else:
                linha.append(1)
        matrizA.append(linha)

    return matrizA

def matrix_man1(nlin, ncol):
    matrizA = []
    for i in range (nlin):
        linha=[]
        for j in range (ncol):
            linha.append(-1)
        matrizA.append(linha)
   
   
    return matrizA

def matrix_rand(lin, col):
    matrizA = []
    for i in range (lin):
        linha=[]
        for j in range (col):
            r=random.randint(0,1)
            if r==0:
                r-=1
            linha.append(r)
        matrizA.append(linha)
    return matrizA

   
def print_matrix(matrix):
    lin=len(matrix)
    col=len(matrix[0])
    for i in range (lin):
        for j in range (col):
            print("%3s"%matrix[i][j],end=" ")
        print()



def energy_2(matrix,i,j,H):
    lin=len(matrix)
    col=len(matrix[0])
   
    elemento=0
    if i==0:
        if j==0:
            elemento=-(matrix[i][j]*matrix[lin-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][col-1])
        if j==col:
            elemento=-(matrix[i][j]*matrix[lin-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][1]+matrix[i][j]*matrix[i][j-1])
        if j!=0 and j!=col-1:
            elemento=-(matrix[i][j]*matrix[lin-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][j-1])
           
               
    if i==lin-1:
        if j==0:
            elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][col-1])
        if j==col-1:
            elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[1][j]+matrix[i][j]*matrix[i][1]+matrix[i][j]*matrix[i][j-1])
        else:
            elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][j-1])
    elif j==0:
        elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][col-1])
    elif j==col-1:
        elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][j-1]+matrix[i][j]*matrix[i][1])
    if i!=0 and i !=lin-1 and j!=col-1 and j!=0:
        elemento=-(matrix[i][j]*matrix[i-1][j]+matrix[i][j]*matrix[i+1][j]+matrix[i][j]*matrix[i][j+1]+matrix[i][j]*matrix[i][j-1])
   
    return elemento

def energy_infinity(matrix,L):
    m=M_tot(matrix)
    elemento=-1/2*m**2/(L**2)
    return elemento,m

def calc_distance(matrix,m,n,i,j):
    r=m-i
    q=n-j
    r=abs(r)
    q=abs(q)
    if r>lin/2:
        r=lin-r
    if q>lin/2:
        q=lin-q
    d=r+q
    return d
 


def rename_spins(i,j,matrix):
    L=len(matrix)
    elemento=i*L+j
    return elemento

def rerename_spins(elemento,matrix):
    L=len(matrix)
    if elemento!=0:
        i=elemento//L
        j=elemento-(i*L)
    else:
        i=0
        j=0
    return i,j

def belong_to_A_but_not_to_B(A,B):
    a=len(A)
    b=len(B)
    C=[]
    if b!=0:
        for i in range(a):
            belong=False
            for j in range(b):
                if A[i]==B[j]:
                    belong=True
            if belong==False:
                C.append(A[i])
    else:
        C=A
    return C


def a_belong_V(V,a):
    L=len(V)
    belong=False
    for i in range(L):
        if V[i]==a:
            belong=True

    return belong

def A_equal_B(A,B):
    m=len(A)
    n=len(B)
    equal=True
    if m==n:
        for i in range(m):
            if A[i]!=B[i]:
                equal=False
    else:
        equal=False
    return equal

def l(i,j,matrix,T):
    l = False
    m,n=rerename_spins(i,matrix)
    t,r=rerename_spins(j,matrix)
    if matrix[m][n]==matrix[t][r]:
        p=random.random()
        if p <= 1-math.exp(-2/T):
            l = True

    return l

def Wolff_criterion(m,n,matrix,A,T):
    if l(m,n,matrix,T)==True and A[m]==-1 and A[n]!=-1 and A[n]!=-2:
        A[n]=-1
    if l(m,n,matrix,T)==True and A[n]==-1 and A[m]!=-1 and A[m]!=-2:
        A[m]=-1
   
def print_matrix_simbolically(matrix):
    n=len(matrix)
    m=len(matrix[0])
    for i in range(n):
        for j in range(m):
            if matrix[i][j]==-1:
                matrix[i][j]="."
                print("%3s" %matrix[i][j],end=" ")
                matrix[i][j]=-1
            if matrix[i][j]==1:
                matrix[i][j]="O"
                print("%3s" %matrix[i][j],end=" ")
                matrix[i][j]=1
        print()
    print()
   
def Wolff_Algorithm_neighbor(matrix,T,H,steps):
    flip=0
    L=len(matrix)
    i=random.randint(0,L-1)
    j=random.randint(0,L-1)
    A=numpy.zeros(L**2)
    a=len(A)
    caboo=False
    for m in range (a):
        A[m]+=m
    A[rename_spins(i,j,matrix)]=-1
    tem=True
    while tem==True:
        tem=False
        for m in range(L**2):
            if A[m]==-1:
                i,j=rerename_spins(m,matrix)
                if i==0:
                    u=rename_spins(L-1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i+1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    if j==0:
                        u=rename_spins(i,L-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j+1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                    if j==L-1:
                        u=rename_spins(i,0,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                    if j!=0 and j!=L-1:
                        u=rename_spins(i,j+1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                if i==L-1:
                    u=rename_spins(0,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i-1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    if j==0:
                        u=rename_spins(i,L-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j-+1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                    if j==L-1:
                        u=rename_spins(i,0,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                    if j!=L-1 and j!=0:
                        u=rename_spins(i,j+1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                        u=rename_spins(i,j-1,matrix)
                        Wolff_criterion(u,m,matrix,A,T)
                if j==L-1 and i!=0 and i!=L-1:
                    u=rename_spins(i,0,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i,j-1,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i+1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i-1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)

                if j==L-1 and i!=0 and i!=L-1:
                    u=rename_spins(i,j+1,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i,L-1,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i+1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i-1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)

                   
                if i!=0 and j!=0 and i!=L-1 and j!=L-1:
                   
                    u=rename_spins(i,j+1,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i,j-1,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i+1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                    u=rename_spins(i-1,j,matrix)
                    Wolff_criterion(u,m,matrix,A,T)
                A[m]=-2
        for n in range (L**2):
            if A[n]==-1:
                tem=True
    for t in range (L**2):
        if (A[t]==-2 or A[t]==-1):
            i,j=rerename_spins(t,matrix)
            matrix[i][j]=-matrix[i][j]
            flip+=1
    return flip
   
def Wolff_Algorithm_full_connected(matrix,T):
    Y=len(matrix)
    i=random.randint(0,L-1)
    j=random.randint(0,L-1)
    A=[]
    C=[]
    A.append(rename_spins(matrix,i,j))
    while A!=C:
        D=belong_to_A_but_not_to_B(A,C)
        d=len(D)
        for m in range (d):
            C.append(D[m])
            for n in range (d):
                i,j=rerename_spins(matrix,D[n])
                for p in range (Y):
                    for q in range (Y):
                        if matrix[p][q]==matrix[i][j] and a_belong_V(A,rename_spins(matrix,p,q))==False :
                            p=random.uniform(0,1)
                            if p<1-mat.exp(-1/T):
                                A.append(rename_spins(matrix,p,q))
               
    L=len(A)
    for m in range (L):
        i,j=rerename_spins(matrix,A[m])
        matrix[i][j]=-matrix[i][j]
    return L
       
def E_tot(matrix,d,H):
    lin=len(matrix)
    col=len(matrix[0])
    E_tot=0
    if d==2:
        for i in range(lin):
            for j in range (col):
                E_tot+=energy_2(matrix,i,j,H)
    if d==8:
        E_tot,s=energy_infinity(matrix,L)
    if H!=0:
        E_tot-=H*M_tot(matrix)
    return E_tot/4

def M_tot(matrix):
    lin=len(matrix)
    col=len(matrix[0])
    M_tot=0
    for i in range(lin):
            for j in range (col):
                M_tot+=matrix[i][j]
    return M_tot


def cor_func_T2(matrix):
    L=len(matrix)
    J=L/2
    J=J/1
    sum=0
    n=0
    cor=[]
    for j in range (10):
        for i in range (11):
            sum+=matrix[j][j]*matrix[i+j][10-i+j]
            n+=1
    sum=sum/n
    sum=sum-M_tot(matrix)**2/L**4
   
    return sum

def flip(T_i,T_f,delta_T,m,steps,L,d,method,IC,CORR,H,board,po):
    if steps==1:
        n_1=int(input("Type the quantity of Monte Carlo steps per site (T between 1.5 and 2): "))
        n_2=int(input("Type the quantity of Monte Carlo steps per site (T between 2.5 and 3): "))
    if steps==2:
        n=int(input("Type the quantity of Monte Carlo steps per site: "))
    x1=[]
    y1=[]
    y2=[]
    x2=[]
    x3=[]
    y11=[]
    y22=[]
    y33=[]
    y44=[]
    y55=[]
    cor=0
    E_F=0
    i=1
    T=T_i
    EP=[]
    if IC==1:
        matrix1=matrix_man1(L,L)
    if IC==2:
        matrix1=matrix_rand(L,L)
    if IC==3:
        matrix1=matrix_man2(L,L)
    if board=="yes":
        for i in range (L):
            matrix1[0][i]=-1
        for j in range (L):
            matrix1[L-1][j]=-1
        for i in range (L):
            matrix1[i][0]=-1
        for j in range (L):
            matrix1[j][L-1]=-1
    L=len(matrix1)
    print_matrix_simbolically(matrix1)
    flips=0
    if steps==1 and method=="Metropolis":
        if T<2:
            n=n_1
        if (T>2 and T<2.5)or T==2 or T==2.5:
            n=5*L**2.166
        if T>2.5:
            n=n_2
        n=n*L**2
    if method=="Wolff":
        n=n*L**2
    while T<=T_f:
        p0=po
        if steps==1 and method=="Metropolis":
            if T<2:
                n=n_1
            if (T>2 and T<2.5) or T==2 or T==2.5:
                n=L**2.166
            if T>2.5:
                n=n_2
            n=n*L**2
        if method=="Wolff" and T==1.5:
            n=n*L**2

        E2=0
        E1=0
        M1=0
        M2=0
        matrix=matrix1
        for p in range (m):
            flips=0
            while flips < n:
                i=random.randint(0,L-1)
                j=random.randint(0,L-1)
                if board=="yes":
                    i=random.randint(1,L-2)
                    j=random.randint(1,L-2)
                if method=='Metropolis':
                    if d==2:
                        E_i=energy_2(matrix,i,j,H)
                        delta_e=-2*E_i+H*matrix[i][j]
                       
                           
                        if delta_e<0:
                            matrix[i][j]=-matrix[i][j]
                        if delta_e>0:
                            p=random.uniform(0,1)
                            if p<math.exp(-delta_e/T):
                                matrix[i][j]=-matrix[i][j]
                        if p0>0:
                            p0=p0-1
                        else:
                            flips+=1
                    if d==8:
                        delta_e=0
                        E_i,MM=energy_infinity(matrix,L)
                        if matrix[i][j]==1 and MM!=0:
                            E_f=E_i*(MM-1)**2/MM**2
                        if matrix[i][j]==-1 and MM!=0:
                            E_f=E_i*(MM+1)**2/MM**2
                        if MM==0:
                            E_f=-1*(L**2)
                        delta_e=E_f-E_i
                        if delta_e<0:
                            matrix[i][j]=-matrix[i][j]
                        if delta_e>0:
                            p=random.uniform(0,1)
                            if p<math.exp(-delta_e/T):
                                matrix[i][j]=-matrix[i][j]
                        if p0>0:
                            p0=p0-1
                        else:
                            flips+=1
                       
                if method=='Wolff':
                    if d==2:
                        fli=Wolff_Algorithm_neighbor(matrix,T,H,flips)
                        if p0>0:
                            p0=0
                        else:
                            flips+=n

                       
                    if d==8:
                        fli=Wolff_Algorithm_full_connected(matrix,T)
                        if p0>0:
                            p0=p0-1
                        else:
                            flips+=1
           
            M11=M_tot(matrix)
            E11=E_tot(matrix,d,H)
            M22=M_tot(matrix)*M_tot(matrix)
            E22=E_tot(matrix,d,H)*E_tot(matrix,d,H)
            E1+=E11
            E2+=E22
            M1+=M11
            M2+=M22


           
        if CORR=="yes":
            C_length=cor_func_T2(matrix)
        E1=E1/m
        E2=E2/m
        M1=M1/m
        M2=M2/m
        if CORR=="yes":
            C_length=abs(C_length)
        S=(M2-M1*M1)/(T*L**2)
        C=E2/(L**2*T**2)-E1**2/(L**2*T**2)
        varm=(M2-M1**2)/L**2
        vare=(E2-E1**2)/L**2
        y11.append(M1/L**2)
        y22.append(E1/L**2)
        y33.append(C)
        y44.append(S)
        if CORR=="yes":
            y55.append(C_length)
        x2.append(T)
        x3.append(T)
        if CORR=="yes":
            print(C_length)
        print("The specific heat: ")
        print(C)
        print("The magnetization per spin and its variance: ")
        print(abs(M1/(L**2)))
        print(varm)
        print("The energy per spin and its variance: ")
        print(E1/L**2)
        print(vare)
        print("The suscebility: ")
        print(S)
        print("The temperature: ")
        print(T)
        T+=delta_T
        print(".")
       
   
    end=time.time()
    print(end-start)
    plt.plot(x2, y11,'ro',color="Blue")
    plt.xlabel('Temperature')
    plt.ylabel('M(Blue)')
    plt.title('Temperature x M')
    plt.show()
    plt.plot(x2, y22,'ro',color="Yellow")
    plt.xlabel('Temperature')
    plt.ylabel('H(yellow)')
    plt.title('Temperature x H')
    plt.show()
    plt.plot(x2, y33,'ro',color="Red")
    plt.xlabel('Temperature')
    plt.ylabel('C (Red)')
    plt.title('Temperature x C')
    plt.show()
    plt.plot(x2, y44,'ro',color="Green")
    plt.xlabel('Temperature')
    plt.ylabel('S (green)')
    plt.title('Temperature x S')
    plt.show()
    if CORR=="yes":
        plt.plot(x2, y55, 'ro',color="Brown")
        plt.xlabel('Temperature')
        plt.ylabel('Correlation length (brown)')
        plt.title('Temperature x Correlation Length')
        plt.show()




def relaxation_function_noneq(T,method,IC,L,steps,H,board,m,d):
    M1=0
    M2=0
    M3=0
    matrix1=matrix_man1(L,L)
    x=[]
    y=[]
    matrix=matrix1
    steps2=steps
    steps=steps*L**2
    stepsi=steps
    M_initial=M_tot(matrix1)
    for u in range (1,steps2):
        if u%1==0:
            print(M_initial)
            print(M_tot(matrix))
            print()
           
            for n in range (m):
                steps=L**2*u
                matrix=matrix1
                while steps>0:
                    i=random.randint(0,L-1)
                    j=random.randint(0,L-1)
                    if board=="yes":
                        i=random.randint(1,L-2)
                        j=random.randint(1,L-2)
                    if method=='Metropolis':
                        if d==2:
                            E_i=energy_2(matrix,i,j,H)
                            delta_e=-2*E_i+H*matrix[i][j]      
                            if delta_e<0:
                                matrix[i][j]=-matrix[i][j]
                            if delta_e>0:
                                p=random.uniform(0,1)
                                if p<math.exp(-delta_e/T):
                                    matrix[i][j]=-matrix[i][j]
                            steps-=1  
                        if d==8:
                            delta_e=0
                            E_i,MM=energy_infinity(matrix,L)
                            if matrix[i][j]==1 and MM!=0:
                                E_f=E_i*(MM-1)**2/MM**2
                            if matrix[i][j]==-1 and MM!=0:
                                E_f=E_i*(MM+1)**2/MM**2
                            if MM==0:
                                E_f=-1*(L**2)
                            delta_e=E_f-E_i
                            if delta_e<0:
                                matrix[i][j]=-matrix[i][j]
                            if delta_e>0:
                                p=random.uniform(0,1)
                                if p<math.exp(-delta_e/T):
                                    matrix[i][j]=-matrix[i][j]
                            steps-=1
                    if method=='Wolff':
                        if d==2:
                            fli=Wolff_Algorithm_neighbor(matrix,T,H,steps)
                            steps-=L**2
                           

                           
                        if d==8:
                            fli=Wolff_Algorithm_full_connected(matrix,T)
                            flips+=fli
                M1+=abs(M_tot(matrix))
                M2+=abs(M_initial)
                M3+=M_tot(matrix)**2
            M1=M1/m
            M2=M2/m
            M3=M3/m
        if u%1==0 and u!=0:
            MM=abs(M1)/abs(M2)
            varm=(M3-M1**2)/L**4
            y.append(abs(MM))
            x.append(u)
            print("Absolute value of the non-linear relaxation function, the variance of the absolute value of the magnetization and its time: ")
            print(abs(MM))
            print(math.sqrt(abs(varm)))
            print(u)
            print()

    TA=0
    leng=len(y)
    for i in range (leng):
        TA+=y[i]
    TA=TA
    print(TA)
    end=time.time()
    print("Time of the run: ")
    print(end-start)
    plt.plot(x, y, 'ro',color="Brown")
    plt.xlabel('Time in MC steps')
    plt.ylabel('Relaxation function')
    plt.title('Time x Relaxation function')
    plt.show()    
       
       

def relaxation_function_eq(T,method,IC,L,steps,H,board,m,d,po):
    VARMET=[]
    VARWOLF=[]
    matrix=matrix_man1(L,L)
    x=[]
    y=[]
    M1=0
    M2=0
    M3=0
    steps2=steps
    steps=steps*L**2
    stepsi=steps
    po=po*L**2
    while po!=0:
        i=random.randint(0,L-1)
        j=random.randint(0,L-1)
        E_i=energy_2(matrix,i,j,H)
        delta_e=-2*E_i+H*matrix[i][j]      
        if delta_e<0:
            matrix[i][j]=-matrix[i][j]
        if delta_e>0:
            p=random.uniform(0,1)
            if p<math.exp(-delta_e/T):
                matrix[i][j]=-matrix[i][j]
        po-=1
    matrix2=matrix
    M_initial=M_tot(matrix)
    for u in range (1,steps2):
        if u%1==0:
            print(M_initial)
            print(M_tot(matrix))
            print()
            for n in range (m):
                steps=L**2*u
                matrix=matrix2
                while steps>0:
                    i=random.randint(0,L-1)
                    j=random.randint(0,L-1)
                    if board=="yes":
                        i=random.randint(1,L-2)
                        j=random.randint(1,L-2)
                    E_i=energy_2(matrix,i,j,H)
                    delta_e=-2*E_i+H*matrix[i][j]      
                    if delta_e<0:
                        matrix[i][j]=-matrix[i][j]
                    if delta_e>0:
                        p=random.uniform(0,1)
                        if p<math.exp(-delta_e/T):
                            matrix[i][j]=-matrix[i][j]
                    steps-=1
                M1+=M_tot(matrix)*M_initial
                M2+=abs(M_tot(matrix))**2
                M3+=abs(M_tot(matrix))
            M1=M1/m
            M2=M2/m
            M3=M3/m
        if u%1==0:
                MM=(M1)/L**4
                y.append(abs(MM))
                x.append(u)
                varm=(M2-M3**2)/L**4
                VARMET.append(varm)
                print("Absolute value of the linear relaxation function, the variance of the absolute value of the magnetization and its time: ")
                print(abs(MM))
                print(math.sqrt(abs(varm)))
                print(u)
                print()
   
    TA=0
    VARFM=0
    leng=len(y)
    lent=len(VARMET)
    for j in range (lent):
        VARFM+=VARMET[j]**2
    for i in range (leng):
        TA+=y[i]
    TA=TA
    print(math.sqrt(VARFM))
    print(TA)
    end=time.time()
    print(end-start)
    plt.plot(x, y, 'ro',color="Brown")
    plt.xlabel('Temperature')
    plt.ylabel('Relaxation function (brown)')
    plt.title('Temperature x Relaxation function(M)')
    plt.show()

    matrix=matrix2
    M1=0
    M2=0
    M3=0
    x=[]
    y=[]
    print('Wolff')
    for u in range (1,steps2):
        if u%1==0:
            print(M_initial)
            print(M_tot(matrix))
            print()
            for n in range (m):
                steps=L**2*u
                matrix=matrix2
                while steps>0:              
                    fli=Wolff_Algorithm_neighbor(matrix,T,H,steps)
                    steps-=L**2              
                       
                M1+=M_tot(matrix)*M_initial
                M2+=abs(M_tot(matrix))**2
                M3+=abs(M_tot(matrix))
            M1=M1/m
            M2=M2/m
            M3=M3/m
           
        if u%1==0:
            MM=(M1)/L**4
            y.append(abs(MM))
            x.append(u)
            varm=(M2-M3**2)/L**4
            VARWOLF.append(varm)
            print("Absolute value of the linear relaxation function, the variance of the absolute value of the magnetization and its time: ")
            print(abs(MM))
            print(math.sqrt(abs(varm)))
            print(u)
            print()

    TA=0
    VARFW=0
    leng=len(y)
    print(leng)
    lent=len(VARWOLF)
    for j in range (lent):
        VARFW+=VARWOLF[j]**2
    for i in range (leng):
        TA+=y[i]
    TA=TA
    print(TA)
    print(math.sqrt(VARFW))
    end=time.time()
    print(end-start)
    plt.plot(x, y, 'ro',color="Brown")
    plt.xlabel('Temperature')
    plt.ylabel('Relaxation function (brown)')
    plt.title('Temperature x Relaxation Function(W)')
    plt.show()


   
def hysteresis(H_0,H_f,delta_H,n,T,method,IC,L):
    flips=0
    x=[]
    y=[]
    H=H_0
    if IC==1:
        matrix1=matrix_man1(L,L)
    if IC==2:
        matrix1=matrix_rand(L,L)
    if IC==3:
        matrix1=matrix_man2(L,L)
    n=n*L**2
    while H<=H_f:
        matrix=matrix1
        flips=0
        while flips<n:
            if method == 'Metropolis':
                i=random.randint(0,L-1)
                j=random.randint(0,L-1)
                E_i=energy_2(matrix,i,j,H)
                delta_e=-2*E_i+H*matrix[i][j]      
                if delta_e<0:
                    matrix[i][j]=-matrix[i][j]
                if delta_e>0:
                    p=random.uniform(0,1)
                    if p<math.exp(-delta_e/T):
                        matrix[i][j]=-matrix[i][j]
                flips+=1
            if method == 'Wolff':
                fli=Wolff_Algorithm_neighbor(matrix,T,H)
                flips+=fli
        x.append(H)
        y.append(M_tot(matrix)/(L**2))
        print(H)
        print(M_tot(matrix)/(L**2))
        print()
        print()
        H+=delta_H
       
    plt.plot(x, y,'ro',color="Blue")
    plt.xlabel('H')
    plt.ylabel('M')
    plt.title('H x M')
    plt.show()
   
def print_matrix_simbolically(matrix):
    n=len(matrix)
    m=len(matrix[0])
    for i in range(n):
        for j in range(m):
            if matrix[i][j]==-1:
                matrix[i][j]="G"
                print("%3s" %matrix[i][j],end=" ")
                matrix[i][j]=-1
            if matrix[i][j]==1:
                matrix[i][j]="."
                print("%3s" %matrix[i][j],end=" ")
                matrix[i][j]=1
        print()
    print()

def average_wolff(Ti,Tf,DT,IC,L,H,board,m,d,eq):
    T=Ti
    y=[]
    x=[]
    if IC==1:
        matrix=matrix_man1(L,L)
    if IC==2:
        matrix=matrix_rand(L,L)
    if IC==3:
        matrix=matrix_man2(L,L)
    eq=eq*L**2
    flips=0
    matrix_i=matrix
    Y=0
    while T<Tf or T==Tf:
        flips=0
        while flips<eq:
            i=random.randint(0,L-1)
            j=random.randint(0,L-1)
            E_i=energy_2(matrix,i,j,H)
            delta_e=-2*E_i+H*matrix[i][j]        
            if delta_e<0:
                matrix[i][j]=-matrix[i][j]
            if delta_e>0:
                p=random.uniform(0,1)
                if p<math.exp(-delta_e/T):
                    matrix[i][j]=-matrix[i][j]
            flips+=1
        for i in range (m):
            matrix=matrix_i
            U=Wolff_Algorithm_neighbor(matrix,T,H)
            Y+=U
        Y=Y/m
        Y=Y/(L**2)
        y.append(Y)
        x.append(T)
        print(Y)
        print(T)
        print()
        T+=DT
    end=time.time()
    print(end-start)
    plt.plot(x, y,'ro',color="Blue")
    plt.xlabel('Temperature')
    plt.ylabel('C/N')
    plt.title('Temperature x C/N')
    plt.show()
   
def thermalization(T,n,method,IC,H,L):
    if IC==1:
        matrix=matrix_man1(L,L)
    if IC==2:
        matrix=matrix_rand(L,L)
    if IC==3:
        matrix=matrix_man2(L,L)
    x=[]
    y=[]
    m=0
   
    for i in range (n):
        i=random.randint(0,L-1)
        j=random.randint(0,L-1)
        E_i=energy_2(matrix,i,j,H)
        delta_e=-2*E_i
        if delta_e<0:
            matrix[i][j]=-matrix[i][j]
        if delta_e>0:
            p=random.uniform(0,1)
            if p<math.exp(-delta_e/T):
                matrix[i][j]=-matrix[i][j]
        m+=1
        if m%10000==0:
            x.append(m)
            y.append(M_tot(matrix)/(L**2))
    plt.plot(x, y,'ro',color="Blue")
    plt.xlabel('Flips')
    plt.ylabel('M(Blue)')
    plt.title('Flips x M')
    plt.show()
print()
print()
print()
print()
print()
print()

process=int(input("Do you want to measure the thermalization, to flip, to measure the hysteresis,non-linear relaxation time,linear relaxation time or average Wolff(1,2,3,4,5,6)? "))
if process==1:
    n=int(input("Type the quantity of Monte Carlo steps: "))
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    method=input("Type the desired Monte Carlo method (Metropolis or Wolff): ")
    L=int(input("Type the length of the square lattice: "))
    d=int(input("Type the dimension of the model (2 or 8): "))
    T=float(input("Type the temperature: "))
    H=float(input("What is the external field? "))
    thermalization(T,n,method,IC,H,L)
   
if process==2:
    T_i=float(input("Type the initial temperature: "))
    T_f=float(input("Type the final temperature: "))
    delta_T=float(input("Type the difference between temperatures: "))
    m=int(input("Type the quantity of groups that will calculate the averages: "))
    steps=int(input("Do you want to change the Monte Carlo steps per site based on the temperature(1,2)? "))
    L=int(input("Type the length of the square lattice: "))
    d=int(input("Type the dimension of the model (2 or 8): "))
    method=input("Type the desired Monte Carlo method (Metropolis or Wolff): ")
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    CORR=input("Do you want to calculate the truncated correlation? It will significantly increase the time run. ")
    board=input("Do you want a fixed board? ")
    H=float(input("What is the external field? "))
    po=int(input("Type the Monte Carlo steps per site that will be used for the thermalization: "))
    flip(T_i,T_f,delta_T,m,steps,L,d,method,IC,CORR,H,board,po)
   
   
if process==3:
    H_0=float(input("Type the initial magnetic field: "))
    H_f=float(input("Type the final magnetic field: "))
    delta_H=float(input("Type the difference between the magnetic fields: "))
    n=int(input("Type the quantity of Monte Carlo steps: "))
    T=float(input("Type the temperature: "))
    method=input("Type the desired Monte Carlo method (Metropolis or Wolff): ")
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    L=int(input("Type the length of the square lattice: "))
    hysteresis(H_0,H_f,delta_H,n,T,method,IC,L)

if process==4:
    T=float(input("Type the temperature: "))
    method=input("Type the desired Monte Carlo method (Metropolis or Wolff): ")
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    L=int(input("Type the length of the square lattice: "))
    steps=int(input("Type the quantity of Monte Carlo steps: "))
    H=float(input("What is the external field? "))
    board=input("Do you want a fixed board? ")
    m=int(input("Type the quantity of groups that will calculate the averages: "))
    d=int(input("Type the dimension of the model (2 or 8): "))
    relaxation_function_noneq(T,method,IC,L,steps,H,board,m,d)

if process==5:
    T=float(input("Type the temperature: "))
    method=input("Type the desired Monte Carlo method (Metropolis or Wolff): ")
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    L=int(input("Type the length of the square lattice: "))
    steps=int(input("Type the quantity of Monte Carlo steps: "))
    H=float(input("What is the external field? "))
    board=input("Do you want a fixed board? ")
    m=int(input("Type the quantity of groups that will calculate the averages: "))
    d=int(input("Type the dimension of the model (2 or 8): "))
    po=int(input("Type the Monte Carlo steps per site that will be used for the thermalization: "))
    relaxation_function_eq(T,method,IC,L,steps,H,board,m,d,po)
   
if process==6:
    Ti=float(input("Type the initial temperature: "))
    Tf=float(input("Type the final temperature: "))
    DT=float(input("Type thh difference between temperatures: "))
    IC=int(input("Choose the initial condition(1,2 or 3): "))
    L=int(input("Type the length of the square lattice: "))
    H=float(input("What is the external field? "))
    board=input("Do you want a fixed board? ")
    m=int(input("Type the quantity of groups that will calculate the averages: "))
    d=int(input("Type the dimension of the model (2 or 8): "))
    eq=int(input("Type the quantity of steps to equilibrate: "))
    average_wolff(Ti,Tf,DT,IC,L,H,board,m,d,eq)
