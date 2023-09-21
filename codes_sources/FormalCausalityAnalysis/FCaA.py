# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 16:06:41 2020

@author: abazin
"""
	
# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)

import subprocess
from subprocess import PIPE
import numpy as np
from oct2py import Oct2Py
from oct2py import octave
import copy
from sklearn import ensemble
import time
import random
if __name__ == "__main__":
    from PCA import properPremises
else:
    from .PCA import properPremises

import re
from sklearn.feature_selection import mutual_info_classif,SelectKBest,chi2,f_regression
from itertools import islice
from random import randint


'''
with open("randomData1000-10.txt") as file: # pareil avec un .txt plein lignes pleines de nombres
    
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = np.array([x.split(",") for x in data]).astype(float)
'''


'''
if __name__ == "__main__":
    with open("rapsodyn_data_tot_new.csv") as file:
        
        data = file.readlines()
        data = [x.replace("\n","") for x in data]
        data = [x.replace(",",".") for x in data]

        # lister les var. 1 à 14 de l'en-tête
        variable_names = np.array([x.split(";") for x in data])[0,4:].astype(str)

        # générer la matrice des valeurs, en float
        data = np.array([x.split(";") for x in data])[1:,4:]#.astype(float)

        # remplacer les cases vides de "11 - seed_nb" et "12 - PMG" par la moyenne
        for X in range(len(data)):
            for Y in range(len(data[0])):
                if data[X,Y] == '':
                    data[X,Y] = '0'
        data = data.astype(float)
        S = 0
        n = 0
        for X in data[:,11]:
            if X != 0:
                S += X
                n += 1
        for x in range(len(data[:,11])):
            if data[x,11] == 0:
                data[x,11] = S/n
        S = 0
        n = 0
        for X in data[:,12]:
            if X != 0:
                S += X
                n += 1
        for x in range(len(data[:,12])):
            if data[x,12] == 0:
                data[x,12] = S/n
    
    
    def loadTestGEENAGE():
        # transformer le .csv transposé en matrice numpy, sans les 1ère colonne
        data = list()
        file = open("STAN_RNASeq_NormalizedCounts_65STADESeq2csv2_corrected.csv","r")
        for line in file:
            row = line.split("\n")[0]
            row = row.split(";")
            data.append(list(row))
        data = np.array(data).transpose()
        data = np.delete(data,0,0)
        
        
        # récuperer les (nom, condition) de chaque entrée
        classes = list()
        file2 = open("STAN_RNAseq_MetaData_65STA_fullLORIA_Alex.csv","r")
        for line in file2:
            row = line.split("\n")[0]
            row = np.array(row.split(";"))[[0,4]]
            classes.append(list(row))
        
        # dans data, remplacer le nom par 0 si sain, 1 si cardiaque
        classes = np.array(classes)
        data[0,0] = 'Class'
        for j in range(data.shape[0]):
            if data[j,0] == '"2H18"':
                    data[j,0] = 0
            else:
                for i in range(classes.shape[0]):
                    if data[j,0] == classes[i,0]:
                        if classes[i,1] == '"cardiac"':
                            data[j,0] = 1
                        else:
                            data[j,0] = 0 
                        
        print(data.shape)
                        
        # la matrice + la liste des 0 1
        data_raw = data[1:,1:].astype(float)
        data_classes = data[1:,0].astype(float)
        
        
        #data_raw_new = SelectKBest(chi2, k = 1000).fit_transform(data_raw,data_classes)
        
        selector = SelectKBest(chi2, k = 100)
        selector.fit(data_raw,data_classes)
        # récupère les indices des 100 features (=parmi les 58k gènes de NormalizedCounts) les plus "significatives"
        cols = selector.get_support(indices=True)
        cols2 = [0]+[x+1 for x in cols]
        data = data[:,cols2]
        
        # renvoie un tableau de 66 lignes-patients, et 100-gènes-les-plus-significatifs
        return data    
'''

'''
with open("pair2.txt") as file:
    
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = np.array([x.split(";") for x in data]).astype(float)
'''

'''
with open("C:\\Users\\abazin\\Documents\\Code\\Data\\sobar-72.csv") as file: # attention on l'a pas lui
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = np.array([x.split(",") for x in data])[1:,:].astype(float)
'''

'''
with open("C:\\Users\\abazin\\Documents\\Code\\Data\\avila-ts.txt") as file:
    # recopier le fichier, mais en remplaçant les lettres dans la dernière colonne par leur ordre dans l'alphabet
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = np.array([x.split(",") for x in data])
    data[:,10] = [ord(A)-65 for A in data[:,10]] 
    data = data.astype(float)
f = open("avila-ts.txt","w")
for L in data:
    for l in range(len(L)-1):
        f.write(str(L[l]))
        f.write(",")
    f.write(str(L[len(L)-1]))
    f.write("\n")
f.close()
'''



'''
with open("C:\\Users\\abazin\\Documents\\Code\\Data\\iris.data") as file:
    # recopier le fichier, mais en remplaçant les noms scientifiques par 0, 1, 2
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = [x.replace("Iris-setosa","0") for x in data]
    data = [x.replace("Iris-versicolor","1") for x in data]
    data = [x.replace("Iris-virginica","2") for x in data]
    data = [x.split(",") for x in data]
    data = np.array(data).astype(float)
f = open("iris.txt","w")
for L in data:
    for l in range(len(L)-1):
        f.write(str(L[l]))
        f.write(",")
    f.write(str(L[len(L)-1]))
    f.write("\n")
f.close()
'''

'''
with open("diabetes.csv") as file:
    # recopier le fichier, mais en enlevant l'en-tête + convertir en float
    data = file.readlines()
    data = [x.replace("\n","") for x in data]
    data = np.array([x.split(",") for x in data])[1:,:].astype(float)
    
with open("diabetes.txt","w") as file:
    for i in data:
        for j in range(len(i)):
            file.write(str(i[j]))
            if j < len(i)-1:
                file.write(";")
        file.write("\n")
'''
    
'''
with open("lsac.data") as file:
        D = file.readlines()
# numériser les données catégrorielles
D = [x.replace("\n","") for x in D]
D = [x.replace("Yes","1") for x in D]
D = [x.replace("No","0") for x in D]
D = [x.replace("Female","1") for x in D]
D = [x.replace("Male","0") for x in D]
D = [x.replace("White","2") for x in D]
D = [x.replace("Black","1") for x in D]
D = [x.replace("Other","0") for x in D]
D = [x.replace("Passed","1") for x in D]
D = [x.replace("Failed_or_not_attempted","1") for x in D]

data = np.array([np.array(x.split(";")) for x in D])[1:,:].astype(float)
list_features = np.array([np.array(x.split(";")) for x in D])[0,:]
'''

'''
with open("C:\\Users\\abazin\\Documents\\Code\\Causality\\wdbc.data") as file:
        D = file.readlines()

D = [x.replace("\n","") for x in D]
 
D = [x.replace("M","1") for x in D]
D = [x.replace("B","0") for x in D]
        
data = np.array([np.array(x.split(",")).astype(float) for x in D])
'''

'''
with open("pair0071.txt") as file:
        D = file.readlines()  
D = [x.replace("\n","") for x in D]
data = np.array([np.array(re.split("\s",x)) for x in D])   
data = np.array([X[:len(X)-1].astype(float) for X in data]) 
'''

'''
with open("Tic2000/ticdata2000.txt") as file:
        D = file.readlines()  
D = [x.replace("\n","") for x in D]
data = np.array([np.array(re.split("\s",x)) for x in D])   
data = np.array([X[:len(X)-1].astype(float) for X in data]) 
''' 



   
oc = Oct2Py()


#Rules = []
#Concepts = []
RulesRandom = []



def testCaus(X,Y):

    # écrit un fichier dont chaque ligne est "les valeurs X puis les valeurs Y"
    with open("inputErgo.txt","w")as fin:
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                fin.write(str(X[i][j]))
                fin.write(";")
            for j in range(Y.shape[1]):
                fin.write(str(Y[i][j]))
                if j < Y.shape[1]-1:
                    fin.write(";")
            fin.write("\n")
            
    args = ['java','-jar','ergo.jar',
            '-FILE_INPUT','inputErgo.txt',
            '-FILE_RUNTIME_OUTPUT','outputErgo.txt',
            '-NUM_ROWS',str(X.shape[0]),
            '-NUM_MEASURE_COLS',str(X.shape[1]+Y.shape[1]),
            '-ALPHA','0.3',
            '-CLUMPS','10',
            '-MAXDIMX',str(X.shape[1]-1)]
    # args = ['java','-cp','C:\\Users\\abazin\\Documents\\Code\\FormalCausalityAnalysis', 'Pairs...

    p = subprocess.Popen(args, stdin=PIPE, stdout=PIPE, text=True)
    #p = subprocess.run(args,capture_output=True)
    #p = subprocess.check_call("java -class Pairs.class -FILE_INPUT pair2.txt -FILE_RUNTIME_OUTPUT rt2.txt -NUM_ROWS 120 -NUM_MEASURE_COLS 8 -ALPHA 0.3 -CLUMPS 10 -MAXDIMX 5")
    p.communicate()

    with open("outputErgo.txt") as fout:
        D = fout.readlines()
        if len(D) > 0:
            R = D[0].replace("\n","")
        else:
            R = ""
    # renvoie TRUE si X infère Y, FALSE sinon
    return R == "->"



def testCorr(file,indX,indY):
    # renvoie le score LCCA entre les colonnes d'indX et d'indY
    
    indX = [int(x+1) for x in indX]
    indY = [int(y+1) for y in indY]
    
    #print("indX = ", indX)
    #if len(indX) > 0:
    #    print("type indX = ",type(indX[0]))
    #print("indY = ", indY)
    
    return oc.LCCA(file, list(indX), list(indY))



def testImpli(data, indX, indY, file):
    # renvoie si les colonnes X ""causent"" les colonnes Y selon nos critères
    
    return testCorr(file, indX, indY) > 0.7 and testCaus(data[:,indX],data[:,indY])



def causClosure(data, X, file):
    # construit la clôture causale de X
    
    R = copy.deepcopy(X)
    R2 = copy.deepcopy(R)
    fin = False
    while not fin:
        for i in range(data.shape[1]):
            if not i in X:
                if testImpli(data, X, [i], file):
                    R += [i]
        if len(R) == len(R2):
            fin = True
        R2 = R
                
    return np.array(R)



def oplusCausality(A, a, data, file):
    # garde dans A tous les élements plus petits que a
    B = set(copy.deepcopy(A))
    for x in A:
        if x > a:
            B.remove(x)
    B.add(a)
    # en renvoie la clôture d'implication (à partir de la liste Rules)
    C = logicalClosure(B)
    return list(C)



def NextCausality(A, data, file):
    for i in reversed(range(data.shape[1])):
        # les elts < i sont ceux "pas encore visités"
        if not i in A:
            B = oplusCausality(A, i, data, file)
            fin = True
            for j in B:
                if j < i and j not in A:
                    fin = False
            if fin:
                return B


#Computes the concepts of a bidimensional context
def NextClosureCausality(data,file,v=0):
    start_time = time.time()
    global Rules
    global Concepts
    Rules = []
    Concepts = []
    A = causClosure(data, [], file) # (tautologies ?)
    while len(A) < data.shape[1]:
        print(A)
        B = causClosure(data, list(A), file)
        if len(B) > len(A):
            Rules += [[A,list(set(B).difference(set(A)))]] # "A implique les autres var de B"
        else:
            Concepts += [A]
        A = NextCausality(A, data, file)
    Concepts += [A]
    
    # on écrit dans un fichier toutes les implications trouvées
#    Q,E = computeFunctionRules()
    f = open("Expe_agro_Zac_v%i.txt" % v,"w")
    print(str(time.time()-start_time)+" seconds")
    f.write(str(time.time()-start_time)+" seconds\n")
    for R in range(len(Rules)):
        for p in Rules[R][0]:
            f.write(str(p)+" ")
        f.write("-> ")
        for c in Rules[R][1]:
            f.write(str(c)+" ")
        #f.write("("+str(E[R])+")")
        f.write("\n")
    f.close()


#Logical / transitive closure of a set by a set of implication rules
# /!\ fonctionne grâce au Rules rempli dans l'algo précédent
def logicalClosure(X):
    global Rules
    S = copy.copy(set(X))
    fin = False    
    while not fin:
        s = len(S)
        for I in Rules:
            if set(I[0]).issubset(S):
                S = S.union(I[1])
        if len(S) == s:
            fin = True         
    return list(S)


def oplusLogical(A, a):
    B = set(copy.deepcopy(A))
    for x in A:
        if x > a:
            B.remove(x)
    B.add(a)
    C = logicalClosure(list(B))
    return list(C)



def NextLogical(A, data):
    for i in reversed(range(data.shape[1])):
        if not i in A:
            B = oplusLogical(A, i)
            fin = True
            for j in B:
                if j < i and j not in A:
                    fin = False
            if fin:
                return B


#Computes the concepts of a bidimensional context
def NextClosureLogical(data):
   
    Closed = []
    A = logicalClosure([])
    while len(A) < data.shape[1]:
        Closed += [A]
        A = NextLogical(A, data)
    Closed += [A]
    return Closed


def writeRules():
    for r in Rules:
        print(r[0],"->",r[1])



def Wild(data): # Geenage?
    global Rules
    
    maxe = []
    # maxe[e] initialisé à la liste [V\{e}]
    for a in range(data.shape[1]):
        maxe += [[set(range(data.shape[1])).difference(set([a]))]] 
    
    for i in range(len(Rules)):
        for e in range(data.shape[1]):
            maxe2 = []
            for X in maxe[e]: # X = un ensemble dans maxe[e]
                if not set(Rules[i][0]).issubset(set(X)) or set(Rules[i][1]).issubset(X):
                    # X pas prémisse, ou X conclusion
                    maxe2 += [X] # on ajoute X
            maxet = []
            for Z in maxe[e]:
                if not Z in maxe2:
                    maxet += [Z]
            maxe3 = maxet # Z est prémisse et n'est pas conclusion

            T = copy.copy(maxe2)
            for a in range(data.shape[1]):
                # on ajoute à T les intersections entre maxe3 et maxe[a]
                if a != e:
                    for X,Y in zip(maxe3,maxe[a]):
                        Z = X.intersection(Y)
                        if not Z in T:
                            T += [Z]
            maxe4 = []
            for A in T: # on garde que les éléments maximaux (pas sous-ensemble d'un autre de T)
                maximum = True
                for B in T:
                    if A!= B and set(A).issubset(set(B)):
                        maximum = False
                        break
                if maximum:
                    maxe4 += [A]
            maxe[e] = maxe4
    
    R = []
    for e in range(data.shape[1]):
        for X in maxe[e]:
            if not X in R:
                R += [X]
    
    return R # les ensembles des maxe sans doublon


def printContext(C,f):
    with open(f,"w") as file:
        for c in C:
            for e in range(len(c)):
                file.write(str(list(c)[e]+1))
                if e < len(c)-1:
                    file.write(",")
            file.write("\n")
            
            
            
def AvCorrCoef(labels,preds):
    R = 0
    
    labels = np.array(labels)
    preds = np.array(preds)
    
    for i in range(labels.shape[1]):
        
        avLabels = np.mean(labels[:,i])
        avPreds = np.mean(preds[:,i])
        
        num = 0        
        denum1 = 0
        denum2 = 0
        for j in range(labels.shape[0]):
            num += (labels[j,i]-avLabels)*(preds[j,i]-avPreds)
            denum1 += (labels[j,i]-avLabels)**2
            denum2 += (preds[j,i]-avPreds)**2
            
            denum = np.sqrt(denum1*denum2)

        if denum != 0:
            R += num/denum
        else:
            R = 0

    return R/labels.shape[1]            
            
            
            
            
def trainFunctionImplication(X,Y):
    C = ensemble.RandomForestRegressor()
    C.fit(X,Y)
    Pred = C.predict(X)
    if len(Y.shape) == 1:
        Y = [[x] for x in Y]
        Pred = [[x] for x in Pred]
    E = AvCorrCoef(Y,Pred)
    return C,E
    
    
def computeFunctionRules(data): # Geenage?
    Cl = []
    Ev = []
    for R in Rules:
        Y = data[:,R[1]]
        if Y.shape[1] == 1:
            Y = np.ravel(data[:,R[1]])
        C,E = trainFunctionImplication(data[:,R[0]],Y)
        Cl += [C]
        Ev += [E]
        print(R[0],"->",R[1],"(",E,")") # print la règle et son score de corrélation
    return Cl, Ev

if __name__ == "__main__":
    
    import ast
    from sys import argv
    
    with open(argv[1], "r") as file: 
        # argv[1] : "matMFCaA.txt" ou "matMFCaA2.txt"
        data = file.read()
    
    M=data.replace('\n','')
    M=ast.literal_eval(M)
    print(M)
    
    NextClosureCausality(np.array(M),argv[2], argv[3])
    # argv[2] : "file25_30.txt" ou "file25_302.txt"
    # argv[3] : "1" ou "2"
    
    
    # M=[[1,1,1,1,1,1,1,1,1]]
    # print(M)
    # NextClosureCausality(np.array(M),"fileExpe.txt")



def generateEffect(Variables):
    #combi = random.randint(1,3)
    combi = 1 # on additionne, multiplie, ou moyenne les variables de chaque ligne ?
    #degrePoly = int(np.floor(random.random()*3+1))
    degrePoly = 1
    coef = [random.random()*3 for i in range(degrePoly)]
    R = []
    stringpoly = ""
    for i in range(Variables.shape[0]):
        S = random.gauss(0,1) #the noise
        if combi == 1:
            combiS = 0
            for j in range(Variables.shape[1]):
                combiS += Variables[i,j]
        if combi == 2:
            combiS = 1
            for j in range(Variables.shape[1]):
                combiS *= Variables[i,j]   
        if combi == 3:
            combiS = 0
            for j in range(Variables.shape[1]):
                combiS += Variables[i,j]
            combiS = combiS/Variables.shape[1]

        for c in range(len(coef)):
            # génère l'image de combiS par un polynôme de deg° "degrePoly" et de coefs "coef"
            CombiS2 = 1
            for i in range(c+1):
                CombiS2 = CombiS2*combiS
            #S += coef[c]*(combiS^(c+1))
            S += coef[c]*(CombiS2)

        R += [S] # pour chaque ligne voilà
    
    for c in range(len(coef)):
        stringpoly += str(coef[c])+"x^"+str(c+1)
        if c < len(coef)-1:
            stringpoly += "+"
    if combi == 1:
            stringpoly += "   x=addition"
    if combi == 2:
            stringpoly += "   x=multiplication"
    if combi == 3:
            stringpoly += "   x=average"
        
    return R,stringpoly # on renvoie le __repr__ du polynôme et les images des variables



def random_chunk(li, min_chunk=1, max_chunk=3):
    # génère successivement li par bouts de 1 à 3 elts
    it = iter(li)
    while True:
        nxt = list(islice(it,randint(min_chunk,max_chunk)))
        if nxt:
            yield nxt
        else:
            break



def generateData(nbObjects, nbVariables):
    V = list(range(nbVariables))
    # on découpe [0,...,nbVariables-1] en des morceaux pas plus grands que 25% de sa taille
    V = list(random_chunk(V,1,np.floor(nbVariables/4)))
    Layers = list(random_chunk(V,1,3)) # on découpe encore en une liste de liste de liste
    Edges = []
    Explication = []
    for i in range(len(Layers)-1):
        for j in Layers[i]:
            for k in Layers[i+1]:
                if random.random() < 0.5:
                    Edges += [[j,k]]
                    # la liste des (j->k), avec j k deux listes de Layers adjacentes
    
    Data = np.zeros((nbObjects,nbVariables))
    # on initialise le contexte
    for L in Layers:
        for v in L:
            P = [] # liste des points allant vers la liste v
            for E in Edges:
                if E[1] == v:
                    P += [E[0]]
            if P == []: # si aucun effet ne pointe vers cet ensemble
                for var in v:
                    mu = random.random()*100-50
                    sigma = random.random()*5
                    for x in range(nbObjects):
                        Data[x,var] = random.gauss(mu, sigma) # on assigne des valeurs au pif partout
            else:
                allEffects = []
                for var in v:
                    for p in P:
                        Eff,Poly = generateEffect(Data[:,p])
                        allEffects += [Eff]
                        Explication += [p,var,Poly] # x, y, poly qui lie x à y
                    allEffects = np.array(allEffects)
                    for x in range(nbObjects):
                        # moyenne de toutes les images des causes
                        Data[x,var] = np.mean(allEffects[:,x])
    # valeurs, groupements d'attributs, arêtes du graphe causal, polynôme de chaque arête
    return Data, V, Edges, Explication


def matrix2File(Matrix,file): # assez explicite
    f = open(file,"w")
    for X in Matrix:
        for x in range(len(X)-1):
            f.write(str(X[x]))
            f.write(",")
        f.write(str(X[len(X)-1]))
        f.write("\n")
    f.close()
    
    
def expli2File(Expli,file):
    f = open(file,"w")
    for X in Expli:
        f.write(str(X))
        f.write("\n")
    f.close()
    
 

def genContext(n,m,p): # liste d'arêtes aléatoires entre [0:n-1] et [0:m-1]
    C = []
    for i in range(n):
        for j in range(m):
            if random.random() < p:
                C += [[i,j]]
    return C
        
        
def generateRules(n,m,p):
    C = genContext(n,m,p)
    return properPremises((C,n,m))
    # X prémisse propre de a, ssi X cause {a} mais qu'aucun de ses sous-ensembles ne cause {a}
    # en l'occurrence de X'', càd les attributs possédés par tous les objets possédant tous les X
    # bref une implication assez triviale





def minCauses(data, file, target):
    
    f = open("results_GEENAGE_0.9.txt","w")
    M = []
    T = set(range(data.shape[1])).difference(set([target]))
    for i in T:
        M += [T.difference(set([i]))] # les T\{i}, pour T les attributs hors target
    
    fin = False
    
    R = []
    I = []
    
    
    while not fin:
        
        print("M est de taille",len(M))
        
        
        finS = True
        m = T
        M2 = copy.copy(M)
        adj = 0
        for x in range(len(M2)):
            cx = testImpli(data, list(M2[x]), [target], file)
            if cx: # si un M2[x] cause target, stop : on va au for suivant
                finS = False
                m = set(M2[x])
                break
            else:
                M.pop(x-adj) # si pas de cause, on retire ledit M2[x] et on reprend
                x = x-1
                adj = adj+1
                
        while not finS:
            finS = True
            for i in m: # on enlève des elts de m jusqu'à ce qu'il soit maximal
                ci = testImpli(data, list(m.difference(set([i]))), [target], file)
                if ci:
                    finS = False
                    m = m.difference(set([i]))
                    break
                
        if len(m) != len(T):
            # on n'est pas passé que par des else, donc par un "if cx"->"while not finS"
            # le m est une prémisse précise de target
            R += [m]
            I += [[m,T]]
            
            f.write(str(m)) # on note la prémisse dans le fichier
            f.write("\n")
            
            print(len(R),"ème Prémisse")
            
            M2 = []
            for x in M:
                if not x.issuperset(m):
                    M2 += [x] # M2 contient M sauf ceux dans lesquels tout m est inclus
                else:
                    for a in m:
                        # M2 contient aussi les sous-ensembles de x
                        x2 = copy.copy(x)
                        x2.remove(a)
                        ajout = True
                        for m2 in M2:
                            if m2 == x2:
                                ajout = False
                        if ajout:
                            M2 += [x2]
                        # bref l'algo parcourt récursivement tous les sous-ensembles pertinents
            M = M2
        else:
            fin = True
            
    f.close()
        
    return R


def data2File(data): # Geenage?
    f = open("Geenage.txt","w")
    for L in data:
        for l in range(len(L)-1):
            f.write(str(L[l]))
            f.write(" ")
        f.write(str(L[len(L)-1]))
        f.write("\n")
    f.close()
    
        
        
def randomData(nbObject,nbAttribute,density):
    # crée un contexte formel binaire aléatoire, de densité donnée
    
    data = []
    for i in range(nbObject):
        L =  []
        for j in range(nbAttribute):
            p = random.random()
            if p < density:
                L += [1]
            else:
                L += [0]
        data += [L]
    return data    
        
        
        
def expesFCA4AI():
    Results = []
    for density in range(1,10): # /!\ pourquoi pas entre 0 et 1 ?
        O = []
        print("density",density)
        for nbObjects in range(1,100,2):
            print("nbobj",nbObjects)
            A = []
            for nbAttributes in range(1,100,2):
                print("nbatt",nbAttributes)
                data = randomData(nbObjects,nbAttributes,density)
                matrix2File(data,"fileExpe.txt")
                start_time = time.time()
                NextClosureCausality(np.array(data),"fileExpe.txt")
                end_time = time.time()
                A += [end_time-start_time]
            O += [A]
        matrix2File(O,"Results_FCA4AI_",nbObjects,"_",nbAttributes,"_",density,".txt")
        Results += [O]
    return Results
                




