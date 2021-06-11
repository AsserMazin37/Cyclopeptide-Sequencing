#1-Weight Map
file = open("weightToaa.txt")
weightDict = {}
for line in file:
    line = line.rstrip()
    word = line.split(':')
    key = int(word[0])
    value = word[1]
    weightDict[key] = value 

#Helper Functions
def get_key(val):
    for key, value in weightDict.items():
         if val == value:
             return key

def count(l):
    dic = {}
    for i in l:
        c=0
        for j in l:
            if i==j:
                c+=1
        dic[i] = c
    return dic

def check_freuency(seq , InitialDic):
    seqDict = {}
    for i in seq:
        c=0
        for j in seq:
            if i==j:
                c+=1
        seqDict[i] = c
        
    for key, value in seqDict.items():
         if InitialDic[key] < value:
             return False
    return True
         
#2-Linear spectrum Function
def linear_spectrum(peptide):
   subPeptideList = ['']
   for i in range(len(peptide)):
     for j in range(i , len(peptide)):
          subPeptideList.append(peptide[i:j+1])
   # Getting sum of spectrum made aa       
   subPeptideWeights = []
   for i in subPeptideList:
       sumWeights =0
       for j in i:
           sumWeights += get_key(j)
       subPeptideWeights.append(sumWeights)
   subPeptideWeights.sort()
   return subPeptideWeights

#3- Is Consistant Function
def isConsistent(subpeptide, spectrum , spectrumDic): # Spectrum is the input spectrum
    linearSpectra = linear_spectrum(subpeptide)
    for p in linearSpectra:
        if p not in spectrum:
            return False
    
    spectraDict = count(linearSpectra)
    # Check for the applicability of the frequencies with the original spectrum
    for key, value in spectraDict.items():
        # Same number occured more times than it was
        if spectrumDic[key] < value:
            # Originally sent
            return False
    # Applicable frequences of spectra with the original
    return True

#4-Initial List Funcion
def Initial_List(spectrum):
    initialList = []
    for i in spectrum:
        # Getting all 1 mers
        if i in weightDict.keys():
          initialList.append(weightDict[i])
    return initialList

def extend(Initial_L , TempList , InitialDic):
    mersList = []
    for i in TempList:
        for j in Initial_L:
            mer = i+j
            applicable = check_freuency(mer , InitialDic) 
            # Making extends that are applicable with the original chars
            if(applicable):
              mersList.append(mer)
    return mersList

# INPUT: 0 97 97 99 101 103 196 198 198 200 202 295 297 299 299 301 394 396 398 400 400 497
if __name__ == "__main__":
    theoriticalSpectrum = list(input("please input the spectrum list: ").split(" "))
    
    # Convert string to int using list comprehension
    theoriticalSpectrum = [int(i) for i in theoriticalSpectrum]
    Initial_L = Initial_List(theoriticalSpectrum)
    
    spectrumDic = count(theoriticalSpectrum)
    
    # Count frequency of every char
    InitialDic = count(Initial_L)
    
    TempList = Initial_L.copy()
    
    # Stop looping when tempList is empty as it carries all the consistant meres
    c = 1
    while c < len(Initial_L): 
     TempList = extend(Initial_L , TempList , InitialDic)
     consistantList = []
     for i in TempList:
        if isConsistent(i, theoriticalSpectrum , spectrumDic):
          consistantList.append(i)     
     # Stop when all are inconsistant
     if len(consistantList) == 0: 
        break
     TempList = set(consistantList.copy())
     c += 1
     
    finalList = TempList
    print(finalList)

