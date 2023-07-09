
print("\n\n ==================================== ASSIGNMENT2 =======================================================")
print("Name: SAURABH KUAMR")
print("ROLL: 2020541")
print("\nGiven Sequence:")
class A_ACID: #class for amino acid and  respective p_alpha and p_beta value
    def __init__(self,amino_acid,pa,pb):
        self.aa = amino_acid
        self.pa = pa
        self.pb = pb

# the sequence of the amino acid
seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"


# list to store the amino acid and corresponding Chou and Fasman parameters.
list = []

list.append(A_ACID("E",1.53,0.26))
list.append(A_ACID("A",1.45,0.97))
list.append(A_ACID("L",1.34,1.22))
list.append(A_ACID("H",1.24,0.71))
list.append(A_ACID("M",1.20,1.67))
list.append(A_ACID("Q",1.17,1.23))
list.append(A_ACID("W",1.14,1.19))
list.append(A_ACID("V",1.14,1.65))
list.append(A_ACID("F",1.12,1.28))
list.append(A_ACID("K",1.07,0.74))
list.append(A_ACID("I",1.00,1.60))
list.append(A_ACID("D",0.98,0.80))
list.append(A_ACID("T",0.82,1.20))
list.append(A_ACID("S",0.79,0.72))
list.append(A_ACID("R",0.79,0.90))
list.append(A_ACID("C",0.77,1.30))
list.append(A_ACID("N",0.73,0.65))
list.append(A_ACID("Y",0.61,1.29))
list.append(A_ACID("P",0.59,0.62))
list.append(A_ACID("G",0.53,0.81))




def p_alpha(aa):    # this function returns the P_alpha of the amino acid aa
    ans = 0;
    for i in list:
        if(i.aa == aa):
            ans = i.pa
    return ans
def p_grt_1(window):       # returns the no of amino acid in window having  P(H) > 1
    sm = 0
    for i in window:
        if(p_alpha(i) >= 1):
            sm = sm +1
    return sm


def p_beeta(aa):   # this function returns the p_beeta of the amino acid aa #
    ans = 0
    for i in list:
        if(i.aa == aa):
            ans = i.pb
    return ans

def b_grt_1(window):  # returns the no of amino acids in window having P(S) > 1
    sm = 0
    for i in window:
        if(p_beeta(i) >= 1):
            sm = sm + 1
    return sm

def sum_palpha_window(window): # returns the sum Pa/P(H) of the amino acids in the given window
    r = 0
    for i in window:
        r = r + p_alpha(i)
    return r

def sum_pbeeta_window(window):   #returns the total sum of P(S) of each residue of the window
    r = 0
    for i in window:
        r = r + p_beeta(i)
    return r

helix = []          # this will store the positions of the candidates for helix
for i in range(len(seq)):
    helix.append("*")

sheet = []          ## this will store the possible candidates for sheets
for i in range(len(seq)):
    sheet.append("*")    # "*" indicate empty position on sheets

def HELIX_extend_r(wndow,ptr):
    if(ptr < len(seq) -1):
        st = seq[ptr-2:ptr+2]   # extending to 1 right amino acid
        l = ptr-2
        j = ptr + 1
        while(sum_palpha_window(st) >= 4):
            helix[j] = "H"
            if(j+1 > len(seq)):
                break
            else:
                j = j + 1
                l = l + 1
                st = seq[l:j + 1]
def HELIX_extend_l(window , ptr):
    if(ptr - 1 >= 0):
        l = ptr-1
        st = seq[l:l+4]
        while(sum_palpha_window(st) >= 4):
            helix[l] = "H"
            if(l-1 >= 0):
                l = l-1
                st = seq[l:l+4]


def PREDICT_HELIX(indx): # the function will predict the helix region
    while (indx != len(seq)):
        window = seq[indx:indx + 6]    #taking window of length 6
        if (p_grt_1(window) >= 4): # means the window can be a helix

            for ind in range(indx, indx + 6):  # this window is in helix
                helix[ind] = "H"
            # left extend
            HELIX_extend_l(window, indx)
            # right extend
            HELIX_extend_r(window, indx + 5)

        indx = indx + 1


def BEETA_extend_l(window, ptr):#this fun will check for extending left during beeta sheet calculation
    if (ptr - 1 >= 0):
        l = ptr - 1
        st = seq[l:l + 4]
        while (sum_pbeeta_window(st) >= 4):
            sheet[l] = "S"
            if (l - 1 >= 0):
                l = l - 1
                st = seq[l:l + 4]


def BEETA_extend_r(window,ptr):   #this fun will check for extending right during beeta sheet calculation
    if (ptr < len(seq) - 1):
        st = seq[ptr - 2:ptr + 2]  # extracting new 4 size window with one extra extend to right
        l = ptr - 2
        j = ptr + 1
        while (sum_pbeeta_window(st) >=4):
            sheet[j] = "S"
            if (j + 1 > len(seq)):
                break
            else:
                j = j + 1
                l = l + 1
                st = seq[l:j + 1]


def PREDICT_SHEETS(indx):        #predicting the beeta sheets region
    while(indx != len(seq)):
        window = seq[indx:indx+5]   # taking window of length 5
        if(b_grt_1(window) >= 3):   #the window cna be sheets
            for ind in range(indx,indx+5):
                sheet[ind] = "S"
            #LEFT EXTEND
            BEETA_extend_l(window,indx)

            #RIGHT EXTEND
            BEETA_extend_r(window,indx+4)

        indx = indx+ 1


# identifying helix
PREDICT_HELIX(0)
PREDICT_SHEETS(0)

ANS = []
for i in range(len(seq)):
    ANS.append("*")      ##-----------------------------------------------------------------------------------------------------------
def NO_OVERLAP():
    for i in range(len(seq)):
        if(helix[i] == "*" and sheet[i] == "S"):
            ANS[i] = "S"
        if(helix[i] == "H" and sheet[i] == "*"):
            ANS[i] = "H"


def CONFLICT(ptr1, ptr2):             # this fn will resolve the conflicting region
    while(ptr1 < len(seq)-1):
        #print("Saurabh")
        if(helix[ptr1] == "H" and sheet[ptr1] == "S"):
            while(helix[ptr2+1] == "H" and sheet[ptr2] == "S"):
                ptr2 = ptr2 + 1
            window = seq[ptr1: ptr2+1]
            if (sum_palpha_window(window) >= sum_pbeeta_window(window)):
                for j in range(ptr1, ptr2 + 1):
                    ANS[j] = "H"
            if (sum_palpha_window(window) < sum_pbeeta_window(window)):
                for j in range(ptr1, ptr2 + 1):
                    ANS[j] = "S"

            ptr1 = ptr2
            #print("hi")
        ptr1 = ptr1+1
        ptr2 = ptr2+1



print(seq)
print()
for i in helix:
    print(i,end ="")

# idendifying sheets
print()

for i in sheet:
    print(i,end = "")

NO_OVERLAP()
CONFLICT(0,0)
print()

print("\nRESULT:")
for i in ANS:
    print(i, end = "")


