
#################### DATA ##################################################################
Mol_wt = {'G' : 57, 'A': 71, 'S': 87,
          'P': 97, 'V' : 99, 'T' : 101,
          'C' : 103 , 'I' : 113, 'L' : 113,
          'N' : 114, 'D' : 115, 'K' : 128, 
          'Q' : 128, 'E' : 129, 'M' : 131,
          'H' : 137, 'F' : 147, 'R' : 156,
          'Y' : 163, 'W' : 186 }

Reverse_MW = {57 : 'G', 71 : 'A',  87 : 'S',
           97 : 'P',  99 : 'V',  101 : 'T',
          103 : 'C'  , 113 : 'I', 113 : 'L',
          114 : 'N', 115 : 'D', 128 : 'K', 
           128 : 'Q', 129 : 'E',  131 : 'M',
           137 : 'H', 147 : 'F', 156 :  'R',
           163 : 'Y', 186 : 'W'}

aa_seq = 'APVE'
aa_size = 4
spec_meter = []
linear_spec = []
Circular_Weight = []
Linear_Weight = []
cal_wt = 0
initial_list = []

input_spec = [0, 97,97,99, 101,103, 196,198, 198, 200, 202, 295, 297,297,299,299,301,394,398,400,400,497]
Tmplist= []
Expand_list= []
peptide_length = 0
#############################################################################################

def Print_Spectrum(list1, list2):
      
    print ("Subpeptide", "  ", "Mass" )  
    for i in range (0, len(spec_meter),1):
        print (spec_meter[i], "           ",Circular_Weight[i])
        print()
    
    print(0)

    return 0

def Circular_Spectrum (seq):
############ Get all mers ############################################
    for i in range (1, len(seq), 1):
    
        for j in range (0, len(seq), 1):
           ite = (len(seq)-j)    
           if(j != len(seq)-1 ):
           
               if(len(seq[j:j+i]) == i):
                    tmp =  seq[j : j+i]
                    spec_meter.append(tmp)

               elif( (len(seq[j:j+1]) != i)):
                   tmp = seq[j:] + seq[:i-ite]
                   spec_meter.append(tmp)
                   ite+=1

           elif(j == len(seq)-1 ):           
               tmp = seq[j:] + seq[:i-1]
               spec_meter.append(tmp)
              
    spec_meter.append(seq)

#################### Calculate Circular Spectrum  #######################

    #for w in range (0, aa_size, 1):
    #    cal_wt = Mol_wt[aa_seq[w]]
    #    weight.append(cal_wt)

    for i in range (0, len(spec_meter), 1):
        tmp = spec_meter[i]
        cal_wt = 0 
        for j in range (0, len(tmp), 1):
        
            cal_wt += Mol_wt[tmp[j]]    
           
        Circular_Weight.append(cal_wt)
       
    Print_Spectrum(spec_meter, Circular_Weight)
    return 0

#print(Circular_Spectrum(aa_seq))

####################### Part 2 ########################################################
def Linear_Spectrum(seq):
    linear_spec=[]
    for i in range (1, len(seq), 1):
        
        for j in range (0, len(seq), 1):
            if(len(aa_seq[j:j+i]) == i and j+i-1 < len(seq)):
                tmp =  seq[j : j+i]
                linear_spec.append(tmp)
    
    linear_spec.append(seq)
    Linear_Weight=[]
########### Calculate Linear Spectrum ########################         
    for i in range (0, len(linear_spec), 1):
        tmp = linear_spec[i]
        cal_wt = 0 
        for j in range (0, len(tmp), 1):
            
            cal_wt += Mol_wt[tmp[j]]
            
        Linear_Weight.append(cal_wt)
    return Linear_Weight

#############################################################

def Initial_List(list):

    length = 0
    for i in range (0 , len(list), 1):
        tmp_var = list[i]
        if(tmp_var in Reverse_MW and Reverse_MW[tmp_var] not in initial_list ):

            initial_list.append(Reverse_MW[tmp_var])
            length +=1 

        elif(tmp_var in Reverse_MW ):
            length +=1

    
    return length , initial_list

#Initial_List(input_spec)
#######################################################################################


def Repetition(in_list):
    unique = set(in_list)
    if len(unique) == len(in_list):
        return True
    else:
        return False

def IsConsistent (seq, list_in):

    tmp_list = Linear_Spectrum(seq)
    count = 0
    done = [0 for i in range(10000)]
    
    for j in range(0, len(tmp_list),1):
            for i in range (len(input_spec)):
                if tmp_list[j]==input_spec[i] and done[i]==0:
                    done[i]=1
                    count +=1
                    break;
    if (count == len(tmp_list)):
        return True
    else:
        return False


def expand(list):

    expanded = []
    for k in range (0, len(list), 1):
        for key in Mol_wt.keys():
             expanded.append( list[k] + key ) 
    return expanded

#Initial_List(input_spec)
#expand(initial_list)

def MainFunction (spectrum):

    peptide_length, Tmplist = Initial_List (spectrum)

    for i in range (0, peptide_length-1, 1):
        Expand_List = expand(Tmplist)
        Tmplist=[]
        for j in range (0, len(Expand_List) ,1):
            if( IsConsistent(Expand_List[j], spectrum) == True):
                Tmplist.append(Expand_List[j])


    for x in range (0, len(Tmplist),1):
        print(Tmplist[x])
        print()

    return Tmplist

#Main_Function(input_spec)

def InputWindow ():

    choice = input("Enter (S) for Theoretical Spectrum or (C) for Cyclopeptide Sequencing : ")

    if (choice == 'S'):
        AA_Seq = input("Enter an Aminoacid Sequence: ")
        AA_Size = len(AA_Seq)
        Circular_Spectrum(AA_Seq)

    elif  (choice == 'C'):

            MainFunction(input_spec)
        

    

    return 0

InputWindow()
