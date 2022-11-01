
'''

@function : Reading A FASTA file :
'''


import os,sys
from urllib.request import urlopen
import urllib

def readFasta (seqFile) :                         
    if os.path.isfile(seqFile):          #test si (seqfile) proposé correspond bien à un fichier
        f = open(seqFile, 'r')
        text = f.read().split('\n')      # lecture du fichier 
        header = text[0][1:]             # Le header
        sequence = ''.join(text[1:])     # La séquence
        f.close()
        return header, sequence                 
    else:
        print("Erreur : le chemin donne '{}' n'est pas valide.".format(seqFile))    #si le fichier n'est pas trouvé message d'erreur  sortie brusque 
        sys.exit()
        

'''
'''
'''
reverse transcribe a  DNA sequence
'''
def reverse_complement (seq):                                            
    """ Computes the reverse complement of the DNA sequence. """
    if not type(seq) is str :                                      
        seq="".join(map(str, seq))                 #devient chaine de caractere 
    seq=seq.replace(" ","")                         #supprime espace
    comp = ""                                          # variable str vide
    for c in seq:
        if c == "A"or c=='a':
            comp = "T" + comp
        elif c == "T" or c=='t':
            comp = "A" + comp
        elif c=="U" or c=='u':
            comp='T' +comp
        elif c == "G" or c=='g':
            comp = "C" + comp
        elif c== "C" or c=='c':
            comp = "G" + comp
        else:
            return("Votre entrée n'est pas une sequence nucleotidique ")
    return comp


'''
Recherche de motif exact
'''
def recherche_All_motif(seq, motif):                         
    res = []
    nombre_occurence = 0               
    for i in range(len(seq) - len(motif) + 1):
        j = 0                                                
        while j < len(motif) and motif[j] == seq[i+j]:    # si valeur motif [0]==seq[0] 
            j += 1                                        #  déplace à motif[1]==seq[1]   
        if j == len(motif):                               # si j = longeur motif ==> motif trouvé 
            res.append(i)
            nombre_occurence += 1
    return res,nombre_occurence


def download_sequence(URL, fileName, URLType):        
    fileObj = open(fileName, 'w')                      # rajoute dans doc
    a = 0
    try :
        response = urlopen(URL)
    
    except IOError:
        print("Erreur : Pas de connexion internet") #si pas d'internet
        a+=1
    
    except urllib.error.URLError :
        print("Erreur : l'ID donne n'est pas valide.")    #si l'ID n'est pas valide, message d'erreur puis sortie brusque
        a+=1
    
    finally:
        if a == 0 :                               
            if URLType == '1':                 #fichier mirbase
                data = response.readlines()[1:-1]    # on renvoie la sequence sans les </pre>
            else :                                  #fichier GB\fasta possede pas </pre>
                data = response.readlines()
                
            if len(data)==0:
                print("Erreur : l'ID donne n'est pas valide.")
                
            else :
                for line in data:
                    line = line.strip()
                    line = line.decode('utf-8')
                    print(line)
                    fileObj.write(line+'\n')
                
            fileObj.close()                          #fermeture
    


def write_on_file(res, verbose=""):             
    if not isinstance(res, str):            
        res = str(res)                      #devient une chaine de caracterere 
    if len(verbose):                        
        print(verbose)
    print(res)
    if not input("Souhaitez vous enregistrer le resultat dans un nouveau fichier ? (oui/o | n)").lower() in ("oui", "o"): 
        return
    file_save = input("Nom du fichier ou le resultat va etre enregistrer : (enregistrer en .fas si aucune extension donner) ")
    if not "." in file_save:                               # ajout de .fas par défault 
        file_save += ".fas"                                
    with open(file_save, "w+", encoding="utf8") as f:        #on ouvre le fichier précedemment crée , on active le mode ecriture dans le fichier précedemment crée , encodé en utf8 pour éviter les erreurs d'encodage en fonction du terminal/systeme de lecture utilisé 
        f.write(res)                                        # on écrit dans le fichier précedemment crée
    print(f"Resultat enregistrer dans {file_save}") 
    
def recherche_positions(gbFile):
    f = open(gbFile, 'r')
    text = f.read()
    f.close()
    
    locus_pos=text.find('LOCUS') #Trouver le début de la partie LOCUS dand le fichier entier
    locus_text=text[locus_pos:] 
    
    features_pos = text.find("FEATURES") #Trouver le début de la partie FEATURES dand le fichier entier
    features_text = text[features_pos:]
    
    CDS_pos = features_text.find("CDS")
    CDS_text = features_text[CDS_pos:]
    

    locus_join_last_line = locus_text.split('\n')[0]       #Trouver la ligne LOCUS


    
    locus_join_last_seq = locus_join_last_line.split()[2]   #Trouver la valeur dans la ligne LOCUS (Fin de séquence)
    
    if "join" not in CDS_text.split("/")[0]:
        CDS_join_line = CDS_text.split("/")[0]                    #Trouver les dernières valeurs de CDS
        CDS_join_last_seq = CDS_join_line.split()[-1]             #Trouver la dernière valeur de CDS
        seq_start = int(CDS_join_last_seq.split('.')[-1])         #Trouver la position du début de la séquence (valeur en entier)
    else:
        CDS_join_line = CDS_text.split("/")[0]                     #Trouver le CDS join
        CDS_join_last_seq = CDS_join_line.split()[-1]              #Trouver les dernières valeurs de CDS join
        seq_start = int(CDS_join_last_seq.split('.')[-1][:-1])     #Trouver la position du début de la séquence (valeur en entier)
    
    seq_end = int(locus_join_last_seq)                             #Trouver la position de fin de la séquence (valeur en entier)
    
    return seq_start, seq_end
    
def comparaison_version(gbFile):
    f = open(gbFile, 'r')
    text = f.read()
    f.close()
    
    version=text.find('VERSION') #Trouver le début de la partie VERSION dand le fichier entier
    version_text=text[version:]  # recuperer la ligne correspondante 
    
    version_join_last_line = version_text.split('\n')[0]                   
    version_join_number = version_join_last_line.split()[1] 
                     
    return version_join_number


