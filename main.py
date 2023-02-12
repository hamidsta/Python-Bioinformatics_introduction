from origin import *
import os,sys
import urllib



while True :
    xt=input("""
    1)Séquence complémentaire
    2)Récupération de la Séquence 3’UTR via une fiche genbank
    3)Recherche de motif exact
    4)Recherche de sites cibles d’un microARN dans un transcrit
    5)Récupération sur internet des séquences
    6) QUITTER

    CHOISIR FONCTION >>>>""")

#Partie 1########################################################################
    if xt=='1':
        print("Partie 1 :")                              
        x=str(input('mettre nom du fichier à analyser:'))
        if os.path.isfile(x)==True:
            header, sequence = readFasta(x)                                                     #lecture fichier fasta , recuperation header, sequence
            seqRev = reverse_complement(sequence)                                               # reverse complémentaire de la sequence
            write_on_file(seqRev,"Voici l'inverse complémentaire de la sequence {}".format(header))       # sauvegard dans un fichier
        else:
            print("fichier '{}' non existant ou pas dans le bon répertoire".format(x))

    #Partie 2#########################################################################
    if xt=='2':
        try:
            print("Partie 2 :")
            x=str(input('mettre nom du fichier GENBANK à analyser:'))
            if os.path.isfile(x)==True:
                seq_start, seq_end = recherche_positions(x)                                     #recupere position du CDS et LOCUS du GENBANK
                version=comparaison_version(x)                                                  #recupere Version du fichier
                y=str(input('mettre nom du fichier FASTA correspondant à analyser:'))
                if os.path.isfile(y)==True:
                    header_2,sequence_2 = readFasta(y)                                         # recuperation sequence, header du fichier fasta
                    header_list=header_2.split()                                               #split du header

                    if version in header_list[0]:                                               # si Version similaire entre FASTA et GENBANK
                        sequence_2 = sequence_2[seq_start:seq_end]                              #recuperation region 3'UTR du FASTA
                        print(header_list[0])
                        print(version)
                        write_on_file(sequence_2,"voici la région 3'UTR de {} " .format(header_2))
                    else:
                        print("chosir le fichier FASTA correspondant à votre fichier genbank")
                else:
                    print("fichier '{}' non existant ou pas dans le bon répertoire".format(y))
            else:
                print("fichier '{}' non existant ou pas dans le bon répertoire".format(x))    
        except IndexError:                                                                       #Si pas d'argument CDS\LOCUS dans la fiche genbank
            print("le fichier '{}' n_est pas un fichier GENBANK ou ne contient pas de CDS".format(x))    
            
        

    
#Partie 3#########################################################################
    if xt=='3':
        x=str(input('mettre nom du fichier MOTIF .fas à analyser:'))
        if os.path.isfile(x)==True:
            motif = readFasta(x)[1]                                                      #recupere deuxieme ligne du fichier (=motif)
            y=str(input('mettre nom du fichier SEQUENCE .fas à analyser:'))
            if os.path.isfile(y)==True:                                              
                seqN= readFasta(y)[1]                                                        #recupere deuxieme ligne du fichier (=sequence)
                res, nombre_occ = recherche_All_motif(seqN, motif)
                if res ==[]:                                                                 # liste vide => pas de resultat
                    print('pas de chance, pas de motif', res) 
                else:                         
                    print(" voici les positions et le nombre d'occurence " , res, nombre_occ)
            else:
                print("fichier '{}' non existant ou pas dans le bon répertoire".format(y))
        else:
            print("fichier '{}' non existant ou pas dans le bon répertoire".format(x))

    #Partie 4################################################################
    if xt=='4':
        try:                                                                                        
            x=str(input('mettre nom du fichier GENBANK à analyser:'))               
            if os.path.isfile(x)==True:
                seq_start, seq_end = recherche_positions(x)
                version=comparaison_version(x)
                y=str(input('mettre nom du fichier FASTA correspondant à analyser:'))
                if os.path.isfile(y)==True:
                    header_2,sequence_2 = readFasta(y)
                    header_list=header_2.split()
                    if version in header_list[0]:
                        sequence_2 = sequence_2[seq_start:seq_end]
                        sequence_inverse = reverse_complement(sequence_2)
                        z=str(input('mettre nom du fichier micro ARN à analyser:'))
                        if os.path.isfile(z)==True:
                            header, sequence = readFasta(z)
                            seed_region=sequence[1:7]
                            res , nombre_occ=recherche_All_motif(sequence_2,seed_region)                 # recupere du 2 eme au 7eme nucleotide
                            if res ==[]:
                                print('pas de chance, pas de motif', res) 
                            else:
                                print("voici le resultat", res )
                                
                        else:
                            print("fichier '{}' non existant ou pas dans le bon répertoire".format(z))
                    else:
                        print("chosir le fichier FASTA correspondant à votre fichier genbank")
                else:
                    print("fichier '{}' non existant ou pas dans le bon répertoire".format(y))
            else:
                print("fichier '{}' non existant ou pas dans le bon répertoire".format(x))    
        except IndexError:           
            print("le fichier '{}' n_est pas un fichier GENBANK ou ne contient pas de CDS".format(x))
       
       
        
    #Partie 5 ##################################################################
    if xt=='5':
        yt = None
        while yt not in ['1', '2'] :
            yt=yt=input("""
                1)MIRBASE
                2)NCBI
        
                CHOISIR PLATEFORME >>>>""")

        x=str(input('mettre ID de la sequence à telecharger:'))
        
        if yt == '1':
            URL = f"http://www.mirbase.org/cgi-bin/get_seq.pl?acc={x}"
            download_sequence(URL, 'fichier_mirbase.fas', yt)
        
        else:
            y=str(input('mettre l\'extension de la sequence à telecharger: fasta ou gb '))
            URL = f"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={x}&rettype={y}&retmode=text"
            if y == 'fasta' :                                               
                download_sequence(URL, 'fichier_fasta.fas', yt)
            elif y == 'gb':
                download_sequence(URL, 'fichier_gb.gb', yt)
            else :
                print("Erreur : Extension non valide")
        
            
    if xt=='6':
            break 
            
    else: 
        print("Vous devez choisir une selection entre 1 et 5")