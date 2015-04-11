# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:26:39 2015

@author: buddingscientist
"""
import re

import os.path

import sys

from itertools import takewhile

import tkFileDialog

def command(refcount, filename, sequenceID, featCount):
    while True:
        userInput = raw_input('R:References S:Sequence M:Motif T:Translate F:Features E:Export Q:Quit >')
        if userInput == 'Q':
            quitOpt = raw_input('Would you like to:\n(1)quit\n(2)load another file? >')
            if quitOpt == '1':            
                quit()
            if quitOpt == '2':
                filename = tkFileDialog.askopenfilename()
                command(refcount, filename, sequenceID, featCount)
                break
            else:
                print 'Not a valid option...'
                command(refcount, filename, sequenceID, featCount)
        if userInput == 'R':
            references(refcount, filename, sequenceID, featCount)
        if userInput == 'S':
            sequence(filename, refcount, sequenceID, featCount)
        if userInput == 'M':
            motifFinder(filename, refcount, sequenceID, featCount)
        if userInput == 'T':
            translate(filename, refcount, sequenceID, featCount)
        if userInput == 'F':
            features(filename, refcount, sequenceID, featCount)
        if userInput == 'E':
            export(filename, refcount, sequenceID, featCount)





def references(refcount, filename, sequenceID, featCount):
    file2 = open(filename, 'r')
    file3 = file2.read()
    file4 = open(filename, 'r')
    j = 0
    reflist = []
    refcount=0
    for line in file4:
            if 'AUTHORS' in line:
                line3 = line.replace('AUTHORS', '')
                firstDot = line3.index('.')
                firstAuthor = line3[:firstDot]
                j+=1
                reflist.append((j, firstAuthor))
                refcount+=1
    print 'There are ' + str(refcount) + ' references for the sequence ' + sequenceID + ':'
    stringList =[]    
    for items in reflist:
        numName = str(items[0])+ ')' + str(items[1])
        stringList.append(numName)
    for item1 in stringList:
        print item1
    
    refDetails = [] 
    toPrint = file3.count('AUTHOR')
    count=0
    index = -1
    endIndex = -1
    while count < toPrint:
        index = file3.index('AUTHOR',index+1)
        endIndex = file3.index('REFERENCE',endIndex+1)
        refDetails.append(endIndex)
        refDetails.append(index)
        count+=1
    refDetails.append(file3.index('FEATURES'))
    
    start=1
    finalList=[]
    while start < len(refDetails):
        refEntry = file3[refDetails[start]:refDetails[start+1]]
        finalList.append(refEntry)
        start+=2
    showReference = raw_input('Enter a reference number to see details, press M to return to the main menu: ')
    try:
        num = int(showReference)-1
        print finalList[num]
    except ValueError:
        if showReference == 'M':
            command(refcount, filename, sequenceID, featCount)
        else:
            print 'Not a valid input...'
            command(refcount, filename, sequenceID, featCount)
    
def sequence(filename,refcount, sequenceID, featCount):
    file4 = open(filename, 'r')
    READ = file4.read()
    if 'ORIGIN' in READ:
        beginning = READ.index('ORIGIN')
        ending = READ.index('//')
        sequence = str(READ[beginning: ending])
        sequence = sequence.replace('ORIGIN','')
        print 'Enter a comma separated range in brackets to see that range. Square brackets will include positions in the range, round brackets will exclude it, first position in sequence is 1, not 0: '
        realSeq=''
        for letter in sequence:
            if letter == 'a':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter =='c':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 'g':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 't':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            else:
                continue
        seqChoice = raw_input('Press M to return to main menu')
        if seqChoice == 'M':
            command(refcount, filename, sequenceID, featCount)
        else:
            ranger = seqChoice.split(',')
            if str(ranger[0]).startswith('('):
                start = int(str(ranger[0]).replace('(',''))
            elif str(ranger[0]).startswith('['):
                start = int(str(ranger[0]).replace('[',''))-1
            else:
                print 'Not a valid range...'
                command(refcount, filename, sequenceID, featCount)
        
            if ')' in str(ranger[1]):
                end = int(str(ranger[1]).replace(')',''))-1
            elif ']' in str(ranger[1]):
                end = int(str(ranger[1]).replace(']',''))
            else:
                print 'Not a valid range...'
                command(refcount, filename, sequenceID, featCount)
            try:    
                userSeq = realSeq[start:end]
                splitIntoChunks = []
                u=59
                v=0
                while u < len(userSeq):
                    chunk = userSeq[v:v+u]
                    splitIntoChunks.append(chunk)
                    v+=59
                    u+=59
                for chunks in splitIntoChunks:
                    print chunks
            except IndexError:
                print 'Not in range...'
                command(refcount, filename, sequenceID, featCount)
            
    else: print 'This Genbank file has no sequence data'

def motifFinder(filename, refcount, sequenceID, featCount):
    askMotif = raw_input('Enter a motif <of at least 5 bases long> to search for it in the sequence (The wildcard character ‘?’ represents any nucleotide in that position, and * represents none or many nucleotides in that position.) You cannot use * as the first character: ')
    listMotif= []    
    letterlist = ['A','C','G','T', 'a', 'c','g','t']
    if len(askMotif) >= 5:    
        for letter in askMotif:
            if letter in letterlist:
                a = letter.capitalize()
                listMotif.append(a)
            if letter == '?':
                listMotif.append('.')
            if letter == '*':
                listMotif.append('.*?')
    else: command(refcount, filename, sequenceID, featCount)
    pattern = ''
    for searcher in listMotif:
        pattern+=searcher
    file4 = open(filename, 'r')
    READ = file4.read()
    if 'ORIGIN' in READ:
        beginning = READ.index('ORIGIN')
        ending = READ.index('//')
        sequence = str(READ[beginning: ending])
        sequence = sequence.replace('ORIGIN','')
        realSeq=''
        for letter in sequence:
            if letter == 'a':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter =='c':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 'g':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 't':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            else:
                continue
    else: print 'This Genbank file has no sequence data'
    try:    
        match = re.search(r'{0}'.format(pattern), realSeq)
        listmatch = match.group()
        length = len(listmatch)
        index = realSeq.index(listmatch)
        length+=index
        printMatch = realSeq[index:length]
        leftside = index-6
        rightside = length+6
        print realSeq[leftside:index-1].lower() +'['+ str(index) + ']' + printMatch + '[' + str(length) + ']' + realSeq[length+1:rightside].lower()
        ask = raw_input('Would you like to\n(1)Enter another motif; or\n(2)Go to the main menu? >')
        if ask == '1':
            motifFinder(filename, refcount, sequenceID, featCount)
        elif ask == '2':
            command(refcount, filename, sequenceID, featCount)
        else:
            print 'Not a valid input...'
            command(refcount, filename, sequenceID, featCount)
    except AttributeError:
        print 'no match'
        command(refcount, filename, sequenceID, featCount)
        
def translate(filename, refcount, sequenceID, featCount):
    file_ = open(filename,'r')
    READ = file_.read()
    if 'ORIGIN' in READ:
        beginning = READ.index('ORIGIN')
        ending = READ.index('//')
        sequence = str(READ[beginning: ending])
        sequence = sequence.replace('ORIGIN','')
        realSeq=''
        for letter in sequence:
            if letter == 'a':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter =='c':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 'g':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 't':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            else:
                continue
    else: print 'This Genbank file has no sequence data'
    try:    
        userRange = raw_input('Enter a range of the seq to see its amino acid translation. OPTIONS: \n1) \'FULL\' will translate the whole sequence. \n2) Specifying the ORF will print only that ORF, specifying no ORF will print all three.')
        if userRange == 'FULL':
            dictionary1 = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
            bases = list(realSeq) 
            bases = [dictionary1[base] for base in bases] 
            string = ''.join(bases)
            complement = string[::-1]
            start1=0
            a=True
            while a == True:
                stop_codons = ('TAA', 'TGA', 'TAG')
                start1 = complement.find('ATG',start1+1)
                trimmed_sequence = complement[start1:]
                codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]
                coding_sequence  =  takewhile(lambda x: x not in stop_codons and len(x) == 3 , codons)
        
                codontable = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
                }
                protein_sequence = ''.join([codontable[codon] for codon in coding_sequence])
                        
                if len(protein_sequence)>0:
                    print protein_sequence
                else:
                    a=False
                
        else:
            ranger = userRange.split(',')
            if str(ranger[0]).startswith('('):
                start = int(str(ranger[0]).replace('(',''))
            elif str(ranger[0]).startswith('['):
                start = int(str(ranger[0]).replace('[',''))-1
            else:
                print 'Not a valid range...'
                command(refcount, filename, sequenceID, featCount)
        
            if ')' in str(ranger[1]):
                end = int(str(ranger[1]).replace(')',''))-1
            elif ']' in str(ranger[1]):
                end = int(str(ranger[1]).replace(']',''))
            else:
                print 'Not a valid range...'
                command(refcount, filename, sequenceID, featCount)
            dictionary1 = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
            bases = list(realSeq[start:end]) 
            bases = [dictionary1[base] for base in bases] 
            string = ''.join(bases)
            complement = string[::-1]
            start1=0
            b=True
            while b==True:
                stop_codons = ('TAA', 'TGA', 'TAG')
                start1 = complement.find('ATG',start1+1)
                trimmed_sequence = complement[start1:]
                codons = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]
                coding_sequence  =  takewhile(lambda x: x not in stop_codons and len(x) == 3 , codons)
        
                codontable = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
                }
                protein_sequence = ''.join([codontable[codon] for codon in coding_sequence])
           
                if len(protein_sequence)>0:
                    print protein_sequence
                else:
                    b=False
            
    except IndexError:
        print 'Not in range...'
        command(refcount, filename, sequenceID, featCount)
    
    
def export(filename, refcount, sequenceID, featCount):
    file2 = open(filename,'r')
    READ =file2.read()
    if 'ORIGIN' in READ:
        beginning = READ.index('ORIGIN')
        ending = READ.index('//')
        sequence = str(READ[beginning: ending])
        sequence = sequence.replace('ORIGIN','')
        realSeq=''
        for letter in sequence:
            if letter == 'a':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter =='c':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 'g':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            elif letter == 't':
                ALLCAPS = letter.capitalize()
                realSeq += ALLCAPS
            else:
                continue
    userFasta = raw_input('Please enter the desired name for your fasta file, including extension> ')
    if userFasta == '':
        command(refcount, filename, sequenceID, featCount)
    elif os.path.exists(userFasta)  == True:
        askOverwrite = raw_input('Filename exists... Overwrite current file? (Y/N) ')
        if askOverwrite == 'Y' or 'y':
            file_ = open(userFasta, 'a')
            file_.write('>')
            file_.write(sequenceID)
            file_.write('\n')
            file_.write(realSeq)
            file_.close()
        elif askOverwrite == 'N' or 'n':
            command(refcount, filename, sequenceID, featCount)
    else:
        file_ = open(userFasta, 'a')
        file_.write('>')
        file_.write(sequenceID)
        file_.write('\n')
        file_.write(realSeq)
        file_.close()
        
def features(filename, refcount, sequenceID, featCount):
    
    
    k=0
    record = open(filename, 'r')
    record2 = open(filename,'r')
    READ = record2.read()
    featureList = []
    featList = []
    print 'There are ' + str(featCount) + ' features... Enter a number to view the features: '
    chooseFeature = input()
    if chooseFeature > featCount:
        print 'not in range...'
        command(refcount, filename, sequenceID, featCount)
   
    else:
        
        for line in record:
            line1 = line.strip()
            if line1.startswith('gene') or line1.startswith('CDS') or line1.startswith('source') or line1.startswith('ORIGIN'):
                append = READ.index(line1)
                
                featureList.append(append)
                
                
        while k< len(featureList)-1:
            featList.append(READ[featureList[k]:featureList[k+1]])
            k+=1
        print featList[chooseFeature-1]
            



def initialize():
    filename = sys.argv[1]
    file1 = open(filename, 'r')
    refcount = 0
    featCount = 0
    for line in file1:
        if line.startswith('LOCUS'):
            header = ','.join(line.split())
            header = header.split(',')
            sequenceID = header[1]
            basepair = header[2]
        if line.startswith('ACCESSION'):
            getAccession = ','.join(line.split())
            getAccession = getAccession.split(',')
            accession = getAccession[1]
        if line.startswith('SOURCE'):
            organism = line.replace('SOURCE','')
            organism = organism.strip()
        if line.startswith('REFERENCE'):
                refcount+=1
        line1 = line.strip()
        if line1.startswith('gene') or line1.startswith('CDS') or line1.startswith('source'):
            featCount+=1
    
        
    print '============================================================'
    print  'GBK Reader -('+filename+') - DNA'
    print '============================================================'
    print 'Sequence: '+ accession + ' (' + sequenceID + ', ' + basepair + ' bp)'
    print 'Description: ' + organism
    print 'Number of references: ' + str(refcount)
    print 'Number of features: ' + str(featCount)
    print '============================================================'
    command(refcount, filename, sequenceID,featCount)
    
initialize()