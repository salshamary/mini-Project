from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import os

import argparse

mydir = os.getcwd()
os.chdir(mydir)
logfile = open('miniProject.log', 'a')
cds_outfile = open('EF999921_CDS.fasta', 'w')
fasta_outfile = open('EF999921.fasta', 'w')


def Initial(SRR):
    wget_1 = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    fastq_files = 'fastq-dump -I --split-files ' + SRR + '.1'
    os.system(wget_1)
    os.system(fastq_files)

def CDS():
    logfile = open('miniProject.log', 'w')
    # outfile for and CDS

    # accessing Entrtez and retrieving the records
    Entrez.email = 'salshamary@luc.edu'
    handle = Entrez.efetch(db='nucleotide', id=' EF999921', rettype='fasta')
    records = list(SeqIO.parse(handle, 'fasta'))
    # create fasta file and writes record to it

    fasta_outfile.write('>' + str(records[0].description) + '\n' + str(records[0].seq))

    # access the genbank recods and find CDS for id
    genCDS = Entrez.efetch(db='nucleotide', id='EF999921', rettype='gb', retmode='text')
    cds = 0
    # for loop to count CDS from genbank
    for i in SeqIO.parse(genCDS, 'genbank'):
        for j in i.features:
            if j.type == 'CDS':
                cds += 1
                cds_outfile.write('>' + str(j.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'",
                                                                                                                  "") + '\n' + str(
                    j.location.extract(i).seq) + '\n')
    cds_outfile.close()
    logfile.write('The HVFC Genome (EF999921) has ' + str(cds) + ' CDS.' + '\n')

    fasta_outfile.close()


def kallisto(SRR):
    kallisto_cmd = 'time kallisto index -i HCMVindex.idx EF999921_CDS.fasta'
    os.system(kallisto_cmd)
    kallisto_run = 'time kallisto quant -i HCMVindex.idx -o ./' + str(SRR) + ' -b 30 -t 4 ' + str(SRR) + '.1_1.fastq ' + str(SRR) + '.1_2.fastq'
    os.system(kallisto_run)'''
#runs sleuth by calling the Rscript and reading the lines into the logfile
def Sleuth():
    runSleuth = 'Rscript sleuth.R'
    os.system(runSleuth)
    output = open('sleuth_output.txt','r')
    #reads in each line of sleuth to the outfile
    outputlines= output.readlines()
    for i in outputlines:
        logfile.write(i)


def bowtie(SRR):
    build_cmd = 'bowtie2-build ./EF999921.fasta EF99992_1'
    os.system(build_cmd)
    bowtie_cmd = 'bowtie2 --quiet --no-unal --al-conc BOW_'+SRR+'.fastq -x EF99992_1 -1 '+SRR+ '.1_1.fastq -2'+SRR+'.1_2.fastq -S '+SRR+ '.sam'
    os.system(bowtie_cmd)
def bowReads(SRR):
    logfile = open("miniProject.log", "a")
    donorNum = ""
    #set donor name based of off SRR
    if SRR == SRR[0]:
        donorNum = "Donor 1 (2dpi)"
    elif SRR == SRR[1]:
        donorNum = "Donor 1 (6dpi)"
    elif SRR == SRR[2]:
        donorNum = "Donor 3 (2dpi)"
    else:
        donorNum = "Donor 3 (6dpi)"
    #look at reads from before bowtie
    SRR_1 = open(str(SRR)+".1_1.fastq")
    SRR_2= open(str(SRR)+".1_2.fastq")
    count1 = 0
    count2 = 0
    #count reads
    for line in SRR_1:
        count1+=1
    for line in SRR_2:
        count2+=1
    #count of before
    count =(count1+ count2)/8
    bowfile1 = open("BOW_"+SRR+".1.fastq")
    bowfile2 = open("BOW_"+SRR+".2.fastq")
    #count reads after bowtie
    afcount1 = 0
    afcount2 = 0
    for line in bowfile1:
        afcount1 +=1
    for line in bowfile2:
        afcount2 +=1
    #the count after
    afcount =(afcount1+ afcount2)/8
    with open("miniProject.log", "a") as logfile:
    	logfile.write((str(donorNum)+ " had " +str(count) + " read pairs before Bowtie2 filtering and "+ str(afcount)+" pairs after." + "\n"))
    print(SRR)

##runs spades on all SRRs and writes the command used to the logfile

def spades(SRR1,SRR2, SRR3,SRR4):
    command = 'spades -k 55,77,99,127 -t 2 --only-assembler --phred-offset 33  --pe1-1 BOW_' + SRR1 + '.1.fastq --pe1-2 BOW_' + SRR1 + '.2.fastq --pe2-1 BOW_' + SRR2 + '.1.fastq --pe2-2 BOW_' + SRR2 + '.2.fastq --pe3-1 BOW_' + SRR3 + '.1.fastq --pe3-2 BOW_' + SRR3 + '.2.fastq --pe4-1 BOW_' + SRR4 + '.1.fastq --pe4-2 BOW_' + SRR4 + '.2.fastq -o Spades/'
    os.system(command)
    with open("miniProject.log", "a") as logfile:
   	 logfile.write(command + "\n")

def contigs():
    contigCount = open("contigs.txt", "w")
    logfile = open("miniProject.log", "a")
    count = 0
    handle = SeqIO.parse("./Spades/contigs.fasta", "fasta")
    #look at the length of each sequence and count how many contigs have more than 1000 bp
    for record in handle:
        recordLen = len(record.seq)
        if recordLen > 1000:
            #appends the id and sequence of the contigs that are greater than 1000 to be looked at in the next function

            count += 1
            contigCount.write('> ' + str(record.id) + '\n' + str(record.seq) + '\n')
    contigCount.close()
    with open("miniProject.log", "a") as logfile:

    	logfile.write("There are " + str(count) + " contigs > 1000 bp in the assembly." + "\n")
def bp_contigs():
    #opens the logfile as well as the file that contains the contigs with more than 1000 bp

    logfile = open('miniProject.log', 'a')
    contigCount = open("contigs.txt", "r")
    handle = SeqIO.parse("contigs.txt", "fasta")
    seqLength = []
    total = 0
    #find the length of each sequence
    for record in handle:
        recordLen = len(record.seq)
        seqLength.append(int(recordLen))
    #for loop to total the lengths of the sequences and find the total bp in te assembly
    for i in range(0, len(seqLength)):
        total = total + seqLength[i]
    with open("miniProject.log", "a") as logfile:
    	logfile.write("There are " + str(total) + " bp in the assembly." + "\n")
def blast():
    logfile = open('miniProject.log', 'w')
    fasta = open("assembly.fasta").read()
    handle = NCBIWWW.qblast("blastn", "nr", fasta, entrez_query='"Betaherpesvirinae "[organism]')
    with open("blast.xml", "w") as outhandle:
        outhandle.write(handle.read())
    outhandle.close()
    blastResult = SearchIO.read("blast.xml", "blast-xml")
    # write the headers 
    logfile.write(('sacc' + '\t' + 'pident' + '\t' + 'length' + '\t' + 'qstart' + '\t' + 'qend' + '\t' + 'sstart' + '\t' + 'send' + '\t' + 'bitscore' + '\t' + 'eval' + '\t' + 'stitle')

if __name__ == '__main__':
    logfile = open("miniProject.log", "a")
    parser = argparse.ArgumentParser(description='SRRs')
    parser.add_argument('SRR', metavar='N', type=str, nargs='+', help='SRR values')
    args = parser.parse_args()
    logfile.write('SRA values tested: ' + str(args.SRR) + '\n')
    CDS()
    for i in args.SRR:
        kallisto(i)
    for i in args.SRR:
        bowtie(i)
    for i in args.SRR:
        bowReads(i)
    spades(args.SRR[0], args.SRR[1], args.SRR[2], args.SRR[3])  
    contigs()
    bp_contigs()
