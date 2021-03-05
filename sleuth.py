def SleuthInput(SRR):
    SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
    cond1 = "2dpi"
    cond2 = "6dpi"
    #initial line in file
    sleuthIN.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to outnput file
    for i in SRR:
        if int(i[3:])%2==0:
            sleuthIN.write(str(i)+ '\t' + cond1 + '\t'+ str(path)+ '\n')
        else:
            sleuthIN.write(str(i)+ '\t' + cond2 + '\t'+ str(path)+ '\n')
