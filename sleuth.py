def SleuthInput(SRR):
    SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
    cond1 = "2dpi"
    cond2 = "6dpi"
    #initial line in file
    sleuthIN.write('sample'+ '\t' + 'condition' + '\n')
    #based on SRR number, write condition to outnput file
    for i in SRR:
        if int(i[3:])%2==0:
            sleuthIN.write(str(i)+ '\t' + cond1 + '\n')
        else:
            sleuthIN.write(str(i)+ '\t' + cond2 + '\n')
