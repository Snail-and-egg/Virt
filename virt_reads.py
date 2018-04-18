def do_reads1(name='EColi_O157',r_len=100,cov=1,dist=[10,20]):
    from scripts import ref
    from math import ceil
    from random import randint
    version = 0.3
    genome = ref(name)
    g_len = len(genome)
    print(g_len)
    check = 0
    min_overlap = max(dist[0],ceil(r_len*0.03))
    max_overlap = max(dist[1],ceil(r_len*0.12))
    with open(name+'_virt1_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'.fasta','w') as f:
        with open(name+'_virt1_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'.metadata.txt','w') as g:
            g.write('version = '+str(version)+'\nname = '+name+'\ngenome legth = '+str(g_len)+'\ncoverage = '+str(cov)+
                    '\nread legth = '+str(r_len)+'\ndistance = ['+str(dist[0])+','+str(dist[1])+']\n')
            for i in range(cov):
                begin = randint(0,max_overlap)
                while begin <= g_len-r_len:
                    check += 1
                    f.write('>#'+str(check)+' | pos='+str(begin)+'\n')
                    f.write(genome[begin:begin+r_len]+'\n')
                    begin = begin + r_len - randint(min_overlap,max_overlap)

def do_reads2(name='EColi_O157',r_len=100,cov=1,error=0.05, paired = True, gap=210, step=1000):
    from scripts import ref, gc_content
    from math import ceil, fabs
    from random import randint, random, normalvariate
    version = 0.2
    genome = ref(name)
    g_len = len(genome)
    print(g_len)
    M = round(g_len/(r_len-r_len*error) * cov)
    comp = {'A':'T','C':'G','G':'C','T':'A'}
    print('Starting read modeling')
    with open(name+'_virt2_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'.metadata.txt','w') as g:
        g.write('version = '+str(version)+'\nname = '+name+'\ngenome legth = '+str(g_len)+
                '\ncoverage = '+str(cov)+'\nread legth = '+str(r_len))
        gc = gc_content(name, genome, win=step, step=step, prog = 'C+G', save='return')
        if paired == False:
            with open(name+'_virt2_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'.fasta','w') as f:
                for t in range(0,g_len//step):
                    M = ceil(step/(r_len-r_len*error) * cov *(1-fabs(0.5-gc[t])*2))
                    for i in range(M):
                        begin = randint(step*t,step*t+step)
                        if round(random()) == 0:
                            f.write('>#'+str(i)+' | pos='+str(begin)+'\n')
                            f.write(genome[begin:begin+r_len]+'\n')
                        else:
                            f.write('>#'+str(i)+' | pos='+str(begin)+'\n')
                            f.write(''.join([comp[i] for i in genome[begin+r_len-1:begin-1:-1]])+'\n')
        else:
            with open(name+'_virt2_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'1.fasta','w') as f:
                with open(name+'_virt2_'+str(r_len)+'_'+str(cov)+'_'+str(version)+'2.fasta','w') as p:
                    for t in range(0,g_len//step):
                        M = ceil(step/(r_len-r_len*error) * cov *(1-fabs(0.5-gc[t])*2))
                        for i in range(M//2):
                            began = round(normalvariate(gap,gap*0.05))
                            #began = gap
                            begin = randint(step*t,step*t+step)
                            f.write('>#'+str(i)+' | pos='+str(begin)+'\n')
                            f.write(genome[begin:begin+r_len]+'\n')
                            p.write('>#'+str(i)+' | pos='+str(begin+began)+'\n')
                            p.write(''.join([comp[i] for i in genome[began+begin+r_len-1:began+begin-1:-1]])+'\n')
                        #p.write(genome[began+begin:began+begin+r_len]+'\n')
                    
if __name__ == '__main__':
    step = 1000
    do_reads2(name='EColi_O157',r_len=200,cov=100,step=step,error=0.05, paired = True, gap=210)
'''
#   Анализ смоделированных ридов
    from scripts import cov_average, graph, l, sam
    name = 'EColi_O157_100_virt2_paired'
    step = 1000
    cov = cov_average(name=name,win=step,step=step,prog='return',repair='yes')
    graph('E_coli O157 virtual2 reads with coverage 100 paired, window '+str(step),'Genome',['Coverage'],
          ['cov_3'],[i for i in range(len(cov))],cov, typ = '-',two_axis='no',
          graphs=1,lines=1,show='save',log='no',figsize=(15,8),dpi=80,linewidth=1,
          color=['green','purple','orange','blue'])

    data = sam(name=name)
    l_full,l100,l_std = l(name, data, step=step, prog = 'return')
    graph('E_coli_O157_100_virt2_paired_win'+str(step)+'_l','Genome',['Read crossing length'],
          ['l'],[i for i in range(len(l_full))],l_full, typ = '-',two_axis='no',
          graphs=1,lines=1,show='save',log='no',figsize=(15,8),dpi=80,linewidth=1,
          color=['green','purple','orange','blue'])
    graph('E_coli_O157_100_virt2_paired_win'+str(step)+'_l100','Length, num',['Amount, percent of max'],
          ['l100'],[i for i in range(1,101)],[l100[i] for i in range(1,101)], typ = '-',two_axis='no',
          graphs=1,lines=1,show='save',log='no',figsize=(15,8),dpi=80,linewidth=1,
          color=['green','purple','orange','blue'])
    graph('E_coli_O157_100_virt2_paired_win'+str(step)+'_lstd','Genome',['Std of read crossing'],
          ['lstd'],[i for i in range(len(l_std))],l_std, typ = '-',two_axis='no',
          graphs=1,lines=1,show='save',log='no',figsize=(15,8),dpi=80,linewidth=1,
          color=['green','purple','orange','blue'])
'''

