def cov_average(name='EColi_O104',win=1000,step=1,prog='return',repair='yes'):
    print('Count average coverage '+name+' '+str(win)+'win '+str(step)+'step')
    cov = []
    with open(name+'.coverage','r') as f:
        for i in f:
            line = i.strip().split('\t')
            if repair == 'yes' and len(cov) < int(line[1])-1:
                cov += [0]*(int(line[1])-len(cov)-1)
            cov += [float(line[2])]
    per_num = []
    for i in range(0,len(cov), step):
        try:
            per_num.append(float(sum(cov[i:i+win]))/win)
        except IndexError:
            per_num.append(float(sum(cov[i:]))/len(cov[i:]))
    if prog == 'save':
        with open(name+'_coverage_per_'+str(win)+'_step_'+str(step)+'.txt','w') as f:
            for i in per_num:
                f.write(str(i) +'\n')
    elif prog == 'return':
        return per_num

def graph(title,x_label,y_label,legend,x,*y, typ = '-',two_axis='no',graphs=1,lines=1,show='save',log='no',figsize=(15,8),dpi=80,linewidth=1,color=['green','purple','orange','blue']):
    print('Graph '+title)
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})
    fig = plt.figure(num=None, figsize=figsize, dpi=dpi, facecolor='w',edgecolor='k')
    plt.subplots_adjust(left=0.07, right=0.93, bottom=0.09, top=0.94)
    if two_axis == 'yes':
        ax0 = fig.add_subplot(111)
        line0, = ax0.plot(x, y[0][0], '-', linewidth=1, color=color[0], \
                 label=legend[0])
        ax1 = ax0.twinx()
        line1, = ax1.plot(x, y[1][0], '-', linewidth=1, color=color[1], \
                 label=legend[1])
        plt.legend((line0,line1),tuple(legend),loc='best')
        ax0.set_title(title)
#        ax0.legend(loc = 'best')
        ax1.set_xlabel(x_label)
        ax1.set_ylabel(y_label[1])
        ax0.set_ylabel(y_label[0])
        ax0.axis([0, len(x), 0,max(y[0][0])*1])
        ax1.axis([0, len(x), 0,max(y[1][0])*4])
#        ax0.legend((line0, line1),tuple(legend), loc = 'best')
    elif graphs==1 and lines==1:
        ax = fig.add_subplot(111)
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label[0])
        line1, = ax.plot(x, y[0], typ, linewidth=1, color='green', label=legend[0])
        ax.legend(loc='lower right', numpoints=1)
#        ax_y=[0, len(x), 0,max(y[0])*0.3]
        ax_y = [min(x),max(x),min(y[0]),max(y[0])]
        ax.axis(ax_y)
    else:
        for i in range(0,graphs):
            exec("ax%s = fig.add_subplot(int(str(graphs)+'1'+str(i+1)))"%i)

            print(ax0)
            for j in range(0,lines):
                exec("line%s%s, = ax%s.plot(x, y[i][j], '-', linewidth=1, color=color[j], \
                     label=legend[j])"% (i,j,i))
                print(j)
            exec("ax%s.legend(tuple(legend),loc='best', numpoints=1); \
                 ax%s.set_title(title); ax%s.set_xlabel(x_label); \
                 ax%s.set_ylabel(y_label[i])"% (i,i,i,i))
            if log == 'yes':
                exec("ax%s.set_yscale('log')"% (i))
    if show == 'save':
        print('Saving...')
        plt.savefig(title)
    elif show == 'show':
        plt.show()
    plt.close()

def ref(name):
    print('Open reference '+name)
    genome = ''
    with open(name+'.fasta','r') as f:
        name = f.readline().strip()
        for i in f:
            genome += i.strip()
    return genome

def gc_content(name, genome, win=10000, step=1, prog = 'GCskew', save='save'):
    print(prog+' '+name+' '+str(win)+'win '+str(step)+'step')
    with open(name+'_'+prog+'_per_'+str(win)+'win_step'+str(step)+'.txt','w') as f: 
        if save == 'save':
            for i in range(win+1,len(genome),step):
                if prog == 'GCskew':
                    avg_g_c = genome.count('С',i-win,i)-genome.count('G',i-win,i)/(genome.count('С',i-win,i)+genome.count('G',i-win,i))
                    f.write(str(avg_g_c/win)+'\n')
                elif prog == 'C-G':
                    avg_g_c = genome.count('C',i-win,i) - genome.count('G',i-win,i)
                    f.write(str(avg_g_c/win)+'\n')
                elif prog == 'C+G':
                    avg_g_c = genome.count('C',i-win,i) + genome.count('G',i-win,i)
                    f.write(str(avg_g_c/win)+'\n')
        elif save == 'return':
            avg_g_c = []
            for i in range(win+1,len(genome),step):
                if prog == 'GCskew':
                    avg_g_c += [genome.count('С',i-win,i)-genome.count('G',i-win,i)/(genome.count('С',i-win,i)+genome.count('G',i-win,i))]
                elif prog == 'C-G':
                    avg_g_c += [(genome.count('C',i-win,i) - genome.count('G',i-win,i))/win]
                elif prog == 'C+G':
                    avg_g_c += [(genome.count('C',i-win,i) + genome.count('G',i-win,i))/win]
            return avg_g_c
                

def autocor(cov=[], title='',prog='file',save='save', step=100):
    print('Autocorrelation '+title+str(step))
    from scipy.stats.stats import pearsonr
    if prog == 'file':
        with open(title+'.txt','r') as f:
            for i in f:
                cov += [float(i.strip())]
    steps = [step*i for i in range(len(cov))]
    if save == 'save':
        with open('Autocor_for_'+file+'_'+str(step)+'step.txt','w') as f:
            for i in steps:
                f.write(str(pearsonr(cov[:len(cov)-i],cov[i:])[0])+'\n')
    elif save == 'return':
        for i in steps:
            print((str(pearsonr(cov[:len(cov)-i],cov[i:])[0])+'\n'))

def l_another(name,data=[],prog='save',step=1000):
    if prog == 'read':
        with open(name+'_L.txt','r') as f:
            l_full = []
            for i in f:
                l_full += [float(i.strip())]
        p = []
        for i in range(step,len(l_full),step):
            p += [sum(l_full[i-step:i])/step]
        return p
    l_full = {}
    print('Generate l_full var')
    for i in range(data[-1][0]+data[-1][1]+1):
        l_full[i] = [0,0]
    print('Count L') 
    for j in range(len(data)):
        n = j+1
        if n >= len(data):
            break
        while data[j][0] + data[j][1] > data[n][0]:
            leng = min(data[j][0]+data[j][1],data[n][0]+data[n][1])-data[n][0]
            for i in range(min(data[j][0]+data[j][1],data[n][0]+data[n][1])-(leng),data[j][0]+leng):
                l_full[i][0] += leng
                l_full[i][1] += 1
            n += 1
            if n >= len(data):
                break
    print('Saving')
    if prog == 'save':
        with open(name+'_L'+'.txt','w') as f:
            for i in range(max(l_full)):
                try:
                    f.write(str(l_full[i][0]/l_full[i][1])+'\n')
                except ZeroDivisionError:
                    f.write('0\n')
    elif prog == 'return':
        ret = []
        for i in range(step,max(l_full),step):
            try:
                ret += [sum([l_full[j][0] for j in range(i-step,i)])/sum([l_full[j][1] for j in range(i-step,i)])]
            except ZeroDivisionError:
                ret += [0]
        return ret
def l_only(name, data, step=10000, prog = 'return'):
    print('l_full '+name+' '+str(step))
    from math import sqrt
    print('Number of reads = ' + str(len(data)))
        
    #L
    print('Count L and 100L')
    l = 0 
    k = 0
    t = 1
    l_full = []
    for j in range(len(data)):
        n = j+1 
        if data[j][0] >= step*t:
            print(str(t*step))
            t += 1
            try:
                l_full += [l/float(k)]
            except ZeroDivisionError:
                print('ZeroDivisionError',str(l),str(k), str(j), data[j])
                l_full += [0]
            l = 0
            k = 0
        if n >= len(data):
            break
        while data[j][0] + data[j][1] > data[n][0]:
                l += data[j][0] + data[j][1] - data[n][0]
                k += 1
                n += 1
                if n >= len(data):
                    break
    return l_full
                
def l(name, data, step=10000, prog = 'return'):
    print('l_full '+name+' '+str(step))
    from math import sqrt
    print('Number of reads = ' + str(len(data)))
        
    #L
    print('Count L and 100L')
    l = 0 
    k = 0
    t = 1
    #L100 + L      
    l100 = {}
    l_full = []
    with open(name+'_L_for_'+str(step)+'step.txt','w') as f:
        for j in range(len(data)):
            n = j+1 
            if data[j][0] >= step*t:
#                print("Avg for " + str(t*step) + " n.u. = " + str(l/float(k)))
                t += 1
                try:
                    f.write(str(l/float(k))+'\n')
                    l_full += [l/float(k)]
                except ZeroDivisionError:
                    print('ZeroDivisionError',str(l),str(k), str(j), data[j])
                    f.write('0\n')
                    l_full += [0]
                l = 0
                k = 0
            if n >= len(data):
                break
            while data[j][0] + data[j][1] > data[n][0]:
                try:
                    l100[data[j][0] + data[j][1] - data[n][0]] += 1
                except KeyError:
                    l100[data[j][0] + data[j][1] - data[n][0]] = 1
                l += data[j][0] + data[j][1] - data[n][0]
                k += 1
                n += 1
                if n >= len(data):
                    break
    m = float(max(l100.values()))
    print('Max 100L = ' + str(m), 'Sum 100L = '+str(sum(l100.values())))
    print('Save 100L norming file')
    with open(name+'_100_L_average_'+str(step)+'step.txt','w') as f:
        for key, value in l100.items():
            f.write(str(value/m)+'\n')
            l100[key] =  value/m

    #L std
    print('Count std')
    t = 0
    l = 0
    l_std = []
    with open(name+'_L_for_'+str(step)+'step.txt','r') as f:
        with open(name+'_Std_for_'+str(step)+'step.txt','w') as g:
            for j in range(len(data)):
                n = j+1 
                if data[j][0] >= step*t:
                    try:
                        d = float(f.readline().strip())
                    except ValueError:
                        continue
                    
                    t += 1
                    try:
                        l_std += [sqrt((1/float(k-1))*l)]
                        g.write(str(sqrt((1/float(k-1))*l))+'\n')
#                        print("STD for " + str(t*step) + " n.u. and " + str(k) + ' intersections = ' + str(sqrt((1/float(k-1))*l)))	
                    except ZeroDivisionError:
                        print('ZeroDivisionError',str(l),str(k),str(j), data[j])
                        g.write(str(sqrt((1/float(k))*l))+'\n')
#                        print("STD for " + str(t*step) + " n.u. and " + str(k) + ' intersections = ' + str(sqrt((1/float(k))*l)))
                    l = 0
                    k = 0
                if n >= len(data):
                    break
                while data[j][0] + data[j][1] > data[n][0]:
                    l += (data[j][0] + data[j][1] - data[n][0] - d)**2
                    k += 1
                    n += 1
                    if n >= len(data):
                        break
    print('Exit')
    if prog == 'return':
        return(l_full,l100,l_std)
def sam(name='EColi_O104'):
    print('Open sam '+name)
    import re
    with open(name+'_all.sam','r') as f:
        f.readline()
        f.readline()
        f.readline()
        p = []
        data = []
        le = 0
        for i in f:
            p = i.strip().split('\t')
            le = len(p[9])
            if p[5] != '101M':
                for j in re.findall(r'[0-9]*[A-Z]',p[5]):
                    if j[-1] == 'D':
                        le += int(j[:-1])
                    elif j[-1] == 'S':
                        le -= int(j[:-1])
                    elif j[-1] == 'N':
                        le += int(j[:-1])
                    elif j[-1] == 'I':
                        le -= int(j[:-1])
            if p[3] != '0':
                data += [[int(p[3])-1, le]]
    data.sort()
    print(data[-1])
    return data

def l_reads(name, data, step=10000, prog='return'):
    print('средняя длина ридов '+name+' '+str(step))
    from math import sqrt
    print('Number of reads = ' + str(len(data)))
    #L
    print('Count reads L')
    l = 0 
    k = 0
    t = 1
    l_all
    #L100 + L
    if prog == 'return':
        for j in range(len(data)):
                if data[j][0] >= step*t:
                    t += 1
                    l_all += [l/float(k)]
                    l = 0
                    k = 0
                l += data[j][1]
                k += 1
        return l_all
    else:
        with open(name+'_L_Reads_for_'+str(step)+'step.txt','w') as f:
            for j in range(len(data)):
                if data[j][0] >= step*t:
                    print("Avg L for " + str(t*step) + " n.u. = " + str(l/float(k)))
                    t += 1
                    f.write(str(l/float(k))+'\n')
                    l = 0
                    k = 0
                l += data[j][1]
                k += 1
    print('Exit')

def plotly_graph(title,x,y,x_name,y1_name,y2_name,y1_range=[40, 60],y2_range=[0, 150]):
    print('plotly_graph')
    import plotly
    from plotly.graph_objs import Scatter, Layout
    layout = Layout(title = title,
                       xaxis = dict(rangeslider=dict(),range=[0,max(x[1])],domain = [0, 1],title = x_name),
                       yaxis = dict(range=y1_range,domain = [0, 0.45],title = y1_name),
                       yaxis2 = dict(range=y2_range,domain = [0.55,1],title = y2_name),
                       showlegend= True)
#    layout = dict(title = 'Average High and Low Temperatures in New York',
#              xaxis = dict(title = 'Month'),
#              yaxis = dict(title = 'Temperature (degrees F)'))
    plotly.offline.plot({
        "data": [Scatter(x=x[0], y=y[0], yaxis = "y1", name = y1_name,
                         line = dict(color = ('rgb(22, 96, 167)'),
                                     width = 1)),
                 Scatter(x=x[1], y=y[1], yaxis = "y2", name = y2_name)],
        "layout": layout
    })

def get_palindrome(tripl):
    compl = {'A':'T','C':'G','G':'C','T':'A'}
#    comp = {}
#    for key, value in tripl.items():    
#        comp[key] = ''.join([compl[i] for i in key[::-1]])
    comp = ''.join([compl[i] for i in tripl[::-1]])
    return comp

def sam_mistakes(name='EColi_O104',leng=5273097,step=10000):
    print('Open sam '+name)
    import re
    with open(name+'_all.sam','r') as f:
        f.readline()
        f.readline()
        f.readline()
        d = []
        i = []
        m = []
        for j in range(0,leng,step):
            d += [0]
            i += [0]
            m += [0]
        for k in f:
            p = k.strip().split('\t')
            if p[5] != '101M':
                for j in re.findall(r'[0-9]*[A-Z]',p[5]):
                    if j[-1] in ['D','N']:
                        d[int(p[3])//step] += int(j[:-1])
                    elif j[-1] in ['S','H']:
                        m[int(p[3])//step] += int(j[:-1])
                    elif j[-1] == 'I':
                        i[int(p[3])//step] += int(j[:-1])
#    cov = cov_average(win=step,step=step,prog='return')
#   for j in range(len(cov)):
#        d[j] = d[j]/cov[j]
#        i[j] = i[j]/cov[j]
#        m[j] = m[j]/cov[j]
#    return name, [d[:len(cov)],i[:len(cov)],m[:len(cov)]]
    return name, [d,i,m]

def count_zeros(name='EColi_O104',cnt=[0]):
    cov = cov_average(name=name,win=1,step=1,prog='return')
    return sum([cov.count(i) for i in cnt])

def cor_slow(name,genome,cov,step=1000,prog='return',cgcp=['CGG','CCG','CTG','CAG','CGC','TGC'],cgcm=[]):
    # Рассчет корреляции замедлителей и покрытия
    from scipy.stats.stats import pearsonr
    print('Подгружаем последовательность генома')
    print('Длина генома '+str(len(genome)))
    print('Считаем покрытие замедляющих триплетов всего в геноме')
    print('Плюс корреляция: ',cgcp)
    print('Минус корреляция: ',cgcm)
    slow = []
    l = len(cgcp[0])
    for i in range(len(genome)-l+1):
        if genome[i:i+l] in cgcp:
            slow += [1]
        elif genome[i:i+l] in cgcm:
            slow += [-1]
        else:
            slow += [0]
    print('Считаем количество замедляющих триплетов в окне')
    k = []
    for i in range(0,len(slow),step):
        try:
            k += [sum(slow[i:i+step])]
        except IndexError:
            k += [sum(slow[i:])]
    print('Подгружаем файл с покрытием и если надо уменьшаем количество данных усреднением')
#    from scripts import cov_average
#    cov = cov_average(name='EColi_O157',win=step,step=step,prog='return')
    print('Считаем корреляцию между покрытием и распределением замедляющих триплетов')
    p = pearsonr(cov,k)[0]
    print(p)
    if prog == 'save':
        with open('correlation_between_cov_k.txt','a') as f:
            f.write(name +' '+ str(step)+ str(p)+ '\n')
    return p

def cor_cov_all(name, genome,cov,step=1000,prog='return'):
    # Корреляция для всех триплетов с покрытием
    from scipy.stats.stats import pearsonr
    print('Подгружаем референс')
    seq = 'E.coli_O157_H7_Sakai'
    print('Cтроим покрытие триплетов с заданным шагом')
    k = {}
    for i in range(len(genome)-2):
        t = i//step
        if 'N' not in genome[i:i+3]:
            try:
                k[genome[i:i+3]][t] += 1
            except KeyError:
                k[genome[i:i+3]] = [0]*t + [1]
            except IndexError:
                k[genome[i:i+3]] += [0]*(t-len(k[genome[i:i+3]]))+[1]
    for i in k:
        k[i] += [0]*(t-len(k[i])+1)
    print('Маленький расчет комплементарных триплетов')
    compl = {'A':'T','C':'G','G':'C','T':'A'}
    comp = {}
    for key, value in k.items():
    #    print(key, str(value),str(tripl[''.join([compl[i] for i in key[::-1]])]))
        comp[key] = ''.join([compl[i] for i in key[::-1]])

    print('Подгружаем файл с покрытием и если надо уменьшаем количество данных усреднением')
    print('Расчитываем корреляцию каждого триплета с покрытием')
    correlation_table = {}
    for key, value in k.items():
        if comp[key] not in correlation_table:
            correlation_table[key] = [pearsonr(cov,k[key])[0],pearsonr(cov,k[comp[key]])[0]]
            print(key, str(correlation_table[key][0]), comp[key], str(correlation_table[key][1]))
    if prog == 'save':
        with open(name+'_correlation_triplets_'+str(step)+'step.csv','w') as f:
            f.write('triplet,cor,compl_triplet,cor\n')
            for key, value in correlation_table.items():
                f.write(key +','+ str("%.3f"%correlation_table[key][0])+','+ comp[key]+','+ str("%.3f"%correlation_table[key][1])+'\n')
    return correlation_table

def epi_count(genome,name,step = 1000, sites = ['CCAGG','CCTGG','GATC']):
    count_sites = {'CCAGG':[],'CCTGG':[],'GATC':[]}
    for i in range(step,len(genome),step):
        for h in sites:
            count_sites[h] += [genome[i-step:i].count(h)]
##    graph(name+'_epi_count_'+str(step),'Genome','Number of sites',['CCAGG','CCTGG','GATC'],[i for i in range(len(genome)//step)],
##          [count_sites['CCAGG'],count_sites['CCTGG'],count_sites['GATC']], typ = '-',two_axis='no',graphs=1,lines=3,show='show',log='no',figsize=(15,8),dpi=80,linewidth=1)
    graph(name+'_epi_count_'+str(step),'Genome','Number of sites',['CCAGG'],[i for i in range(len(genome)//step)],
          count_sites['CCAGG'], typ = '-',two_axis='no',graphs=1,lines=1,show='show',log='no',figsize=(15,8),dpi=80,linewidth=1)
    graph(name+'_epi_count_'+str(step),'Genome','Number of sites',['CCTGG'],[i for i in range(len(genome)//step)],
          count_sites['CCTGG'], typ = '-',two_axis='no',graphs=1,lines=1,show='show',log='no',figsize=(15,8),dpi=80,linewidth=1)
    graph(name+'_epi_count_'+str(step),'Genome','Number of sites',['GATC'],[i for i in range(len(genome)//step)],
          count_sites['GATC'], typ = '-',two_axis='no',graphs=1,lines=1,show='show',log='no',figsize=(15,8),dpi=80,linewidth=1)
    graph(name+'_epi_count_all_'+str(step),'Genome','Number of sites',['all'],[i for i in range(len(genome)//step)],
          [count_sites['GATC'][i]+count_sites['CCTGG'][i]+count_sites['CCAGG'][i] for i in range(len(genome)//step)],
          typ = '-',two_axis='no',graphs=1,lines=1,show='show',log='no',figsize=(15,8),dpi=80,linewidth=1)
##    
    return count_sites
    
if __name__ == '__main__':
#    name, genome = ref('E.coli_0104')
#    step=1000
#    gc = gc_content('E.coli_0104', genome, step=step,save='return',prog='C+G')
#    cov = [3,10,15,36,8,4,2,8,0,13,4,56,10,21,4]
    name = 'EColi_O104'
#    cov = cov_average(win=step,step=step,prog='return')
    
#    cov = []
#    step = 600
#    with open(name+'_per_'+str(step)+'.coverage','r') as f:
#        for i in f:
#            cov += [float(i.strip().split('\t')[2])]
#    print(len(cov))
#    name, data = sam(name='EColi_O104')
#    l_full = l_another(name, data, step=step,prog='save')
#    plotly_graph('L',x=[[i*(step/1000) for i in range(len(l_full))],
#                        [i*(step/1000) for i in range(len(cov))]],y=[l_full,cov], step=step)
#    a = autocor(cov, prog='var',save='return',step=1)
#    name, dim = sam_mistakes(name=name,step=10000)
#    graph('Mistakes from SAM','Genome',['Deletions,num','Insertions,num','Mutations,num'],
#          ['D','I','M'],[i for i in range(len(dim[1]))],[dim[0]],[dim[1]],[dim[2]],graphs=3,lines=1)
#    O104_zero = count_zeros(name='EColi_O104',cnt=[0,1,2,3,4,5])
#    O157_zero = count_zeros(name='EColi_O157',cnt=[0,1,2,3,4,5])

#    name, genome = ref('E.coli_O157_H7_Sakai')
#    name = 'EColi_O157'
#    gplusc = gc_content(name,genome,win=1000,step=1000,prog='C+G',save='return')
#    cov = cov_average(name=name,win=1000,step=1000,prog='return')
#    graph(name+'_Coverage and GC-content','Genome, bp',['Coverage, num','G+C, num'],['Coverage','G+C'],
#          [i for i in range(len(cov))],[cov],[gplusc[:len(cov)]],two_axis='yes',figsize=(12,9))

#    name, genome = ref('E.coli_O157_H7_Sakai')
#    name = 'EColi_O157'
#    step=1000
#    gplusc = gc_content(name,genome,win=step,step=step,prog='C+G',save='return')
#    cov = cov_average(name=name,win=step,step=step,prog='return')
#    graph(name+'_Coverage and GC-content step'+str(step),'Genome, bp',['Coverage, num','G+C, num'],['Coverage','G+C'],
#          [i for i in range(len(cov))],[cov],[gplusc],two_axis='yes',figsize=(12,9))
#    from scipy.stats.stats import pearsonr
#    print(pearsonr(cov,gplusc)[0])

#    cov = cov_average(name=name,win=1000,step=1000,prog='return')
#    graph(name+' coverage with step 1000','Genome, bp',['Coverage, num'],['Coverage'],
#          [i for i in range(len(cov))],cov)
#    cov_with = cov_average(name=name,win=1,step=1,prog='return')
#    cov_without = cov_average(name=name,win=1,step=1,prog='return',repair='no')

#    name = 'EColi_O157'
#    step=1
#    cov5 = cov_average(name=name,win=step,step=step,prog='return')
#    z5 = count_zeros(name)
#    print(z5)
#    name = 'EColi_O104'
#    z4 = count_zeros(name)
#    cov4 = cov_average(name=name,win=step,step=step,prog='return')
#    print(z4)
    
#    x = [i for i in range(10)]
#    y1 = [i for i in range(0,-10,-1)]
#    y2 = [i for i in range(0,-20,-2)]
#    print(y1)
#    plotly_graph('Test',[x,x],[y1,y2],'x','y1','y2',y1_range=[-10,0],y2_range=[-20,0])

#    cov = cov_average(win=1000,step=1000)
#    name,data=sam()
#    print(name)
#    name='EColi_O104'
#    l_gr = l(name,data,step=1000)
#    print(len(cov),len(l_gr))
#    with open('l.csv','w') as f:
#        f.write('x,coverage,read_length\n')
#        for i in range(len(cov)):
#            f.write(','.join(i,cov[i],l_gr[0][i])+'\n')

    '''
    # Различные корреляции: с GC, G-C, с замедлителями, со всеми триплетами
    name = 'EColi_MGH107'
    step = 1000
    cov = cov_average(name=name,win=step,step=step,prog='return',repair='yes')
    n, genome = ref(name)
    gcp = gc_content(name, genome, win=step, step=step, prog = 'C+G', save='return')
    gcm = gc_content(name, genome, win=step, step=step, prog = 'C-G', save='return')
    gcmm = gc_content(name, genome, win=step, step=step, prog = 'C-G', save='return')
    from scipy.stats.stats import pearsonr
    cov_gcp = pearsonr(cov[:-1],gcp)[0]
    cov_gcm = pearsonr(cov[:-1],gcm)[0]
    cov_gcmm = pearsonr(cov[:-1],[abs(i) for i in gcmm])[0]
    print('корреляция покрытия и g+c', cov_gcp)
    print('корреляция покрытия и g-c',cov_gcm)
    print('корреляция покрытия и модуля g-c',cov_gcmm)
    p = cor_slow(name,genome,cov,step=step,prog='return',cgcp=['ACG','CGG','GGC','GCG','CCG','GCC','CGC'],
                 cgcm=['ATT','ATA','TTA','AAT','TAT','TAA','AAA'])
    cor_table = cor_cov_all(name, genome,cov,step=step,prog='save')
    '''
    '''
#   Интерактивные графики для покрытия и средней длины пересечения ридов
    name = 'EColi_O104_3'
    step = 1000
    cov = cov_average(name=name,win=step,step=step,prog='return',repair='yes')
    name, data = sam(name=name)
    name = 'EColi_O104_3'
    le = l_only(name, data, step=step, prog='return')
    print(len(cov),len(le))
    x = [i for i in range(len(le))]
    plotly_graph('Покрытие и среднее пересечение модельных ридов E.coli О104 покрытие 3, окно '+str(step),
                 [x,x],[le, cov[:len(le)]],'Genome, '+str(step)+'b.p.',
                 'L','Coverage',y1_range=[35,55],y2_range=[0, 13])
    with open('EColi_O104_3_'+str(step)+'.csv','w') as f:
        f.write('x,coverage,l\n')
        for i in range(x[-1]):
            f.write(str(i)+','+str(cov[i])+','+str(le[i])+'\n')
            
    '''
    '''
# 1  
# График зависимости покрытия от GC-контента
    name = 'EColi_O157_PCRfree'
    step = 1000
    print('Set name = '+name+'\nSet step = '+str(step))
    print('Downloading and repairing of coverage data...')
    cov = cov_average(name=name,win=step,step=step,prog='return',repair='yes')
    print('Downloading of reference genome...')
    genome = ref('EColi_O157_EDL933')
    print('Counting G+C data...')
    gcp = gc_content(name, genome, win=step, step=step, prog = 'C+G', save='return')
    from scipy.stats.stats import pearsonr
    print('Counting correlation between coverage and GC-content...')
    cov_gcp = pearsonr(cov[:-1],gcp)[0]
    print(cov_gcp)
    print('Show data at a graph')
    graph(name+'_'+str(step),'GC-content',['Coverage'],['Coverage of GC'],gcp,cov[:-1], typ = '.',
          two_axis='no',graphs=1,lines=1,show='show',log='no',figsize=(15,8),
          dpi=80,linewidth=1,color=['green','purple','orange','blue'])
    print('Writing  coverage/gc-content data into a file...')
    with open(name+'_gc_cov_'+str(step)+'.csv','w') as f:
        f.write('cov,gc-cont\n')
        for i in range(len(gcp)):
            f.write(str(cov[i])+','+str(gcp[i])+'\n')

# Линейная регрессия графика зависимости покрытия от GC-контента
    from sklearn import linear_model
    print('Build a linear regression model of coverage data')
    clf = linear_model.LinearRegression()
    gcp = [[i] for i in gcp]
    clf.fit(gcp,cov[:-1])
    print(clf.coef_)
    print('Predicting a coverage data on the base of model')
    y = []
    y = clf.predict(gcp)
    print('Write predicting data into a file')
    with open(name+'_predict_'+str(step)+'.txt', 'w') as f:
        for i in y:
            f.write(str(i)+'\n')
    print('Plot a graph of predicting data')
    from matplotlib import pyplot as plt
#    plt.plot(gcp, cov[:-1], color='b',label='real')
#    plt.plot(gcp,y, color='r',label='predict')


    x = [i for i in range(len(gcp))]
#    graph(name+'_'+str(step)+'_predict','GC-content',['Coverage'],['Real','Predicted'],gcp,[cov[:-1],y], typ = '.',
#          two_axis='no',graphs=1,lines=2,show='save',log='no',figsize=(15,8),
#          dpi=80,linewidth=1,color=['blue','green','orange','purple'])
    plt.plot(x,cov[:-1],color='b',label='real_cov')
    plt.plot(x,y,color='r',label='predict_cov')
#    plt.plot(x,[cov[i]/y[i] for i in range(len(y))])
    plt.legend(loc=2)
    plt.show()

    # Модделирование ридов на основе полученной модели GC-покрытия
    from virt_reads import do_reads_with_ready_cov
    do_reads_with_ready_cov(name='EColi_O157_EDL933', gc = y, r_len=100,error=0.05, paired = True, gap=210, step=step)
    
    '''

    #2 Кореляция двух покрытий
    name1 = 'EColi_O157_PCRfree'
    name2 = 'EColi_O157_virt3'
    step=1000
    cov1 = cov_average(name=name1,win=step,step=step,prog='return',repair='yes')
    cov2 = cov_average(name=name2,win=step,step=step,prog='return',repair='yes')
    from scipy.stats.stats import pearsonr
    cor = pearsonr(cov1,cov2)
    print(cor)
    
    
    
    '''
    name = 'E.coli_0104'
    genome = ref(name)
    step = 100000
    sites = ['GATC']
#    count_sites = epi_count(genome,name,step = step, sites = sites)
    name = 'EColi_O104_3'
    cov = cov_average(name=name,win=step,step=step)
    cor = cor_slow(name,genome,cov,step=step,prog='return',cgcp=sites)
    '''

    
