import numpy as np

for k in ['3e2h','4e2h','4e1h','4h']:
    for l in ['b','g','c']:
        for t in range(1,4,1):
            for j in ['r1','r2','r3','r4','r5']:
                try :
                    f = open("./result/" + str(k) + "_" + str(l) + str(t) + "_" + j + "_cn.txt")
                    data=[]
                    data = f.readlines()[0:]
                    n_ca = int(data[0])
                    lines = data[1:]
                    nf = len(lines[0::n_ca])
                    ic=0
                    for i in range(n_ca):
                        ic+=len(lines[0:n_ca][i].split())-2
                    nc = [ic]
                    for frame in range(0,nf):
                        nc.append(ic)
                        for index in range(n_ca):
                            for i in lines[index::n_ca][0].split():
                                if i in lines[index::n_ca][frame].split():
                                    pass
                                else: nc[frame]-=1
                    native_contact=[]
                    for i in range(nf):
                        native_contact.append([i*25, float(nc[i])/float(nc[0])])
                    native_contact=np.array(native_contact)
                    with open("./result/" + str(k) + "_" + str(l) + str(t) + "_" + j + "-cn.txt", 'w') as f:
                        for i in range(len(native_contact)):
                            f.write("%s %s \n" %(native_contact[i][0], native_contact[i][1]))
                except:
                    pass
