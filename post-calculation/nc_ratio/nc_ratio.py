import numpy as np
import matplotlib.pyplot as plt

def compute_nc_ratio(filename):
    """
    The file format for input (raw native contacts) should be 
    
    number-of-indices
    timeframe ref-index contact-index-1 contact-index-2 contact-index-3 ...
    timeframe ref-index contact-index-1 contact-index-2 contact-index-3 ...
    timeframe ref-index contact-index-1 contact-index-2 contact-index-3 ...
    ...
    
    It calculates the time evolution of overall native contact ratio R = NC(t)/NC(t0),
    where NC stands for total number of native contacts; returns a np.array of time 
    series of contact ratio.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

        data = {}
        data["num_indices"] = int(lines[0])
        data["num_frame"] = (len(lines)-1)//data["num_indices"]
        data["dt"] = int(lines[1+data["num_indices"]].split()[0]) - int(lines[1][0])
        data["contact_indices"] = []

        for line in lines:
            line_split = line.split()
            if len(line_split) > 1:
                contacts_each_index = np.array([int(i) for i in line_split[2:]])
                data["contact_indices"].append( contacts_each_index )

    initial_contacts=np.sum( [len(i) for i in data["contact_indices"][0:data["num_indices"]]] )

    nc_ratio = np.full(shape=(data["num_frame"],2), fill_value=1.)
    for frame in range(data["num_frame"]):
        nc_ratio[frame][0]=frame*data["dt"]
        for index in range(data["num_indices"]):
            for i in data["contact_indices"][index]:
                if i not in data["contact_indices"][frame*data["num_indices"]+index]:
                     nc_ratio[frame][1]-=1/initial_contacts
    return nc_ratio

def output_nc_ratio(filename, nc_ratio):
    with open(filename, 'w') as f:
        for i in range(len(nc_ratio)):
            f.write(str(nc_ratio[i][0]) + " " + str(nc_ratio[i][1]) + "\n")
            
nc_ratio = compute_nc_ratio("4h_b2_cn.txt")
output_nc_ratio("4h_b2_cn.out", nc_ratio)