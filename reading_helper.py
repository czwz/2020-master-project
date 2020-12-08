import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class toolbox(object):
    
    def __init__(self):
        
        self.samples = ['1', '2', '3', '4','5']
        self.types = ['3E2H','4E1H','4E2H','4H']
        self.labels = ['b', 'g', 'c']
        self.mutations = ['1', '2', '3', '4']
        self.panel_data = None
        self.file_path = "./processed_output/"
    
    def create_time_series_data(self,ptype):
        """
        Load the output data from post-processing.
        FilePath Default: ./processed_output/ptype.xxx
        Make sure the naming of file is in the fashion: r_[sameple]_[type]_[label]_[metation].[rgyr/prmsd/ermsd/rmsf/nc].
        The output contains 8 time series: rgyr, pprmsd, eprmsd, sprmsd, eermsd, sermsd, pnc, enc  
        -----
        input(str)
        output(pd.DataFrame object)  
        """
        if os.path.isfile( self.file_path + ptype + '.rgyr'):
            
            # defined path for files with extension of nc and non-nc
            rec1 = glob.glob(self.file_path + ptype + ".*[!nc]", recursive=True)
            rec2 = glob.glob(self.file_path + ptype + ".*[nc]", recursive=True)

            # read the files
            df1 = (pd.read_csv(f, sep=" ", index_col=False, header=None).iloc[:,1:] for f in rec1)
            df2 = (pd.read_csv(f, sep=" ", index_col=False, header=None).iloc[:,1:] for f in rec2)
            f1 = pd.concat(df1, axis=1)
            f2 = pd.concat(df2, axis=1)
            f1.index = f1.index * 0.2
            f2.index = f2.index * 2

            # create dataframe
            big_dataframe = pd.concat([f1, f2], axis=1)
            big_dataframe = big_dataframe.dropna()
            big_dataframe.columns = ["eermsd", "eprmsd", "pprmsd", "rgyr", "sermsd", "sprmsd", "enc", "pnc"]
            big_dataframe["rgyr"] = big_dataframe["rgyr"] / big_dataframe["rgyr"][0]
            big_dataframe["time"] = big_dataframe.index
            big_dataframe.index = [ptype] * len(big_dataframe)
            big_dataframe.index.name = "type"

            return big_dataframe

    def conver_to_panel_data(self, filename):
        """
        read all the trj data into one panel data and save it as filename
        -----
        input(str)
        """
        df = pd.DataFrame({})

        for s in self.samples:
            for t in self.types:
                for l in self.labels:
                    for m in self.mutations:
                        ptype = 'r_' + s + '_' + t + '_' + l + '_' + m
                        df = pd.concat([df, self.create_time_series_data(ptype)])
        
        df.to_csv(filename)
        
    def conver_to_average_panel_data(self, filename, samplingFrequency):
        """
        precondition: self.panel_data has been loaded by self.read_panel_data
        read all the trj data into one panel data, avergae over every defined interval, and save it as filename
        -----
        input(str, int)
        """
        
        if self.panel_data is None:
            raise Exception("Error: no file been read yet.")

        output = pd.DataFrame({})

        for s in self.samples:
            for t in self.types:
                for l in self.labels:
                    for m in self.mutations:
                        ptype = 'r_' + s + '_' + t + '_' + l + '_' + m

                        ptype_df = self.panel_data[self.panel_data.index.str.contains(ptype)]
                        ptype_df = ptype_df.groupby(ptype_df["time"] // samplingFrequency).mean()
                        ptype_df["type"] = [ptype]*len(ptype_df)

                        if l == 'b':
                            ptype_df["label"] = [0]*len(ptype_df)
                        else:
                            ptype_df["label"] = [1]*len(ptype_df)

                        output = pd.concat([output, ptype_df], axis=0)

        output.to_csv(filename, index=False)       

                    
    def read_panel_data(self, path):
        """
        read the existed panel data file at path to attribute self.panel_data
        -----
        input(str)
        """
        if os.path.isfile(path):
            self.panel_data = pd.read_csv(path, index_col=0)
        else:
            raise Exception("Error: File not found.")
            
    def read_average_panel_data(self, path):
        """
        read the existed average panel data file at path to attribute self.panel_data
        -----
        input(str)
        """
        if os.path.isfile(path):
            self.average_panel_data = pd.read_csv(path, index_col=None)
        else:
            raise Exception("Error: File not found.")