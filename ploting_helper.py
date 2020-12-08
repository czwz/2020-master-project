import seaborn as sns
import matplotlib.pyplot as plt

def plot_2d_density(df, ptype, param_x, param_y, level):
    """
    plot the density scatter plot based on every 10 steps (20 ns)
    Provided parameter pair:
    1. rgyr- pprmsd
    2. eprmsd- pprmsd
    3. pnc- pprmsd
    4. sermsd- sprmsd
    -----
    input(pd.DataFrame, str, str, int)
    output(plt.plot)
    """
    x = df[df.index.str.contains(ptype+"_b", regex=True)][param_x].values
    y = df[df.index.str.contains(ptype+"_b", regex=True)][param_y].values
    sns.kdeplot(x[0:-1:10],y[0:-1:10], cmap="Reds", levels=level)

    x = df[df.index.str.contains(ptype+"_g", regex=True)][param_x].values
    y = df[df.index.str.contains(ptype+"_g", regex=True)][param_y].values
    sns.kdeplot(x[0:-1:10],y[0:-1:10], cmap="Blues", levels=level)

    if ptype == "4H" or ptype == "4E1H":
        x = df[df.index.str.contains(ptype+"_c", regex=True)][param_x].values
        y = df[df.index.str.contains(ptype+"_c", regex=True)][param_y].values
        sns.kdeplot(x[0:-1:10],y[0:-1:10], cmap="Greens", levels=level)

    if param_x == "rgyr" and param_y == "pprmsd":
        if ptype == "4E1H":
            plt.ylim(0.5, 5.5)
            plt.xlim(0.95, 1.05)
        elif ptype == "4H":    
            plt.ylim(1, 10)
            plt.xlim(0.94, 1.1)
        elif ptype == "3E2H":
            plt.ylim(1, 7)
            plt.xlim(0.95, 1.2)
        elif ptype == "4E2H":
            plt.ylim(1, 7)
            plt.xlim(0.95, 1.15)
    elif param_x == "eprmsd" and param_y == "pprmsd":
        if ptype == "4E1H":
            plt.ylim(0.5, 5)
            plt.xlim(0.3, 3)
        elif ptype == "4H":
            plt.ylim(2, 10)
            plt.xlim(1.5, 10)
        elif ptype == "3E2H":
            plt.ylim(1, 6.5)
            plt.xlim(0, 4.5)
        elif ptype == "4E2H":
            plt.ylim(1.5, 7)
            plt.xlim(0, 12)            
    elif param_x == "pnc" and param_y == "pprmsd":
        if ptype == "4E1H":
            plt.ylim(0, 6)
            plt.xlim(0.68, 0.81)
        elif ptype == "4H":    
            plt.ylim(2, 10)
            plt.xlim(0.68, 0.81)
        elif ptype == "3E2H":
            plt.ylim(2, 8)
            plt.xlim(0.68, 0.81)
        elif ptype == "4E2H":
            plt.ylim(1, 7)
            plt.xlim(0.68, 0.81)
    elif param_x == 'sermsd' and param_y == "sprmsd":
        if ptype == "4E1H":
            plt.ylim(1, 7)
            plt.xlim(0, 5)
        elif ptype == "4H":    
            plt.ylim(2, 10)
            plt.xlim(1, 8)
        elif ptype == "3E2H":
            plt.ylim(2.5, 10)
            plt.xlim(0, 4)
        elif ptype == "4E2H":
            plt.ylim(2, 8)
            plt.xlim(1, 7)
    else:
        raise Exception("Error: The chosen parameter pair is not provided (" + param_x + " & " + param_y + ")")
            
    plt.title(ptype, fontsize=15)
    plt.xlabel(param_x, fontsize=15)
    plt.ylabel(param_y, fontsize=15)
    
        
