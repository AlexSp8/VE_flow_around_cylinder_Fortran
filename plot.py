import pandas as pd
from matplotlib import pyplot as plt
import os
import numpy as np

working_dir = os.getcwd()


class Plot(object):

    def __init__(self, filenames: list, title: str, plot_directory_name="Graphs", plot_filename="plot_data"):

        self.filenames = filenames
        self.title = title
        self.dir_name = plot_directory_name
        self.plt_name = plot_filename
        self.plot_dir = os.path.join(working_dir, self.dir_name)

    ###############################################################
    def make_plot_directory(self):

        if not os.path.exists(self.plot_dir):
            os.mkdir(
                self.plot_dir
            )
            print("Folder generated")

        return None

    ###############################################################
    def plot(self):
        self.make_plot_directory()
        self.set_general_plot_confs()

        title = self.title
        plot_filename = self.plt_name

        for filename in self.filenames:

            if not os.path.isfile(filename):
                # print(f"File {filename} not found")
                return

            print(filename)

            if filename.__contains__("Fd"):
                self.plot_Drag(filename)
            elif filename.__contains__("Ap_"):
                self.plot_Ap(filename)
            elif filename.__contains__("Line"):
                self.plot_Line_period(filename, plot_filename)
            elif filename.__contains__("Txx"):
                self.plot_Txx(filename)
            elif filename.__contains__("Eigenvalues/Leading"):
                self.plot_Leading_eigenvalue(filename)
            elif filename.__contains__("Eigenvalues/Multi"):
                self.plot_Multi_Continuation_data(filename)
            elif filename.__contains__("Eigenvalues/"):
                self.plot_Complex_plane(filename)
            else:
                print(f"No plot function found for filename: {filename}")

        self.finalize()

    ###############################################################
    def finalize(self):
        # plt.show()
        plt.savefig(f"{self.plot_dir}/{self.plt_name}.png", dpi=300)
        # plt.savefig(f"{self.plot_dir}/{self.plt_name}.svg", dpi=1000, format="svg")
        plt.close()

    ###############################################################
    def set_general_plot_confs(self):
        plt.title(self.title)

    ###############################################################
    def plot_Drag(self, filename):

        if filename.__contains__(".dat"):
            col_names = ['var', 'Fd']
            df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
            plt.plot(df['var'].to_numpy()[:, None], df['Fd'].to_numpy()[:, None], 
                color='r', label="Numerical", linewidth=1)
        elif filename.__contains__("data/"):
            col_names = ['var', 'Fd']
            df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
            plt.scatter(df['var'], df['Fd'], label="Experiment", linewidth=1)
        else:
            print("Wrong data for Drag")

        # plt.yscale("log")
        # plt.xscale("log")
        plt.legend()

    ###############################################################
    def plot_Ap(self, filename):

        col_names = ['var', 'Ap_2D', 'Ap_3D']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        plt.plot(df['var'].to_numpy()[:, None], df['Ap_2D'].to_numpy()[:, None], 
            color='b', label="Ap_2D", linewidth=1)
        plt.plot(df['var'].to_numpy()[:, None], df['Ap_3D'].to_numpy()[:, None], 
            color='r', label="Ap_3D", linewidth=1)

        plt.legend()

    ###############################################################
    def plot_Line_period(self, filename, plot_filename):

        col_names = ['z', 'Ux', 'Uy', 'Uz', 'P', 'Txx', 'Txy', 'Tyy', 'Txz', 'Tyz', 'Tzz']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        variable = plot_filename.split()[0]
        # print(variable)
        plt.plot(df['z'].to_numpy()[:, None], df[variable].to_numpy()[:, None], 
            color='r', label=variable, linewidth=1)

        plt.legend()

    ###############################################################
    def plot_Txx(self, filename):

        col_names = ['s', 'Txx']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        plt.plot(df['s'].to_numpy()[:, None], df['Txx'].to_numpy()[:, None], 
            color='r', label='Txx', linewidth=1)

        plt.legend()

    ###############################################################
    def plot_Leading_eigenvalue(self, filename):

        col_names = ['var1', 'Re(l)', 'Im(l)', 'target']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        plt.plot(df['var1'].to_numpy()[:, None], df['Re(l)'].to_numpy()[:, None], 
            color='r', label='Numerical', linewidth=1)

        plt.legend()

    ###############################################################
    def plot_Complex_plane(self, filename):

        col_names = ['Re(l)', 'Im(l)']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        plt.scatter(df['Re(l)'].to_numpy()[:, None], df['Im(l)'].to_numpy()[:, None], 
            color='r', label='Numerical', linewidth=1)

        plt.legend()

    ###############################################################
    def plot_Multi_Continuation_data(self, filename):

        col_names = ['var2', 'Q', 'var1','Re(l)', 'Im(l)', 'var1c']
        df = pd.read_csv(filename, sep=r'\s+', names=col_names, header=None)
        plt.scatter(df['Q'].to_numpy()[:, None], df['var1c'].to_numpy()[:, None], 
            color='r', label='Numerical', linewidth=1)

        plt.legend()

###############################################################
data_folder = 'data/'
root_folder = 'Results/'
Problem_Main = 'Cylinder_3D'
Problem_Stability = 'Cylinder_Stability_3D'
Cvar1_name = "WiN"
# Cvar2_name = ["BvN", "VeN", "ReN"]
Cvar2_name = ["Dum"]

variables = ["Uz"]
# variables = ["Ux", "Uy", "Uz", "P", "Txx", "Txy", "Tyy", "Txz", "Tyz", "Tzz"]

Nshifts = 1
Neigen = 50

###############################################################
def setVariableVals(var0,varf,dvar):
    Nsteps = round(1+(varf-var0)/(dvar+1e-8))
    vals = []
    for i in range(0,Nsteps):
        val = var0+i*dvar
        vals.append(val)
    return vals

dBvN = 0.05
BvN0 = 0.1
BvNf = 0.8
BvN_vals = setVariableVals(BvN0,BvNf,dBvN)

dVeN = 0.005
VeN0 = 0.01
VeNf = 0.03
VeN_vals = setVariableVals(VeN0,VeNf,dVeN)

dReN = 0.2
ReN0 = 0.0
ReNf = 5.0
ReN_vals = setVariableVals(ReN0,ReNf,dReN)

dDum = 0.0
Dum0 = 0.1
Dumf = 0.1
Dum_vals = setVariableVals(Dum0,Dumf,dDum)

dQbN = 0.05
QbN0 = 0.0
QbNf = 0.5
QbN_vals = setVariableVals(QbN0,QbNf,dQbN)

dcont = 0.1
cont0 = 3.5
contf = 3.5
Cvar1_vals = setVariableVals(cont0,contf,dcont)

###############################################################
plotter = Plot(
    filenames=[],
    plot_directory_name='Graphs/',
    plot_filename="",
    title=""
)
plotter.make_plot_directory()
del plotter
###############################################################
#Base
for i in range(0,len(Cvar2_name)):

    name2 = Cvar2_name[i]

    Cvar2_vals = []

    if (name2=="BvN"):
        Cvar2_vals = BvN_vals
    elif (name2=="VeN"):
        Cvar2_vals = VeN_vals
    elif (name2=="ReN"):
        Cvar2_vals = ReN_vals
    elif (name2=="Dum"):
        Cvar2_vals = Dum_vals

    folder1 = root_folder+"Base/"+name2+"/"+name2

    for j in range(0,len(Cvar2_vals)):

        val2 = round(Cvar2_vals[j],4)

        folder = folder1+"_{val_2:.4f}/"

        #Fd
        file = folder+"Flow_data/Fd_"+Problem_Main+".dat"
        files = [file.format(val_2 = val2)]
        plotter = Plot(
            filenames=files,
            plot_directory_name='Graphs/Base',
            plot_filename="Fd_"+Problem_Main,
            title="Fd_"+Problem_Main
        )
        plotter.plot()
        del plotter

        #Ap
        file = folder+"Flow_data/Ap_"+Problem_Main+".dat"
        files = [file.format(val_2 = val2)]
        plotter = Plot(
            filenames=files,
            plot_directory_name='Graphs/Base',
            plot_filename="Ap_"+Problem_Main,
            title="Ap_"+Problem_Main
        )
        plotter.plot()
        del plotter
            
        #Line
        for k in range(0,len(Cvar1_vals)):

            val1 = round(Cvar1_vals[k],4)

            file = folder+"Line/"+Problem_Main+"_"+Cvar1_name+"_{val_1:.4f}.dat"
            files = [file.format(val_2 = val2, val_1 = val1)]
            for v in range(0,len(variables)):
                plotter = Plot(
                    filenames=files,
                    plot_directory_name='Graphs/Base',
                    plot_filename=(variables[v]+" period: "+name2+" = "
                        +str(val2)+", "+Cvar1_name+" = "+str(val1)),
                    title=(variables[v]+" period: "+name2+" = "
                        +str(val2)+", "+Cvar1_name+" = "+str(val1)),
                )
                plotter.plot()
                del plotter

            #Txx
            file = folder+"Txx/"+Problem_Main+"_"+Cvar1_name+"_{val_1:.4f}.dat"
            files = [file.format(val_2 = val2, val_1 = val1)]
            for v in range(0,len(variables)):
                plotter = Plot(
                    filenames=files,
                    plot_directory_name='Graphs/Base',
                    plot_filename=("Txx: "+name2+" = "
                        +str(val2)+", "+Cvar1_name+" = "+str(val1)),
                    title=("Txx "+name2+" = "
                        +str(val2)+", "+Cvar1_name+" = "+str(val1)),
                )
                plotter.plot()
                del plotter

###############################################################
#Stability
files = []
for i in range(0,len(Cvar2_name)):

    name2 = Cvar2_name[i]

    Cvar2_vals = []

    if (name2=="BvN"):
        Cvar2_vals = BvN_vals
    elif (name2=="VeN"):
        Cvar2_vals = VeN_vals
    elif (name2=="ReN"):
        Cvar2_vals = ReN_vals
    elif (name2=="Dum"):
        Cvar2_vals = Dum_vals

    folder1 = root_folder+"Stability/"+name2+"/"+name2

    for j in range(0,len(Cvar2_vals)):

        val2 = round(Cvar2_vals[j],4)

        folder2 = folder1+"_{val_2:.4f}/"

        for k in range(0,len(QbN_vals)):

            Q = round(QbN_vals[k],4)

            folder = folder2+"QbN_{Q_val:.4f}/"

            #Leading Eigenvalues
            file = folder+"Eigenvalues/Leading_eigenvalues.dat"
            files = [file.format(val_2 = val2, Q_val = Q)]
            plotter = Plot(
                filenames=files,
                plot_directory_name='Graphs/Stability',
                plot_filename="Leading_eigenvalues: Q = "+str(Q),
                title="Leading_eigenvalues: Q = "+str(Q)
            )
            plotter.plot()
            del plotter

            for i2 in range(0,len(Cvar1_vals)):

                val1 = round(Cvar1_vals[i2],4)

                #Complex plane
                file = folder+"Eigenvalues/"+Cvar1_name+"_{val_1:.4f}.dat"
                files = [file.format(val_2 = val2, Q_val = Q, val_1 = val1)]
                plotter = Plot(
                    filenames=files,
                    plot_directory_name='Graphs/Stability',
                    plot_filename=("Complex plane: "+name2+" = "+str(val2)
                        +", QbN = "+str(Q)+", "+Cvar1_name+" = "+str(val1)),
                    title=("Complex plane: "+name2+" = "+str(val2)
                        +", QbN = "+str(Q)+", "+Cvar1_name+" = "+str(val1)),
                )
                plotter.plot()
                del plotter

                for j2 in range(1,Neigen):

                    for ishift in (0,Nshifts):

                        ival = j2 + ishift*1000

                        #Line
                        folder3 = folder+"Line/"+Cvar1_name+"_{val_1:.4f}/"
                        file = folder3+"Kr_{ival_1:}.dat"
                        files = [file.format(val_2 = val2, Q_val = Q, val_1 = val1, ival_1 = ival)]
                        for v in range(0,len(variables)):
                            plotter = Plot(
                                filenames=files,
                                plot_directory_name='Graphs/Stability',
                                plot_filename=(variables[v]+" period stability: "
                                +name2+" = "+str(val2)+", QbN = "+str(Q)+", "
                                +Cvar1_name+" = "+str(val1)+", i = "+str(ival)),
                                title=(variables[v]+" period stability: "+name2+" = "+
                                str(val2)+", QbN = "+str(Q)+", "
                                +Cvar1_name+" = "+str(val1)+", i = "+str(ival)),
                            )
                            plotter.plot()
                            del plotter

# ###############################################################
# #Multi Continuation data
# file = root_folder+"Eigenvalues/Multi_data.dat"
# plotter = Plot(
#     filenames=[file],
#     plot_directory_name='Graphs',
#     plot_filename="Multi_data.dat",
#     title="Multi_data.dat"
# )
# plotter.plot()
# del plotter

