import mbfit
import os
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde

from mbfit.training_set import TrainingSet
from mbfit.training_set import TrainingSetElement
from mbfit.polynomials import MoleculeSymmetryParser
from mbfit.utils import SettingsReader

class PCA_analysis:

    def __init__(self, setting_path, train_file_path, scan_file_path, sym_string):
        """
        Generate a new object of PCA_analysis for a molecule with given information.
        Args:
            setting_path            - The path of the setting.ini with all relevant settings.
            train_file_path         - The path of the training xyz file.
            scan_file_path          - The path of the testing xyz file or the scanning data file.
            sym_string              - The symmetrical string for the molecule.

        Returns:
            Nothing
        """
        scan_lst = []
        setting = SettingsReader(setting_path)
        train = TrainingSet.get_training_set_from_xyz_file(train_file_path, setting, ['weight', 'energy'])
        
        if type(scan_file_path) is str:
            scan = TrainingSet.get_training_set_from_xyz_file(scan_file_path, setting, ['weight', 'energy'])
            scan_lst.append(scan)

        else:
            for file_path in scan_file_path:
                scan = TrainingSet.get_training_set_from_xyz_file(file_path, setting, ['weight', 'energy'])
                scan_lst.append(scan)

        X_train = []
        y_train = []
        X_test = []
        y_test = []
        self.sym_string = sym_string
        self.setting = setting

        for mole in train.get_elements():
            distances = self.mole_dis(mole.get_molecule(), sym_string)
            X_train.append(distances)
            y_train.append(mole.get_energy("energy"))

        for scan in scan_lst:
            for mole in scan.get_elements():
                distances = self.mole_dis(mole.get_molecule(), sym_string)
                X_test.append(distances)
                y_test.append(mole.get_energy("energy"))
    
        pca = PCA()
        pca.fit(X_train, y_train)
    
        self.pca = pca
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test
        self.all_data = pca.transform(X_train)
        self.all_scan_data = pca.transform(X_test)

    
    def get_pca(self):
        """
        Get the sklearn PCA object.

        Documentation: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html

        Returns:
            A Sklean PCA object
        """
        return self.pca
    

    def get_train_energies(self):
        """
        Get the energies of the training set.

        Returns:
            A list with all the energies of molecules in training set.
        """
        return self.y_train
        

    def get_test_energies(self):
        """
        Get the energies of the testing set or the scanning molecules.

        Returns:
            A list with all the energies of molecules in testing set. 
        """
        return self.y_test

    
    def pca_graph(self, save, scan_color="red"):
        """
        Generate the PCA scatter plots. It will save figures files to local based on input.

        Args:
            save        - True or False determining if this method save the figures files
            scan_color  - Cutomize the color of the scanning points on the graph

        Returns:
            Nothing
        """
        y_train = self.y_train
        all_data = self.all_data

        y_test = self.y_test
        all_scan_data = self.all_scan_data

        for n in range(len(all_data[0])):
            comp = self.to_comp(all_data, n)
            s_comp = self.to_comp(all_scan_data, n)
            title = "plot for {} component".format(n + 1)
            s_title = "plot_for_{}_component.png".format(n + 1)
            print(title)
            
            if save:
                self.density_scatter_save(comp, y_train, s_comp, y_test, scan_color, s_title)
            
            else:
                self.density_scatter_display(comp, y_train, s_comp, y_test, scan_color)


    
    def get_train_comps(self):
        """
        Get the PCA components of the training set.

        Args:
            Nothing

        Returns:
            A list storing all the PCA components of the training set.
        """
        comp_lst = []
        pca = self.pca
        X_train = self.X_train
        all_data = pca.transform(X_train)

        for n in range(len(all_data[0])):
            comp = self.to_comp(all_data, n)
            comp_lst.append(comp)

        return comp_lst


    def get_test_comps(self):
        """
        Get the PCA components of the testing set.

        Args:
            Nothing
        
        Returns:
            A list storing all the PCA components of the testing set.
        """
        comp_lst = []
        pca = self.pca
        X_test = self.X_test
        all_scan_data = pca.transform(X_test)

        for n in range(len(all_scan_data[0])):
            s_comp = self.to_comp(all_scan_data, n)
            comp_lst.append(s_comp)

        return comp_lst




    def mole_dis(self, mole, string):
        """
        Helper method Computes the distances between each atoms
        """
        out = {}
        sym = MoleculeSymmetryParser(string)    
        lst = list(sym.get_variables())
        dic = {}
        for tup in lst:
            key = (tup[0], tup[3])
            value = tup[-1]
            dic[key] = value
            
        
        num_atom = mole.get_num_atoms()
        
        atom_lst = mole.get_atoms() # Get the list of the atoms
        
        
        for index_fir in range(num_atom):
            for index_sec in range(index_fir + 1, num_atom):
                #loop through all combinations of atoms
                
                atom_fir = atom_lst[index_fir]
                atom_sec = atom_lst[index_sec]
                dis = atom_fir.distance(atom_sec) #calculate the distances between two atoms
                
                sym_fir = atom_fir.get_symmetry_class()
                sym_sec = atom_sec.get_symmetry_class()
                
                if (sym_fir, sym_sec) not in dic.keys():
                    #if this combination is not existed in the dictionary from symmetric pairs, skip
                    continue
            
                sym_term = dic[(sym_fir, sym_sec)]
                
                if sym_term not in out.keys():
                    out[sym_term] = [dis]
                    
                else:
                    dis_lst = out[sym_term]
                    dis_lst.append(dis)
        
        final_out = []
        
        for sym_term in out.keys():
            dis_lst = out[sym_term]
            dis = np.mean(dis_lst)
            final_out.append(dis)
        
        return final_out


    def to_comp(self, array, n_comp):
        """
        Helper method that categorizes PCA data 

        """
        comp = []
        for lst in array:
            comp.append(lst[n_comp])
            
        return comp


    def density_scatter_display(self, x, y, scan_x, scan_y, color):
        """
        Helper method used to display a density scatter plot.
        """
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)

        fig, ax = plt.subplots()
        ax.scatter(x, y, c=z, s=50)
        ax.scatter(scan_x, scan_y, c=color, s=50)
        plt.show()

    
    def density_scatter_save(self, x, y, scan_x, scan_y, color, save_path):
        """
        Helper method used to display and save a density scatter plot.
        """
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)

        fig, ax = plt.subplots()
        ax.scatter(x, y, c=z, s=50)
        ax.scatter(scan_x, scan_y, c=color, s=50)
        plt.show()
        plt.savefig(save_path)
