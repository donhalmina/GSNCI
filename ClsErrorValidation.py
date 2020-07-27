# -*- coding: utf-8 -*-
"""
Created on Thu Dec 7 11:08:55 2017

@author: Darlington Mensah
"""
import pandas as pd
import numpy as np
from collections import defaultdict
import time
import ClsSemivariogram as sv
import ClsKriging as kg
from ClsModel import Model as fnc
from pylab import plt
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import warnings


class ErrorValidation:
    """
    ErrorValidation class for performing the interpolation error validation.
    INPUT VALUE:
        semi_filepath: Filepath of file containing data (if elevation, add boundary)
        krige_filepath: Filepath of file containing data (if elevation, DO NOT add boundary) 
        prediction_filepath: Filepath of file containing interpolation grid
    OUTPUT VALUE:
        Returns multiple excel file such as minimum distance and error per each blanking radii file,
        prediction data on grid and useful multiple plot. 
    """

    def __init__(self):
        """
        Instantiation of the class
        """        
        pass

    def blanking(self, semi_data, int_data, nug, R):
    # def blanking(self, semi_data, int_data, prediction, nug, R):
        # =============================================================================
        # Defining and Initializing Variables
        # =============================================================================
        predictor = defaultdict(list)
        mindistance = defaultdict(list)
        color = defaultdict(list)
        frequency = []
        mean = []
        std = []
        lag = []
        classed_resi = []
        classed_dist = []
        col2del = []     
        
        covariance = sv.Semivariogram(semi_data).isotropy(nug)
        # input("Press Enter string to continue")

        # =============================================================================
        # 150m = max. radius for Elevation (Z) in 2001
        # 350m = max. radius for Elevation (Z) in 2011
        # 250m = max. radius for Synthetic "Profile A" Ice thickness and Profile A&B
        # 350m = max. radius for Synthetic "Profile B" Ice thickness
        # =============================================================================
        
        # =============================================================================
        # R = maximum radius between the two closest data points in the dataset.
        # It is obtained by direct observation of the distribution of the dataset.
        # =============================================================================
#        R = 350
        
        C = 10
        r0 = R / 100
        sep = np.linspace(R / C, R, C)
        blanking_radii = np.hstack((0, r0, sep))
        # blanking_radii = np.hstack((0, r0))
        
#        # =============================================================================
#        # Estimating the interpolated value for the entire grid of points
#        # =============================================================================
#        inter = []
#        dist_min = []
#        for k in range(len(prediction)):
#            inter.append(kg.Kriging().ordinary(covariance, semi_data, prediction[k, :2], 0, 0))
#            dist_min.append(np.min(cdist(int_data[:, :2], prediction[k, :2][None])))
#            print(str(k) + ' ' + str(len(prediction)))
#        inter = np.hstack(inter)
#        dist_min = np.hstack(dist_min)
#        krige_mindist = pd.DataFrame(np.column_stack((prediction, inter, dist_min)))
#        w_krige_mindist = pd.ExcelWriter(str(time.localtime()[0]) + str(time.localtime()[1]) + str(time.localtime()[2])+'_' + 
#                    str(time.localtime()[3]) + str(time.localtime()[4]) + str(time.localtime()[5]) + 
#                    '_Prediction.xlsx', engine='xlsxwriter')
#        krige_mindist.to_excel(w_krige_mindist, sheet_name='Prediction')
#        w_krige_mindist.save()
        
        # =============================================================================
        # Blanking data inside defined radius prior to kriging to obtain interpolation 
        # with its error
        # =============================================================================
        for i in range(len(blanking_radii)):
            krige = []
            min_dist = []
            for j in range(len(int_data)):
                unblanked = semi_data[((semi_data[:, :2] - int_data[j, :2])**2).sum(1) > blanking_radii[i]**2]
                min_dist.append(np.min(cdist(unblanked[:, :2], int_data[j, :2][None])))
                krige.append(np.hstack(kg.Kriging().ordinary(covariance, unblanked, int_data[j, :2])))
                print(str(i) + ' ' + str(j))
             
            predictor[i] = np.hstack(krige)
            mindistance[i] = np.hstack(min_dist)
        del krige
        del min_dist
        
        mindistance = pd.DataFrame(mindistance)
        predictor = pd.DataFrame(predictor)
        color = pd.DataFrame(predictor)
        
        predictor = predictor.T.drop_duplicates().T
        mindistance = mindistance.T.drop_duplicates().T
        color = color.T.drop_duplicates().T

        predictor = predictor.apply(pd.Series.drop_duplicates, axis=1)
        mindistance = mindistance.apply(pd.Series.drop_duplicates, axis=1)
        color = color.apply(pd.Series.drop_duplicates, axis=1)
                
        # =============================================================================
        # Get the interpolator error
        # =============================================================================
        error = (np.array(predictor).transpose() - int_data[:, 2].transpose()).transpose()
 
        # =============================================================================
        # Scatter plot of minimum distance between points versus its interpolator error
        # =============================================================================       
        mindistance = np.array(mindistance)
        color = np.array(color)
        
        # color = np.random.rand(len(predictor.columns))
        for i in range(len(predictor.columns)):
            for k in range(len(predictor)):
                color[k, i] = blanking_radii[i]

        plt.scatter(mindistance, error, c=color)
        plt.xlim(0, 450)
        plt.ylim(-250, 300)
        plt.savefig('Scatter Plot.pdf', fmt='pdf', dpi=200)
        plt.show()

        for i in range(len(predictor.columns)):
            plt.scatter(mindistance[:, i], error[:, i])
            plt.xlim(0, 450)
            plt.ylim(-250, 300)
            plt.savefig('Scatter-'+str(blanking_radii[i])+'.pdf', fmt='pdf', dpi=200)
            plt.show()
        
        vecresi = np.array(error).ravel()
        vecdist= np.array(mindistance).ravel()
        sep = np.linspace(R / C, R, C)
        lags = (sep[1:] + sep[:-1]) / 2
        lags = np.hstack((0, r0, R / (2 * C), lags, 2*lags[-1]-lags[-2]))
        
        count = -1
        for ilag in lags[:-1]:
            count = count + 1
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                frequency.append(np.sum((vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])))
                classed_resi.append(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])])
                classed_dist.append(np.average(vecdist[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
                mean.append(np.average(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
                std.append(np.std(vecresi[(vecdist[:] >= ilag) & (vecdist[:] < lags[count + 1])]))
        classed_error = pd.DataFrame(classed_resi).transpose()
        
        lag = np.hstack((0, r0, sep))
        iclassed_error = pd.DataFrame(classed_resi).transpose()
        count = -1
        for i in range(len(classed_error.columns)):
            count = count + 1
            if np.count_nonzero(~np.isnan(np.array(classed_error)[:, i])) < 100:
                col2del.append(count)
                iclassed_error.drop(i, axis=1, inplace=True)
        lag = np.delete(lag, col2del)
        mean = np.delete(mean, col2del)
        std = np.delete(std, col2del)
        classed_dist = np.delete(classed_dist, col2del)
        frequency = np.delete(frequency, col2del)
        
        # =============================================================================
        # Export Interpolator Error grouped in classes as excel file
        # =============================================================================
        error = pd.DataFrame(error)
        w_err_mindist = pd.ExcelWriter(str(time.localtime()[0]) + str(time.localtime()[1]) + str(time.localtime()[2])+'_' + 
                    str(time.localtime()[3]) + str(time.localtime()[4]) + str(time.localtime()[5]) + 
                    '_Error.xlsx', engine='xlsxwriter')
        error.to_excel(w_err_mindist, sheet_name='Error')
        
        # =============================================================================
        # Export Minimum Distance grouped in classes as excel file
        # =============================================================================
        mindistance = pd.DataFrame(mindistance)
        mindistance.to_excel(w_err_mindist, sheet_name='Minimum Distance')
        w_err_mindist.save()
        
        # =============================================================================
        # Scatter plot of the error gropued in classes defined by the vector "lag"
        # =============================================================================
        iclassed_error = np.array(iclassed_error)
        plt.scatter(np.tile(classed_dist, len(iclassed_error)), iclassed_error.flatten())
        plt.plot(classed_dist, mean, 'o', c='k')
        plt.plot(classed_dist, std, 'o', c='r')
        plt.xlabel('Distance')
        plt.ylabel('Error')
        plt.title('DBF and DEF')
#        plt.savefig(str(time.localtime()[0]) + str(time.localtime()[1]) + str(time.localtime()[2])+'_' + 
#                    str(time.localtime()[3]) + str(time.localtime()[4]) + str(time.localtime()[5]) + 
#                    '_Validation.png', fmt='png', dpi=200)
        plt.show()
        
        # weight = np.sqrt(std) ** 2
        weight = ((frequency * classed_dist) / np.sum(frequency * classed_dist))
        # weight = np.sqrt(blanking_radii) ** 2
        
        # =============================================================================
        # Least square fit function to obtain coefficient of the indeterminate
        # =============================================================================
        paramsMean = curve_fit(fnc.fitfunction, classed_dist, mean, sigma=weight, bounds=((-np.inf, -np.inf, -np.inf, 0),
                                                            (np.inf, np.inf, np.inf, 0.000001)))
        paramsStD = curve_fit(fnc.fitfunction, classed_dist, std, sigma=weight, bounds=((-np.inf, -np.inf, -np.inf, 0),
                                                          (np.inf, np.inf, np.inf, nug)))
        
        [m1, m2, m3, m4] = paramsMean[0]
        [s1, s2, s3, s4] = paramsStD[0]
        
        classed_dist = np.hstack((0, classed_dist))
        mean = np.hstack((0, mean))
        std = np.hstack((s4, std))
        frequency = np.hstack((0, frequency))
        x_int = np.linspace(np.min(classed_dist), np.max(classed_dist), 200)
        f_mean = np.poly1d(paramsMean[0])
        f_std = np.poly1d(paramsStD[0])
        mean_int = f_mean(x_int)
        std_int = f_std(x_int)
        _, ax = plt.subplots()
        plt.plot(x_int, mean_int)
        plt.plot(x_int, std_int)
#        plt.plot(classed_dist, mean, 'o', c='k', label=str(m1) + "x^3 " + str(m2) + "x^2 ")# + str(m3) + "x")
        plt.plot(classed_dist, mean, 'o', c='k', label=str(m1) + "x^3 " + str(m2) + "x^2 " + str(m3) + "x")
        plt.plot(classed_dist, std, 'x', c='r', label=str(s1) + "x^3 " + str(s2) + "x^2 " + str(s3) + "x " + str(s4))
        for i, txt in enumerate(frequency):
            ax.annotate(txt, (classed_dist[i], mean[i]))
            ax.annotate(txt, (classed_dist[i], std[i]))
        plt.legend(loc=2, fontsize='xx-small', borderaxespad=0.)
        plt.savefig('Validation Fit.pdf', fmt='pdf', dpi=200)
        plt.show()
