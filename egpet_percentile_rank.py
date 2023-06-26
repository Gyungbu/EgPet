import os, datetime
import pandas as pd
from scipy.stats import percentileofscore, pearsonr
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy
from skbio.stats.composition import multiplicative_replacement, clr

#-------------------------------------------------------
# Common Function
#-------------------------------------------------------
def WriteLog(functionname, msg, type='INFO', fplog=None):
    #strmsg = "[%s][%s][%s] %s\n" % (datetime.datetime.now(), type, functionname, msg)
    #if( DEBUGMODE ): 
    #    print(strmsg)
    #else:
    #    if( fplog != None ):
    #        fplog.write(strmsg)
    
    head = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    writestr = f"[{head}][{functionname}] {msg}\n"
    #if( DEBUGMODE ):
    if( True ):
        #print(message)
        writestr = f"[{functionname}] {msg}\n"
        print(writestr)
        
    if( fplog != None ):
        fplog.write(writestr)
        fplog.flush()

###################################
# MainClass
###################################
class EgPetAnalysis:
    def __init__(self, path_exp, outdir=None, fplog=None):
        """
        Initializes a EgPetAnalysis object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.
        """
        self.path_exp = path_exp
        self.species = None
        self.__fplog=fplog
        if( os.path.basename(self.path_exp).startswith('PD') ):
            self.species = 'dog'
        elif( os.path.basename(self.path_exp).startswith('PC') ):
            self.species = 'cat'

        if( self.species is None ):
            print("The species should be dog[PD] or cat[PC]")
            print("Please check the path_exp")

        ## Reference csv files
        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_beta = f"{curdir}/input/phenotype_microbiome_{self.species}.csv"
        self.path_healthy = f"{curdir}/input/healthy_profile_{self.species}.csv"
        self.path_mrs_db = f"{curdir}/input/egpet_mrs_db_{self.species}.csv"
        self.path_percentile_rank_db = f"{curdir}/input/egpet_percentile_rank_db_{self.species}.csv"
        self.path_dysbiosis = f"{curdir}/input/dysbiosis_microbiome_{self.species}.csv"
        self.path_db = f"{curdir}/input/db_abundance_{self.species}.csv"
        
        ###output
        if( outdir is not None ):
            self.outdir = outdir
        else:
            self.outdir = f"{curdir}/output"

        #self.path_egpet_percentile_rank_output = f"{self.outdir}/{os.path.basename(self.path_exp).replace('_report.txt','')}.csv"
        #self.path_egpet_eval_output = f"{self.outdir}/{os.path.basename(self.path_exp).replace('_report.txt','_eval')}.csv"
        #self.path_egpet_scatterplot_output = f"{self.outdir}/{os.path.basename(self.path_exp).replace('_report.txt','_scatterplot')}.png"
        
        self.path_egpet_percentile_rank_output = f"{curdir}/output/egpet_percentile_rank_{self.species}.csv"
        self.path_egpet_eval_output = f"{curdir}/output/egpet_eval_{self.species}.csv"
        self.path_egpet_scatterplot_output = f"{curdir}/output/egpet_scatterplot_{self.species}.png"

        ## Dataframes read by the ReadDB process
        self.df_beta = None
        self.df_healthy = None
        self.df_exp = None
        self.df_mrs_db = None
        self.df_dysbiosis = None
        self.df_percentile_rank_db = None
        self.df_db = None
        
        ## Dataframes to calculate
        self.df_mrs = None
        self.df_percentile_rank = None
        self.df_eval = None
        
        ## Lists used for calculation
        self.li_diversity = None
        self.li_observed = None
        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None

    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):
        """
        Read the Dataframes.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:
            self.df_beta = pd.read_csv(self.path_beta, encoding='cp949')
            self.df_dysbiosis = pd.read_csv(self.path_dysbiosis, encoding='cp949')
            self.df_healthy = pd.read_csv(self.path_healthy, encoding='cp949')
            self.df_exp = pd.read_csv(self.path_exp, encoding='cp949')
            self.df_mrs_db = pd.read_csv(self.path_mrs_db, index_col=0, encoding='cp949') 
            self.df_exp = pd.read_csv(self.path_exp, encoding='cp949')
            self.df_percentile_rank_db = pd.read_csv(self.path_percentile_rank_db, encoding='cp949')
            self.df_db = pd.read_csv(self.path_db, encoding='cp949')

            self.df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_beta = self.df_beta[["phenotype", "ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_beta['beta'] = self.df_beta['beta'].replace({'유해': 1, '유익': -1})

            self.df_dysbiosis.rename(columns = {"NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_dysbiosis = self.df_dysbiosis[["ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_dysbiosis['beta'] = self.df_dysbiosis['beta'].replace({'유해': 1, '유익': -1})
            
            print(self.df_exp)
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(self.df_db['taxa'][0:2]) == ['diversity', 'observed']):
                self.li_diversity = list(self.df_exp.iloc[0,1:]) # li_diversity : Alpha-Diversity list 
                self.li_observed = list(self.df_exp.iloc[1,1:]) # li_observed : Number of Microbiome list
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_db = self.df_db.iloc[2:,:]
                            
            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
            print(self.df_beta)
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            
        return rv, rvmsg


    def CalculateMRS(self): 
        """
        Calculate the MRS (Microbiome Risk Score).

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:                
            # df_mrs : Data frame of MRS corresponding to specific phenotype and sample
            self.df_mrs = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            self.df_mrs = self.df_mrs.fillna(0) 

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    condition_phen = (self.df_beta.phenotype == self.li_phenotype[j])   
                    mrs = 0

                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows():
                        condition_micro = (self.df_exp.taxa == row_beta['microbiome'])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0] 
                        
                        li_micro_sub = []

                        if pd.isna(row_beta['microbiome_subtract']) is False:
                            li_micro_sub = row_beta['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_exp.taxa == micro_sub)

                                if len(self.df_exp[condition_sub]) > 0:
                                    abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]     
                            
                        mrs += row_beta['beta'] * math.log10(100*abundance + 1) 

                    mrs /= len(self.df_beta[condition_phen])       
                    self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = -mrs

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
    
        return rv, rvmsg

    def CalculateDysbiosis(self): 
        """
        Calculate the Dysbiosis.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_mrs['Dysbiosis'] = 0
            self.li_microbiome = list(dict.fromkeys(self.df_dysbiosis['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                dysbiosis = 0
                
                for j in range(len(self.li_microbiome)):
                    condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) 
                    condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) 
                    
                    if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]    
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]                                  
                                        
                            dysbiosis += math.log10(100*abundance + 1)            
                            
                    elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     
                                
                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]       
                                        
                            dysbiosis -= math.log10(100*abundance + 1)      
                            
                self.df_mrs.loc[self.li_new_sample_name[i], 'Dysbiosis'] = -dysbiosis
                         
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
    
        return rv, rvmsg
     
    def CalculateHealthyDistance(self): 
        """
        Calculate the Healthy Distance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """           
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_mrs['HealthyDistance'] = 0     
            self.df_mrs['Diversity'] = self.li_diversity
            
            np_healthy_abundance = self.df_healthy['RA'].to_numpy()
            np_healthy_abundance = np.append(np_healthy_abundance, 100-np_healthy_abundance.sum())

            np_healthy_abundance = clr(np_healthy_abundance)
            
            # Subtract the abundance - df_exp_healthy
            for idx in range(len(self.li_new_sample_name)): 
                df_exp_one = self.df_exp[['taxa', self.li_new_sample_name[idx]]]
                df_exp_one = df_exp_one[df_exp_one[self.li_new_sample_name[idx]] != 0]
                np_abundance = np.array([], dtype=np.float64).reshape(0,1)
                np_abundance_others = np.ones((1,1), dtype=float)                
                
                for idx_healthy, row_healthy in self.df_healthy.iterrows(): 
                    li_micro_sub = []
                    li_micro = row_healthy['microbiome'].split('\n')
                    np_abundance_temp = np.zeros((1,1), dtype=float)

                    for micro in li_micro:
                        condition_append = (df_exp_one.taxa == micro)

                        if len(df_exp_one[condition_append]) > 0:
                            np_abundance_temp += df_exp_one[condition_append].to_numpy()[:,1:].astype(np.float64)
                            np_abundance_others -= df_exp_one[condition_append].to_numpy()[:,1:].astype(np.float64)

                    if pd.isna(row_healthy['microbiome_subtract']) is False:
                        li_micro_sub = row_healthy['microbiome_subtract'].split('\n')

                        for micro_sub in li_micro_sub:
                            condition_sub = (df_exp_one.taxa == micro_sub)

                            if len(df_exp_one[condition_sub]) > 0:
                                np_abundance_temp -= df_exp_one[condition_sub].to_numpy()[:,1:].astype(np.float64)
                                np_abundance_others += df_exp_one[condition_sub].to_numpy()[:,1:].astype(np.float64)

                    np_abundance = np.concatenate((np_abundance,np_abundance_temp),axis=0)

                np_abundance = np.concatenate((np_abundance,np_abundance_others),axis=0)
                np_abundance = np_abundance.transpose()

                # Apply multiplicative replacement and CLR transformations
                np_abundance = multiplicative_replacement(np_abundance)
                np_abundance = clr(np_abundance)   
            
                # Calculate healthy distance for each new sample
                healthy_dist = np.linalg.norm(np_abundance - np_healthy_abundance)  
                
                self.df_mrs.loc[self.li_new_sample_name[idx], 'HealthyDistance'] = -healthy_dist
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    
    
    def CalculatePercentileRank(self):
        """
        Calculate the Percentile Rank and Save the Percentile Rank data as an csv file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:      
            # Append the Dysbiosis, HealthyDistance, Diversity, TotalRiskScore to phenotype list
            self.li_phenotype += ['Dysbiosis', 'HealthyDistance', 'Diversity']

            # Create an empty data frame with the same index and columns as the df_mrs data frame
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs_db[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                 
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95.0      
                
            self.df_percentile_rank['TotalScore'] = (self.df_percentile_rank['Dysbiosis']*1.1 + self.df_percentile_rank['HealthyDistance']*1.1 + self.df_percentile_rank['Diversity']*0.8)/3
            
            self.df_percentile_rank['TotalScore'] = self.df_percentile_rank['TotalScore'].astype(float).round(1)
                        
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')

            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_egpet_percentile_rank_output, encoding="utf-8-sig", index_label='serial_number')
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg


    def DrawScatterPlot(self):
        """
        Draw a scatter plot using the values of Diversity, Dysbiosis, and HealthyDistance

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:  
            '''
            #create regplot
            p = sns.regplot(data=self.df_percentile_rank_db, x=self.df_percentile_rank_db['Diversity'], y=(self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2)

            #calculate slope and intercept of regression equation
            slope, intercept, r, p, sterr = scipy.stats.linregress(x=p.get_lines()[0].get_xdata(),
                                                                   y=p.get_lines()[0].get_ydata())

            # generate x and y values for the line
            x_vals = np.linspace(start=self.df_percentile_rank_db['Diversity'].min(), stop=self.df_percentile_rank_db['Diversity'].max(), num=100)
            y_vals = intercept + slope * x_vals                   
            '''
            sns.scatterplot(x=self.df_percentile_rank_db['Diversity'], y=(self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2, hue = self.df_percentile_rank_db['TotalScore'] , data=self.df_percentile_rank_db)
            
            # add new points to the scatter plot
            sns.scatterplot(x=self.df_percentile_rank['Diversity'], y=(self.df_percentile_rank['Dysbiosis'] + self.df_percentile_rank['HealthyDistance'])/2, data=self.df_percentile_rank, color='g')            
            
            #plt.plot(x_vals, y_vals, '--', color='lightgray', label=f'y = {slope:.2f}x + {intercept:.2f}')
            plt.xlabel('DiversityScore')
            plt.ylabel('avg(DysbiosisScore, HealthSimilarityScore)')
            plt.legend()
                                 
            plt.axhline(y=60/1.1, xmin=0, xmax=1, color='red', linestyle='--')    
            plt.axvline(x=60/0.8, ymin=0, ymax=1, color='red', linestyle='--')
            
            '''
            E_data = self.df_percentile_rank_db[(self.df_percentile_rank_db['Diversity'] >= 60/0.8) & ((self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2 >= 60/1.1)]
            B_data = self.df_percentile_rank_db[(self.df_percentile_rank_db['Diversity'] < 60/0.8) & ((self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2 >= 60/1.1)]
            D_data = self.df_percentile_rank_db[(self.df_percentile_rank_db['Diversity'] < 60/0.8) & ((self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2 < 60/1.1)]
            I_data = self.df_percentile_rank_db[(self.df_percentile_rank_db['Diversity'] >= 60/0.8) & ((self.df_percentile_rank_db['Dysbiosis'] + self.df_percentile_rank_db['HealthyDistance'])/2 < 60/1.1)]

            E_percent = len(E_data) / len(self.df_percentile_rank_db) * 100
            B_percent = len(B_data) / len(self.df_percentile_rank_db) * 100
            D_percent = len(D_data) / len(self.df_percentile_rank_db) * 100
            I_percent = len(I_data) / len(self.df_percentile_rank_db) * 100
            
            print(f"<{self.species}>")
            print("Percentage of samples in E: ", E_percent, '%')
            print("Percentage of samples in B: ", B_percent, '%') 
            print("Percentage of samples in D: ", D_percent, '%')
            print("Percentage of samples in I: ", I_percent, '%')  
            '''
            
            # save the scatter plot
            plt.savefig(self.path_egpet_scatterplot_output , dpi=300, bbox_inches='tight')          
              
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    
    
    def EvaluatePercentileRank(self):
        """
        Evaluate based on percentile rank value and Save the Evaluation data as an csv file

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                  
            # Define the conditions and corresponding values
            conditions = [
                self.df_percentile_rank > 90,
                (self.df_percentile_rank > 70) & (self.df_percentile_rank <= 90),
                (self.df_percentile_rank > 50) & (self.df_percentile_rank <= 70),
                (self.df_percentile_rank > 30) & (self.df_percentile_rank <= 50),
                self.df_percentile_rank <= 30
            ]
            values = [1, 2, 3, 4, 5]

            # Apply the conditions and values using np.select()
            self.df_eval = pd.DataFrame(np.select(conditions, values),
                                        index=self.df_percentile_rank.index,
                                        columns=self.df_percentile_rank.columns)
            self.df_eval = self.df_eval.iloc[:, :-1]

            # Type E, B, I, D
            conditions = [
                (self.df_percentile_rank['Diversity'] >= 60/0.8) & ((self.df_percentile_rank['Dysbiosis'] + self.df_percentile_rank['HealthyDistance'])/2 >= 60/1.1),
                (self.df_percentile_rank['Diversity'] < 60/0.8) & ((self.df_percentile_rank['Dysbiosis'] + self.df_percentile_rank['HealthyDistance'])/2 >= 60/1.1),
                (self.df_percentile_rank['Diversity'] >= 60/0.8) & ((self.df_percentile_rank['Dysbiosis'] + self.df_percentile_rank['HealthyDistance'])/2 < 60/1.1),
                (self.df_percentile_rank['Diversity'] < 60/0.8) & ((self.df_percentile_rank['Dysbiosis'] + self.df_percentile_rank['HealthyDistance'])/2 < 60/1.1)
            ]
            values = ['E', 'B', 'I', 'D']

            self.df_eval['Type'] = np.select(conditions, values)
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    
       
    def CalculateMicrobiomeRatio(self): 
        """
        Calculate the Beneficial Microbiome Ratio & Harmful Microbiome Ratio 
        Calculate the Number of Beneficial & Harmful Microbiome

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.li_microbiome = list(dict.fromkeys(self.df_dysbiosis['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                harmful_abundance = 0
                beneficial_abundance = 0
                
                harmful_number = 0
                beneficial_number = 0
                
                for j in range(len(self.li_microbiome)):
                    condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) 
                    condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) 
                    
                    if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]    
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]                                  
                                        
                        harmful_abundance += abundance  

                        if abundance > 0:
                            harmful_number += 1
                            
                    elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     
                                
                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]       
                                        
                        beneficial_abundance += abundance   

                        if abundance > 0:
                            beneficial_number += 1
                            
                self.df_eval.loc[self.li_new_sample_name[i], 'harmful_abundance[%]'] = harmful_abundance * 100
                self.df_eval.loc[self.li_new_sample_name[i], 'beneficial_abundance[%]'] = beneficial_abundance * 100
                #self.df_eval.loc[self.li_new_sample_name[i], 'other_abundance[%]'] = 100 - 100 * (harmful_abundance + beneficial_abundance)
                
                self.df_eval.loc[self.li_new_sample_name[i], 'num_harmful_species'] = harmful_number
                self.df_eval.loc[self.li_new_sample_name[i], 'num_beneficial_species'] = beneficial_number

                #self.df_eval.loc[self.li_new_sample_name[i], 'num_total_species'] = self.li_observed[i]
                #self.df_eval.loc[self.li_new_sample_name[i], 'num_other_species'] = self.li_observed[i] - harmful_number - beneficial_number
                              
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg  

    def CalculateAverageMicrobiomeRatio(self): 
        """
        Calculate the average Beneficial Microbiome Ratio & Harmful Microbiome Ratio

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 

            harmful_mean_abundance = 0
            beneficial_mean_abundance = 0

            for j in range(len(self.li_microbiome)):
                condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) 
                condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) 

                if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                    condition_micro = (self.df_db.taxa == self.li_microbiome[j])
                    abundance_mean = 0

                    if (len(self.df_db[condition_micro]) > 0):      
                        abundance_mean += self.df_db[condition_micro].mean(axis=1, numeric_only=True).values[0]      
                        li_micro_sub = []
                        if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                            li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]       

                        harmful_mean_abundance += abundance_mean       

                elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                    condition_micro = (self.df_db.taxa == self.li_microbiome[j])
                    abundance_mean = 0

                    if (len(self.df_db[condition_micro]) > 0):      
                        abundance_mean += self.df_db[condition_micro].mean(axis=1, numeric_only=True).values[0]     
                        li_micro_sub = []
                        if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                            li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]     

                        beneficial_mean_abundance += abundance_mean   

            self.df_eval.loc[:,'harmful_mean_abundance[%]'] = harmful_mean_abundance * 100
            self.df_eval.loc[:,'beneficial_mean_abundance[%]'] = beneficial_mean_abundance * 100
            
            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_egpet_eval_output, encoding="utf-8-sig", index_label='serial_number')
                    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg      
    
####################################
# main
####################################
if __name__ == '__main__':
    
    #path_exp = 'input/PDmirror_output_dog_1629.csv'
    #path_exp = 'input/PCmirror_output_cat_1520.csv'
    
    #path_exp = 'input/PD_dog_one_sample.csv'
    path_exp = 'input/PC_cat_one_sample.csv'
    
    egpetanalysis = EgPetAnalysis(path_exp)
    egpetanalysis.ReadDB()
    egpetanalysis.CalculateMRS()    
    egpetanalysis.CalculateDysbiosis()    
    egpetanalysis.CalculateHealthyDistance()
    egpetanalysis.CalculatePercentileRank()
    egpetanalysis.DrawScatterPlot()    
    egpetanalysis.EvaluatePercentileRank()    
    egpetanalysis.CalculateMicrobiomeRatio()
    egpetanalysis.CalculateAverageMicrobiomeRatio()
    
    print('Analysis Complete')
    