# EgPet : Calculate the index and percentile ranks of pets.

![Python](https://img.shields.io/badge/Python-v3.9.0-blue.svg?style=flat&logo=python)&nbsp;
![Pandas](https://img.shields.io/badge/pandas-v2.0.1-blue.svg?style=flat&logo=pandas)&nbsp;
![Numpy](https://img.shields.io/badge/NumPy-v1.24.3-blue.svg?style=flat&logo=numpy)&nbsp;
![Scipy](https://img.shields.io/badge/SciPy-v1.10.1-blue.svg?style=flat&logo=scipy)&nbsp;
![scikit-bio](https://img.shields.io/badge/scikit_bio-grey.svg?style=flat&logo=scikit-bio)&nbsp;

## Installation

You can install the EgPet with following command.
	
	git clone https://github.com/Gyungbu/EgPet.git
 
The list of required packages for `script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded. (version : `python 3.9.0`)
	
	conda create -n env_egpet
	conda activate env_egpet
	conda install pip  
	conda install python=3.9.0
	pip install -r ./EgPet/requirements.txt 

# egpet_update_mrs : (Optional) Update Reference Data
## How to use

### 1. Prepare Merged Proportion File
Place the csv file of your Merged Proportion File in the `./EgPet/input/` folder.

Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.
3. In the case of dog, the name of the Merged Proportion File should start with 'PD', and in the case of cat, the name of the Merged Proportion File should start with 'PC'.

### 2. Run egpet_update_mrs
To run egpet_update_mrs,
 
Run the command below:
  
    python ./EgPet/egpet_update_mrs.py {path_exp}
    ### ex) python egpet_update_mrs.py "/home/kbkim/EgPet/input/PDmirror_output_dog_1629.csv"
    ### ex) python egpet_update_mrs.py "/home/kbkim/EgPet/input/PCmirror_output_cat_1520.csv"
   
    
    
When egpet_update_mrs is executed as above, the file `db_abundance_{self.species}.xlsx`, `egpet_mrs_db_{self.species}.xlsx`, `egpet_percentile_rank_db_{self.species}.csv` will be created or modified in the `./EgPet/input/` folder (where, {self.species} : dog or cat).
And the file `mrs_hist_{self.species}g.png` will be created or modified in the `./EgPet/output/` folder (where, {self.species} : dog or cat).


# egpet_percentile_rank : Calculate the index and percentile ranks of pets.
## How to use

### 1. Prepare Input data
Place the csv file of your proportion file in the `./EgPet/input/` folder.
Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.
3. In the case of dog, the name of the Merged Proportion File should start with 'PD', and in the case of cat, the name of the Merged Proportion File should start with 'PC'.
4. In egpet_percentile_rank.py, under if __name__ == '__main__': enter the path of the proportion file you want to analyze in the 'path_exp' value and save it.

### 2. Run egpet_percentile_rank
To run egpet_percentile_rank,
 
Run the command below:

    python ./EgPet/egpet_percentile_rank.py
    

When egpet_percentile_rank is executed as above, the file `egpet_eval_{self.species}.csv`, `egpet_percentile_rank_{self.species}.csv`, `egpet_harmful_{self.species}.csv`, `egpet_beneficial_{self.species}.csv` and `egpet_scatterplot_dog.png` will be created or modified in the `./EgPet/output/` folder.


