## Basic Info  
In the occasions that you need updated versions of a software that I have used in this repository, or any extra softwares, you will need to install your own.
Follow this instruction to learn how to __set up anaconda or python virtue environment__ , __download wanted software__ and correctly __update your slurm scripts__.
### This instruction is suitable for using a shared HPC system managed by Georgetown University
Our HPC uses a module system and have different versions of anaconda and python installed, and you don't have authority to change change any of these modules. However, you can set up virtue environments using conda or python, and then download your __isolated__ set of sftwares needed for yout project.
* Remember to set up different virtue enviroments if your projects involve using different versions of the same software and dependencies.
* This is to advoid conflicts.

__To check the anaconda/python versions available in your shared facility__     
```module avail```

__To check what environment you already have    
`conda info --envs`__
* The output will list all the environemnt you already have, including their paths.
* and the * symbol highlights the env you are currently in.
* __ base__ : base is the default environment that comes with Anaconda.      
  * When you open a new terminal or a new session, by default, the base environment is activated unless you specify otherwise.  
  * The path /home/share/apps/python/anaconda3-3.9 shows where this environment is located on your file system. This path contains the Python interpreter and installed packages for the base environment.
 
#### Additional Notes:
Using Conda is recommended

__Purpose of Multiple Environments:__    

The main advantage of having multiple conda environments is to maintain different sets of dependencies and packages for different projects or purposes. For instance, one environment can have Python 3.8 with certain libraries for a specific project, while another can have Python 3.9 with a different set of libraries for another project.    
Switching Between Environments:    

To switch to a different environment, you use the command conda activate <env-name>. For example, to switch to the htseq-clip environment, you would use __`conda activate htseq`__    
Managing Environments:    

You can create, remove, or modify environments using various conda commands. This allows you to manage your project dependencies effectively and avoid conflicts between different projects' requirements.    

__How do I check What software version I have in an environment?__
To check the version of htseq or any other Python package installed in a specific virtual environment, you need to first activate that environment and then use a package inspection command.        
__```conda activate htseq```__     
Use __`pip`__ to show version     
__```pip show htseq```__     

__How to upgrade the htseq version in conda__
If you installed HTSeq via conda     
```conda activate htseq   
conda update htseq
```  

If you installed via __pip__        
```conda activate htseq-clip     
pip install --upgrade HTSeq
```
       
* *It won't work if your conda env is managed by others or HPC shared. Scroll down for __install your own virtue env__

__How do I know which software management tool I used to install htseq?__
```conda activate htseq
conda list
```
check if htseq is in the list.
```pip list```
check if htseq is in this list.
* if both list has it, then you probably has installed using one and updated using the other. Be sure to stick with one in the future to manage the same software.

### Create Your Own Virtue Conda Environment
Be sure to check the python version or other dependency versions __required__ for the software you wanted.         
Get off the environment you are currently in if not the __base__: __`conda deactivate`__         
```conda create -n new_htseq_env python=3.8       
conda activate new_htseq_env      
conda install -c bioconda htseq       
conda deactivate
```
### How to Direct to Your Virtue Env Software in Slurm Scripts
#### In a module system like our HPC
```#!/bin/bash
#SBATCH --job-name=htseq-job
#SBATCH --output=htseq_job_%j.out
#SBATCH --error=htseq_job_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Load the module for Anaconda
module load anaconda3

# Activate your conda environment
conda activate htseq-env

# Your HTSeq command
htseq-count [options] <alignment_file> <gff_file>

# Deactivate conda environment
conda deactivate
```

#### Direct to your software in Slur Scripts, in a none-module-managed system
Find the path your Conda is located:
__```which conda```__

```#!/bin/bash
#SBATCH --job-name=htseq-job
# ... other SLURM directives ...

# Load Conda
source /path/to/conda.sh

# Activate Conda environment
conda activate htseq-env

# Run HTSeq
htseq-count [options] <alignment_file> <gff_file>

# Deactivate environment
conda deactivate
```










