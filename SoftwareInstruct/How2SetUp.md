## Basic Info  
In the occasions that you need updated versions of a software that I have used in this repository, or any extra softwares, you will need to install your own.
Follow this instruction to learn how to __set up anaconda or python virtue environment__ , __download wanted software__ and correctly __update your slurm scripts__.
### This instruction is suitable for using a shared HPC system managed by Georgetown University
Our HPC uses a module system and have different versions of anaconda and python installed, and you don't have authority to change change any of these modules. However, you can set up virtue environments using conda or python, and then download your __isolated__ set of sftwares needed for yout project.
* Remember to set up different virtue enviroments if your projects involve using different versions of the same software and dependencies.
* This is to advoid conflicts.

To check the anaconda/python versions available in your shared facility `module avail '

To check what environment you already have `conda info --envs`
* The output will list all the environemnt you already have, including their paths.
* and the * symbol highlights the env you are currently in.
* base:  
  base is the default environment that comes with Anaconda.
  When you open a new terminal or a new session, by default, the base environment is activated unless you specify otherwise.  
  * The path /home/share/apps/python/anaconda3-3.9 shows where this environment is located on your file system. This path contains the Python interpreter and installed packages for the base environment.




