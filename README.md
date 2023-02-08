# iskay2

Installation:
1. Install modified corrfunc
2. clone this repo to your local system
3. make folder iskay2/iskay2 discoverable in your pythonpath and your system should also be able to find executables in iskay2/misc

You can typically do this in a startup script with the following lines
```
export PYTHONPATH=$PYTHONPATH:$HOME/code/iskay2
export PATH=$PATH:$HOME/code/iskay2/misc
```
4. Run iskay2_gen_iskay2rc.py in the home directory and iskay2_gen_paramfile.py in a project directory. Check there's a .iskay2rc in home and a params.json in your project directory
