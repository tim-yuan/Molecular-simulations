# Molecular-simulations

##     Generate initial position           

```
python gen-pos.py
```

This commend will generate a 8x8x8 box with 500 LJ particles

```
python gen-pos.py -h
```

This will help you generate the configurations you want


##           Energy Minimization           

In the EM folder, you should find 2 files, "EM.f" and "EM.py"
First you need to compile the "EM.f" to a library that can be read by "EM.py", type
```
f2py -c -m EM EM.f
```
and you should see a file named "EM...so"
If you have created an initial configuration using gen-pos.py, type
```
python EM.py -i ../initial.gro
```
Type
```
python EM.py -h 
```
for more options

##           MD with PYTHON                

Type
```
python MD.py -i config.gro
```
to run, or for more information, type
```
python MD.py -h
```

##           MD with PYTHON FORTRAN        

Similar to Energy minimization, we first compile the FORTRAN code to a PYTHON library

```
f2py -c -m mdfun mdfun.f
```

To run the code, type 
```
python MD.py -i config.gro
```

or 

```
python MD.py -h
```

for more information

##     MD with CPP (work in progress)      #
Type 

```
source compile.sh
./main
```

to run the code

##                   MC                   #
Type 

```
python MC.py
```

to run, you will be asked to make some input


