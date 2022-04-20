import os
import shutil as sh

name_dir_vec = ['NoAF']
name_dir = name_dir_vec[0]

Nf = 8

for i in range(Nf):
    i += 1
    #sh.rmtree(name_dir + str(i))
    
    os.mkdir(name_dir + str(i))

    os.system('cp -r Codes' + name_dir + ' ' + name_dir + str(i))

    os.system('mkdir ' + name_dir + str(i) + '/data')
    os.chdir(name_dir + str(i) + '/Codes' + name_dir)
    os.system('ipython cell_geom.py')
    os.chdir('../../')
    
    f = open(name_dir + str(i) + '/Codes' + name_dir + '/compile', 'w')
    f.write('gfortran -Ofast -Wall -O3 geom-python.f90 random.f90 currents.f90 LCC_monte.f90 RyR_monte.f90 -o ' + name_dir + str(i))
    f.close()

    os.chdir(name_dir + str(i) + '/Codes' + name_dir)
    os.system('./compile')
    os.system('nohup ./' + name_dir + str(i) + '&')
    os.chdir('../../')
