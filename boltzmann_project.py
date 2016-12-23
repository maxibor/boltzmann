import subprocess
import numpy as np

###############################################################
# Parsing du fichier de trajectoires pour les composantes des vitesses
# de la forme
# v[ 6950]={-5.55062e-01, -6.19157e-01,  5.54301e-01}
# et ecriture dans le fichier speed.csv
# sous la forme speed_x, speed_y, speed_z
###############################################################

cmd = "gmx dump -f ./data/md_1BTA_100ns_dt2ns.trr"
res = subprocess.check_output(cmd, shell = True)
res = res.split("\n")

speed_x = []
speed_y = []
speed_z = []

i = 0
for line in res[2:] :
    if "v" in line and len(line) > 20 :
        line = line.split("{")
        line = line[1].split(",")
        speed_x.append(float(line[0]))
        speed_y.append(float(line[1]))
        z = float(line[2].split("}")[0])
        speed_z.append(z)
        i += 1


with open("./data/speed.csv","w") as speed_file :
    speed_file.write("x,y,z\n")
    for xs,ys,zs in zip(speed_x,speed_y,speed_z) :
        speed_file.write(str(xs)+","+str(ys)+","+str(zs)+"\n")

###############################################################
# Calcule de la norme de la vitesse
# Et ecriture dans speed_vector.txt
# avec une norme par ligne
###############################################################
speed_vector = []
with open("./data/speed.csv","r") as speed :
    next(speed)
    for line in speed :
        line = line.rstrip()
        line = line.split(",")
        x = float(line[0])**2
        y = float(line[1])**2
        z = float(line[2])**2
        speed_vector.append(np.sqrt(x+y+z))


with open("./data/speed_vector.txt","w") as fileout:
    for elem in speed_vector :
        fileout.write(str(elem)+"\n")




###############################################################
# Parsing des masses des atomes dans le fichier pdb
# et ecriture dans le fichier mass.txt
# avec une masse par ligne pour les 18699 atomes
###############################################################

atom_mass = {"C" : 12, "N" : 14, "O" : 16, "H" : 1, "S" : 32}

mass = []
with open("./data/md_start.pdb", "r") as pdb :
    for line in pdb :
        if "ATOM" in line[0:6] :
            line = line.split("\n")[0]
            line = line[12:16]
            line = line.rstrip()
            if "NA" in line :
                mass.append(23) #Ions Sodium
                continue
            if "NA" not in line  :
                for key in atom_mass.keys() :
                    if key in line :
                        mass.append(atom_mass[key])
                        break

with open("./data/mass.txt", "w") as mf :
    for elem in mass :
        mf.write(str(elem)+"\n")


print "Parsing finished, R analysis starting, please wait..."
###############################################################
# Lancement de l'analyse en R
###############################################################

cmd = "R CMD BATCH boltzmann.r"
subprocess.check_output(cmd, shell = True)
