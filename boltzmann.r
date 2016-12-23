Navo = 6.022140857*10**23 #nb avogadro
KB = 1.38064852*10^-23 #constante de Boltzmann
NDEF = 37236 #Degrees of freedom in MD

NB_ATOM = 18699 # number of atoms in MD

mass = read.table("./data/mass.txt") # masse molaire pour chaque atome
mass = mass * 10**-3  # conversion en kg par atome

##################################
# Vitesse sur chaque composante
##################################
speeds = read.csv("./data/speed.csv") 
speeds = speeds*1000 #conversion en m/s
speeds = cbind(speeds, rep(mass[,1],51))


##################################
#Norme des vitesses
##################################
speed_norm = read.table("./data/speed_vector.txt")
speed_norm = cbind(speed_norm[,1]*1000, rep(mass[,1],51))




##################################
# Limites des frames
##################################
counter = 0
frame_range = c()
for (j in 1:51) {
    counter = counter+NB_ATOM
    frame_range = append(frame_range,counter)
}


##################################
# Energie par frame
##################################

Ekinetik = c()
for (i in frame_range){
    start = i-NB_ATOM+1
    stop = i
    print(paste("start:",start," - stop:",stop))
    Ekin = c()
    for (j in start:stop) {
        masse = (speed_norm[j,2]) / Navo
        vitesse = speed_norm[j,1]
        Ekin = append(Ekin,(0.5 * masse * (vitesse)**2)) 
    }
    print(sum(Ekin))
    Ekinetik = append(Ekinetik, sum(Ekin))
}
png(filename = "./plot/kinetik.png", width = 1440, height = 1440, units = "px", pointsize = 12)
plot(seq(0,100000, length.out = 51), Ekinetik, main = "Kinetik energy by frame", ylab = "J", xlab = "time (ps)", type = "l" )
dev.off()
#################################
##### Température
#################################

temperature = c()
for (i in Ekinetik){
    temperature = append(temperature, (i / (NDEF* KB))*(2))
}



#############################################
# Equation distribution Maxwell-Boltzmann
#############################################
pv = function(mass, vit, temp){
  return((((mass)/(2*pi*KB*temp))**(1/2)) * exp(- ((mass * (vit**2)) / (2*KB*temp)) ))
}

#############################################
# Verification que l'union des composantes correspond a chaque composante
#############################################

comp_union = union(speeds[,1], union(speeds[,2],speeds[,3]))


#########################################
# Equation Distribution Gaussienne
#########################################
gauss = function(x, mu, sigma) {
    return((1/(sigma*sqrt(2*pi)))*exp(-((x-mu)**2)/(2*sigma**2)))
}


#########################################
# Distribution théorique pour du Carbone a 300K
#########################################

dist_prob = c()
for (i in -3000:3000){
  dist_prob = append(dist_prob, pv( mass= ((12/1000)/Navo), vit = i, temp = 300))
}

png(filename = "./plot/speed_carbon_theo.png", width = 1440, height = 1440, units = "px", pointsize = 12)
plot(-3000:3000, dist_prob, type = "l", xlab = "speed (m/s)", ylab = "probability", main = "Carbon, 300K")


######################################
# Least Square on Gauss
######################################
LS_gauss = function(par, xobs, yobs){
    #par[1] = moyenne
    #par[2] = ecart-type
    sum((((1/(par[2]*sqrt(2*pi)))*exp(-((xobs-par[1])**2)/(2*par[2]**2)))-yobs)**2)
}


######################################
# Application sur distribution théorique Carbon 300K
######################################
xobs = -3000:3000
yobs = pv(mass= ((12/1000)/Navo), vit = -3000:3000, temp = 300)

res = optim(c(10,10), LS_gauss, method = "BFGS", xobs = xobs, yobs = yobs)
points(xobs, gauss(xobs,mu = res$par[1], sigma = res$par[2]), col = "red")
dev.off()

######################################
# Selon maxwell-Boltmzan, ecart type de la gaussiene :
# sigma = sqtr(KB*t/m)
######################################


sigma_temp = function(sigma,mass){
    ((sigma**2)/KB)*(mass/Navo)
}

############################################################################
# Section du point de départ de l'algo de minimisation de LS (BFGS) sur la 1ere frame
############################################################################

xobs = c()
yobs = c()
index = c()
count = 0
for (i in 1:NB_ATOM){
    if (speeds[i,4] == 0.012){
        xobs = append(xobs, speeds[i,1])
        yobs = append(yobs, pv(mass = speeds[i,4]/Navo, vit = speeds[i,1], temp = temperature[1]))
        index = append(index, i)
        count = count + 1
    }
}

best_res = c()
for (i in 200:600){
    res = optim(c(10,i), LS_gauss, method = "BFGS", xobs = xobs, yobs = yobs)
    print(res$value)
    best_res = append(best_res, res$value)
}
best_start = 200+which.min(best_res)
print(best_start)

###############################################
# Application sur les 51 frames pour le carbone
###############################################
counter = 1
all_sd = c()
for (j in frame_range){
    start = j-NB_ATOM+1
    stop = j
    print(paste("start:",start," - stop:",stop))
    xobs = c()
    yobs = c()
    index = c()
    for (i in start:stop){
        if (speeds[i,4] == 0.012){
            xobs = append(xobs, speeds[i,1] )
            yobs = append(yobs, pv(mass = speeds[i,4]/Navo, vit = speeds[i,1], temp = temperature[counter]))
            index = append(index, i)
            count = count + 1
        }
    }

    best_res = c()
    for (i in seq(200,600, by=0.1)){
        res = optim(c(10,i), LS_gauss, method = "BFGS", xobs = xobs, yobs = yobs)
        best_res = append(best_res, res$value)
    }
    best_start = 200+(which.min(best_res)*0.1)
    all_sd = append(all_sd, best_start)
    counter = counter + 1
    print(counter)
}
all_sd



########################
# Temperature Plot
########################
temp_from_sigma = sigma_temp(all_sd,0.012)
png(filename = "./plot/temp_time.png", width = 1440, height = 1440, units = "px", pointsize = 12)
plot(seq(0,100000, length.out = 51), temperature, main = "Temperature variation", ylab = "K", xlab = "time (ps)", type = "l", col = "blue", ylim = c(294,302))
lines(seq(0,100000, length.out = 51), sigma_temp(all_sd, 0.012), col = "orange")
legend("bottomright",legend = c("equipartition","MB-SD"), pch = "-", col = c("blue", "orange"))
dev.off()
##############################
# Speed plot - Gauss and MB
##############################



xobs = c()
yobs = c()
index = c()
count = 0
for (i in 1:NB_ATOM){
    if (speeds[i,4] == 0.012){
        xobs = append(xobs, speeds[i,1] )
        yobs = append(yobs, pv(mass = speeds[i,4]/Navo, vit = speeds[i,1], temp = temperature[1]))
        index = append(index, i)
        count = count + 1
    }
}

best_res = c()
for (i in 200:600){
    res = optim(c(10,i), LS_gauss, method = "BFGS", xobs = xobs, yobs = yobs)
    print(res$value)
    best_res = append(best_res, res$value)
}
best_start = 200+which.min(best_res)


png(filename = "./plot/prob_speed.png", width = 1440, height = 1440, units = "px", pointsize = 12)
hist(xobs, freq = F, prob = T, main="", xlab = "speed", ylab = "probability")
points(xobs,yobs, col = "green", pch = 3)
points(xobs, gauss(xobs, res$par[1], sigma = all_sd), col = "red", pch = 4)
legend("topleft", legend = c("Boltzmann", "Gaussian"), col = c("green", "red"), pch = c(3,4))
dev.off()

