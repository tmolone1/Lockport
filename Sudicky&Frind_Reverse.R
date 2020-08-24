

f<-function (m,x,b) {
  y=m*x+b
  return(y)
}

f(1,2,3)


seepage_velocity<-function (i,K,ne) {
  v=(i*K)/ne
  return(v)
}
fracture_porosity<- function (b,B) {
  thetaf=(b/B)
  return(thetaf)
}
hydrodynamic_dispersion<- function (aL, v, Dstar) {
  D=(aL*v+Dstar)
  return(D)
}
matrix_diffusion<-function (tau, Dstar) {
  Dprime=(tau*Dstar)
  return(Dprime)
}
SF_G<-function(Rprime,Dprime){
  G=(Rprime/Dprime)^(1/2)
  return(G)
}
SF_s<-function(G,B,b){
  s=G*(B-b)
  return(s)
}
SF_g<-function(v,lambda, D, R){
  g=(v^2)/(lambda*D*R)
  return(g)
}
SF_beta<-function(thetam, Rprime, Dprime, lambda, b, R, s){
  beta=(thetam*(Rprime*Dprime)^(1/2))/(lambda^(1/2)*b*R)*tanh(s*lambda^(1/2))
  return(beta)
}
SF_Eq43a<-function(v,g,beta){
  vA=v*((g/2)*(-1+(1+(4/g)*(1+beta))^(1/2)))^(-1)
  return(vA)
}
SF_Eq44<-function(SWSL,R,lambda, z, vA){
  c0=SWSL/(exp(-R*lambda*(z/3.28084)/vA))
  return(c0)
}

SWSL<-1.6e-4 #benzo(a)anthracene, mg/L
z<-seq(0,100,10)


K<-1.9e-3 # Hydraulic Conductivity (cm/s)
i<-2e-3 # Hydraulic gradient (cm/cm)
ne<-0.2 # effective porosity (dimensionless)
b<-4.48e-4 # fracture aperture (m)
B<-4.95e-2 # fracture spacing (m)

v<-seepage_velocity(i,K,ne)/100  # divided by 100 to convert from cm/s to m/s

R<-1 # Face Retardation Coefficient
Rprime<-1 # Matrix Retardation Coefficient
rb<-1.5 #Bulk Density of Porous Matrix, (g/cm3)
Km<-0 # Porous Matrix Distribution Coefficient, (m3/g)
Kf<-0 #Fracture Distribution Coefficient, m
aL<-0.1 #Longitudinal Dispersivity, m
tau<-0.1 #Matrix Tortuosity, (unitless)
thetam<-0.01 #Matrix Porosity (cm3/cm3soil)
thetaf<-fracture_porosity(b,B)

lambda<-5.9027e-9 # First Order Degradation Constant, (sec-1)
Dstar<-9e-10 # Molecular Diffusion Coefficient for a Solute in Free Solution, D* (m2/s)
D<- hydrodynamic_dispersion(aL,v,Dstar) #Hydrodynamic Dispersion Coefficient, D (m2/s)
Dprime<- matrix_diffusion(tau,Dstar)# Diffusion Coefficient, D'
(m2/s)

G<- SF_G(Rprime,Dprime) # S&F equation paramters
s<- SF_s(G,B,b) # S&F equation paramters
g<- SF_g(v,lambda,D,R) # S&F equation paramters
beta<-SF_beta(thetam, Rprime, Dprime, lambda, b, R, s) # S&F equation paramters
vA<-SF_Eq43a(v,g,beta) # SF Equation 43a
c<-SF_Eq44(SWSL,R,lambda, z, vA) # SF Equation 43a
plot(z,c, ylim= c(1e-4, 1e-2))

#High fraacture porosity
B<-0.04957848/2
b<-0.001519/2
thetaf<-fracture_porosity(b,B)
s<- SF_s(G,B,b) # S&F equation paramters
beta<-SF_beta(thetam, Rprime, Dprime, lambda, b, R, s) # S&F equation paramters
vA<-SF_Eq43a(v,g,beta) # SF Equation 43a
c<-SF_Eq44(SWSL,R,lambda, z, vA) # SF Equation 43a
lines(z,c, pch = 24)


#Low fracture porosity 
B<-0.165385/2
b<-0.0005/2
thetaf<-fracture_porosity(b,B)
s<- SF_s(G,B,b) # S&F equation paramters
beta<-SF_beta(thetam, Rprime, Dprime, lambda, b, R, s) # S&F equation paramters
vA<-SF_Eq43a(v,g,beta) # SF Equation 43a
c<-SF_Eq44(SWSL,R,lambda, z, vA) # SF Equation 43a
lines(z,c, pch=25)

# Low K
b<-4.48e-4 # fracture aperture (m)
B<-4.95e-2 # fracture spacing (m)
thetaf<-fracture_porosity(b,B)
K<-8.94e-4 # cm/s, 1st quartile of GWIA F dataset
v<-seepage_velocity(i,K,ne)/100  # divided by 100 to convert from cm/s to m/s
D<- hydrodynamic_dispersion(aL,v,Dstar) #Hydrodynamic Dispersion Coefficient, D (m2/s)
G<- SF_G(Rprime,Dprime) # S&F equation paramters
s<- SF_s(G,B,b) # S&F equation paramters
g<- SF_g(v,lambda,D,R) # S&F equation paramters
beta<-SF_beta(thetam, Rprime, Dprime, lambda, b, R, s) # S&F equation paramters
vA<-SF_Eq43a(v,g,beta) # SF Equation 43a
c<-SF_Eq44(SWSL,R,lambda, z, vA) # SF Equation 43a
lines(z,c, lty=3)

# High K
K<-7.9e-3 # cm/s, 3rd quartile of GWIA F dataset
v<-seepage_velocity(i,K,ne)/100  # divided by 100 to convert from cm/s to m/s
D<- hydrodynamic_dispersion(aL,v,Dstar) #Hydrodynamic Dispersion Coefficient, D (m2/s)
G<- SF_G(Rprime,Dprime) # S&F equation paramters
s<- SF_s(G,B,b) # S&F equation paramters
g<- SF_g(v,lambda,D,R) # S&F equation paramters
beta<-SF_beta(thetam, Rprime, Dprime, lambda, b, R, s) # S&F equation paramters
vA<-SF_Eq43a(v,g,beta) # SF Equation 43a
c<-SF_Eq44(SWSL,R,lambda, z, vA) # SF Equation 43a
lines(z,c, lty=3)

