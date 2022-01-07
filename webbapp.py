import streamlit as st
import plotly.express as px
import numpy as np
import math
#def main():

st.title("Spontaneous Parametric Down-Conversion (SPDC) Simulator ")




url = ' https://drive.google.com/file/d/1lktC_ilHUILRbJFqZH9o6vM3i6RBX6au/view?usp=sharing'
path = 'https://drive.google.com/uc?export=download&id='+url.split('/')[-2]
st.sidebar.image(path, width=310)
st.sidebar.header('Input parameter values')
wavep = st.sidebar.number_input("Wavelength of pump beam (nm)",0.0,2000.0,405.0)
thetap = st.sidebar.number_input("Phase matching angle (degree)",0.0,360.0,28.64)
wp= st.sidebar.number_input("Beam waist (um)",0.0,100000.0,388.0)
L= st.sidebar.number_input("Crystal Thickness (mm)",0.0,100.0, 2.0)
#distz =st.sidebar.number_input("Crystal plane to screen distance(cm)",100) 

distz=100*10000.0
L=L*1000.0
wavep = 0.001*wavep  # wavelength of pump field in um
thetap = np.radians(thetap)  # thetap in radians 28.649



def noo(n):
    global wavep
    n2 = 2.7405 + ((0.0184) / ((n ** 2) - 0.0179)) - (0.0155 * n ** 2)
    rio = round(math.sqrt(n2), 4)
    return rio
 

def neo(n):
    global wavep
    n2 = 2.3730 + ((0.0128) / ((n ** 2) - 0.0156)) - (0.0044 * n ** 2)
    rio = round(math.sqrt(n2), 4)
    return rio
   

def two_photon_wavefunction(xs,ys,xi,yi,wp):
    global L, distz,alphap,betap,gammap,alphas,betas,gammas,alphai,betai,gammai
    ss= (np.square(xs)+ np.square(ys)) + (np.square(xi)+ np.square(yi))+ (2*np.sqrt( (np.square(xs)+ np.square(ys))*(np.square(xi)+ np.square(yi)))*np.cos(np.arctan2(ys,xs)-np.arctan2(yi,xi)))
    kp = (2 * (np.pi)) / wavep
    kpz = (kp * etap) + (alphap * (xs + xi)) - ((1 / (2 * kp * etap)) * (np.square(betap * (xs + xi)) + np.square(gammap * (ys + yi))))
    ksz = (kp * nobar / 2) - ((1 / (kp * nobar)) * (np.square(xs) + np.square(ys)))
    kiz = (kp * nobar / 2) - ((1 / (kp * nobar)) * (np.square(xi) + np.square(yi)))
    tmp = ((ksz + kiz - kpz) * (L / 2))
#-------------------------------------------------------------------------------------
    ksze=((kp/2.0) * etas) + (alphas * (xs)) - ((1 / ( kp * etas)) * (np.square(betas * xs ) + np.square(gammas * ys)))
    tmpeo = ((ksze + kiz - kpz) * (L / 2))
#----------------------------------------------------------------------------------------------
    kize=((kp/2.0) * etai) + (alphai * (xi)) - ((1 / ( kp * etai)) * (np.square(betai * xi) + np.square(gammai * yi)))
    tmpoe = ((ksz + kize - kpz) * (L / 2))
#--------------------------------------------------------------------------------------------
    pumpfield =np.exp(-(wp**2 + ((2j)*(distz/(kp*etap))))*(ss/4.0))
    fullfunction =pumpfield*((np.sinc(tmpeo/(np.pi))*np.exp(-1.0j *tmpeo))+(np.sinc(tmpoe/(np.pi))*np.exp(-1.0j *tmpoe))+np.sinc(tmp/(np.pi))*np.exp(-1.0j *tmp))
    return fullfunction 
#-----------------------------------------------------
#***************************************************
rho=2.0
grdpnt=200
dx= (2*rho/grdpnt)
xs,ys = np.meshgrid(np.linspace(-rho, rho,grdpnt), np.linspace(-rho,rho,grdpnt))
# alpha ,beta ,gamma, eta ***********************
no = noo(wavep)
ne = neo(wavep)
den = np.square(no * np.sin(thetap)) + np.square(ne * np.cos(thetap))
alphap = ((np.square(no) - np.square(ne)) * (np.sin(thetap)) * (np.cos(thetap))) / den
betap = (no * ne) / den
gammap = no / math.sqrt(den)
etap = ne * gammap
nobar = noo(2 * wavep)
#------------------8888888888888------------------------
nos = noo(2*wavep)
nes = neo(2*wavep)
thetas= thetap - np.arctan2(ys,xs)
dens = np.square(nos * np.sin(thetap)) + np.square(nes * np.cos(thetap))
alphas = ((np.square(nos) - np.square(nes)) * (np.sin(thetap)) * (np.cos(thetap))) / dens
betas = (nos * nes) / dens
gammas = nos / np.sqrt(dens)
etas = nes * gammas
nobar = noo(2*wavep)




funcmat=np.zeros([grdpnt,grdpnt])
for k in np.arange(-5,5):
  for m in np.arange(-5,5):

      noi = noo(2*wavep)
      nei = neo(2*wavep)
      thetai= thetap + np.arctan2(-ys + m*dx,-xs + k*dx)
      deni = np.square(noi * np.sin(thetap)) + np.square(nei * np.cos(thetap))
      alphai = ((np.square(noi) - np.square(nei)) * (np.sin(thetap)) * (np.cos(thetap))) / deni
      betai = (noi * nei) / deni
      gammai = noi / np.sqrt(deni)
      etai = nei * gammai
      nobar = noo(2 * wavep)

      funcmat=funcmat+ (np.square(np.absolute(two_photon_wavefunction(xs,ys,-xs+ k*dx,-ys+ m*dx, wp))))*dx*dx





  
#col2.image("spdc_diagram.PNG",width=400,caption='Diagram of SPDC set up')

st.subheader('simulated rings formed by down-converted photons')

fig= px.imshow(funcmat,0.00005,0.00038 )
st.plotly_chart(fig)

st.write('-->Design & Developed by:') 
st.write("Anshul Dhanka")
st.write("Indian Institute of Technology Madras")
st.write("Department of Physics")
