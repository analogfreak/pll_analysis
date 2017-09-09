from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import math
#the function below adds two transfer fucntions
def add_tr(tf1,tf2):
    len_tf1=len(tf1);
    len_tf2=len(tf2);
    sum_tf=np.zeros(max(len_tf1,len_tf2));
    if len_tf1>len_tf2:
        for i in range(0,len_tf1-len_tf2):
            sum_tf[i]=tf1[i];
        for i in range(len_tf1-len_tf2,len_tf1):
            sum_tf[i]=tf1[i]+tf2[i-(len_tf1-len_tf2)]
    elif len_tf2 > len_tf1:
        for i in range(0,len_tf2-len_tf1):
            sum_tf[i]=tf2[i];
        for i in range(len_tf2-len_tf1,len_tf2):
            sum_tf[i]=tf2[i]+tf1[i-(len_tf2-len_tf1)]
    else:
        for i in range(0,len_tf1):
            sum_tf[i]=tf1[i]+tf2[i];
    return sum_tf;


#Defining the parameters here
Icp=3.5172e-6;     #charge pump current
k1=520e6;          #fine path gain (in Hz/V)
k2=8.8e9;          #coarse path gain (in Hz/V)
N=260;             #division factor
c1=4.416e-12;      #Loop filter capacitor (C1)
c2=0.324e-12;      #Loop filter capacitor (C2)
r1=176300;         #Loop filter resistor
#gm_en=0
gm_en=1;
gm=310e-9;         #gain of the transconductance stage
                   #if you don't want to analyse a dual loop PLL
                   #set this parameter to 1e-12
gm_cc=21.9e-12;    #cap at the output of the gm stage
rout_gm=100e9;     #Output impedance of the gm stage


pi=3.141592653589793
f_vec=np.logspace(-1,10,1000); #run the sweep from 0.1Hz to 10GHz
w_vec=2*pi*f_vec;              #convert the corresponding f_vec into its w_vec equivalent
#convert the VCO gains into rads/s-V
k1=k1*2*pi;
k2=k2*2*pi;

#calculaye the various time constants
tau3=gm_cc*rout_gm;
tau1=r1*c1;
tau2=tau1*(c2/(c1+c2));

#define the low gain path transfer fucntion
num_fine=np.array([tau1,1])*Icp*k1; 
den_fine=np.array([tau2,1,0,0])*2*pi*N*(c1+c2)

#define the numerator and denominator of GH(s)
num_gh=np.array([tau1*tau3*k1,(tau1+tau3)*k1,k1+(gm*rout_gm*k2)])*Icp;
den_gh=np.array([tau2*tau3,tau2+tau3,1,0,0])*2*pi*(c1+c2)*N;

#define the numerator and denomiantor of the fine loop
num_gh_fine=np.array([tau1,1])*Icp*k1;
den_gh_fine=np.array([tau2,1,0,0])*2*pi*(c1+c2)*N;
pll_fine_loop=signal.lti(num_gh_fine,den_gh_fine)

#define the numerator and denominator of the coarse loop
num_gh_coarse=gm*rout_gm*k2*Icp;
den_gh_coarse=np.array([tau2*tau3,tau2+tau3,1,0,0])*2*pi*(c1+c2)*N;
pll_coarse_loop=signal.lti(num_gh_coarse,den_gh_coarse);

#Create the PLL open loop transfer function using signal.lti
#signal.lti returns a continious time system object
pll_open_loop=signal.lti(num_gh,den_gh);

#Find the magnitude and phase response of the PLL
w,mag,phase=signal.bode(pll_open_loop,w_vec);
w,mag_fine,phase_fine=signal.bode(pll_fine_loop,w_vec);
w,mag_coarse,phase_coarse=signal.bode(pll_coarse_loop,w_vec)
#Use the intepolate function to make a function w_mag(p) which returns the frequency (in rads/s) corresponding to the point in the
#PLL open loop magnitude curve whose value is 'p' in dB
w_mag=interp1d(mag,w);
w_UGB=w_mag(0); #this finds the frequency(in rads/s) corresponding to UGB

#create two subplots corresponding to magnitude and phase
fig,ax=plt.subplots(nrows=2,ncols=1);
fig.suptitle('PLL Open Loop Response',fontsize=20)
#plot the magnitude response in the first sub-plot
ax_magPlot=plt.subplot(2,1,1);
plt.semilogx(f_vec,mag,linewidth=2,label='Dual path');
plt.legend(['Dual path response'])
plt.semilogx(f_vec,mag_fine,'-',label='Fine path',linewidth=0.5);
plt.legend(['Fine loop response']);
plt.semilogx(f_vec,mag_coarse,'-',label='Coarse path',linewidth=0.5);
plt.legend(['Coarse loop response'])
plt.grid()
plt.legend()
#add a 'x' marker at the UGB point
plt.semilogx(w_UGB/(2*pi),0,'bx',ms=12.5)
ax_magPlot.set_xlabel('Frequency in Hz')
ax_magPlot.set_ylabel('Magnitude in dB')
string_UGB='UGB(MHz): %f'%(w_UGB/(2*pi*1e6))
ax_magPlot.annotate(string_UGB,xy=(w_UGB/(2*pi),0),xytext=(1.5*w_UGB/(2*pi),0))
print string_UGB


ax_phasePlot=plt.subplot(2,1,2,sharex=ax_magPlot);
plt.semilogx(f_vec,phase,label='Dual Path');
ax_phasePlot.set_xlabel('Frequency in Hz')
ax_phasePlot.set_ylabel('Phase in degrees')
w_phase=interp1d(w,phase)
phase_at_UGB=w_phase(w_UGB)
plt.semilogx(w_UGB/(2*pi),phase_at_UGB,'x',ms=12.5)
plt.grid()
plt.legend()

if phase_at_UGB >=0 and phase_at_UGB <180:
    phase_margin=180-phase_at_UGB
elif phase_at_UGB >=180 and phase_at_UGB < 360:
    phase_margin=phase_at_UGB-180;
elif phase_at_UGB >=360 and phase_at_UGB <540:
    phase_margin=540-phase_at_UGB
elif phase_at_UGB < 0:
    phase_margin=phase_at_UGB+180
string_pm="Phase margin = %2.2f"%phase_margin
print string_pm
ax_phasePlot.annotate(string_pm,xy=(w_UGB/(2*pi),phase_at_UGB),xytext=(w_UGB*1.5/(2*pi),phase_at_UGB))

#now the closed loop part.
#forward path transfer function
num_g=num_gh;
den_g=den_gh/N; #den_gh is the denominator of the open loop transfer fucntion, which includes the multiplication by N
#feedback path transfer function
num_h=np.array([1.0])
den_h=np.array([N*1.0])
#the closed loop transfer function
num_pll_cltf=np.convolve(num_g,den_h);
den_pll_cltf=add_tr(np.convolve(den_h,den_g),np.convolve(num_g,num_h));
pll_cltf=signal.lti(num_pll_cltf,den_pll_cltf);
w,mag_cltf,phase_cltf=signal.bode(pll_cltf,w_vec)
plt.figure() #lets create a new plot
plt.semilogx(w/(2*pi),mag_cltf)
dc_value=20*(math.log(N,10));
peaking_in_dB=np.max(mag_cltf)-dc_value;
string_peaking="Peaking= %f dB"% peaking_in_dB
print(string_peaking)
w_peak=np.argmax(mag_cltf);
plt.plot(w[w_peak]/(2*pi),np.max(mag_cltf),'bx',ms=12.5)
#cltf_plot_ax.annotate(string_peaking,xy=(w_peak/(2*pi),np.max(mag_cltf)),xytext=((w_peak/(2*pi))*1.5,np.max(mag_cltf)))
plt.grid()
plt.show()








