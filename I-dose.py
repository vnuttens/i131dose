import os
import matplotlib.pyplot as plt
import numpy as np

import lmfit
#--------------------------------------------------------------------------
# Starting values (fill in)

A0 = 100                    # start
AF = 0.5                    # activity reduction factor
bloods = [10,3,2,1]         # blood sample [MBq]
T1 = 42                     # time dialysis 1 [h]
T2 = 90                     # time dialysis 2 [h]
C1 = 0.6                    # MAX clearance after dialysis 1 (60.1%–71.5%) [%]
CL1 = 0.715                 # MIN clearance after dialysis 1 (60.1%–71.5%) [%]
CA1 = (C1+CL1)/2            # clearance after dialysis 1 (60.1%–71.5%) [%]
C2 = 0.3                    # clearance after dialysis 2 (30.4%–67.8%) [%]
CL2 = 0.678
CA2 = (C2+CL2)/2

if not os.path.exists('I-dose resultaten'):
  os.mkdir('I-dose resultaten')
#--------------------------------------------------------------------------
# Constants

# https://jnm.snmjournals.org/content/49/9/1445

T_phys = 192 # half life I-131 [h]
T_eff = 10.5 # eff. half life I-131 for patient with normal kidneys and human thyroid-stimulating hormone [h]
T_eff = 15.7 # eff. half life I-131 for patient with normal kidneys and thyroid hormone withdrawal [h]

#--------------------------------------------------------------------------
# Calculations
a = 0
b = 120
N = 121
x0 = np.linspace(a, b, N)
x1 = np.linspace(a, T2-T1, (T2-T1)+1)
x2 = np.linspace(a, b-T2, (b-T2)+1)
xb = np.linspace(1, len(bloods), len(bloods))

A_phys = A0*(0.5)**(x0/T_phys)
A_normal = A0*(0.5)**(x0/T_eff)

A_dialysis = A0*(1-AF)*(0.5)**(x0/T_phys)
A_dialysisL = A0*(1-AF)*(0.5)**(x0/T_phys)
A_dialysisA = A0*(1-AF)*(0.5)**(x0/T_phys)

A_dialysis[T1:T2+1] = A_dialysis[T1]*(0.5)**(x1/T_phys)*(1-C1)
A_dialysisL[T1:T2+1] = A_dialysisL[T1]*(0.5)**(x1/T_phys)*(1-CL1)
A_dialysisA[T1:T2+1] = A_dialysisA[T1]*(0.5)**(x1/T_phys)*(1-CA1)

A_dialysis[T2:b+1] = A_dialysis[T2]*(0.5)**(x2/T_phys)*(1-C2)
A_dialysisL[T2:b+1] = A_dialysisL[T2]*(0.5)**(x2/T_phys)*(1-CL2)
A_dialysisA[T2:b+1] = A_dialysisA[T2]*(0.5)**(x2/T_phys)*(1-CA2)

AC_check = np.sum(A_phys)*(b-a)/N # must be equal to https://www.wolframalpha.com/input?i=integral%28100*%280.5%29%5E%28x%2F192%29%2C0%2C120%29
AC_dialysis = np.sum(A_dialysis)*(b-a)/N
AC_dialysisL = np.sum(A_dialysisL)*(b-a)/N
AC_dialysisA = np.sum(A_dialysisA)*(b-a)/N
AC_normal = np.sum(A_normal)*(b-a)/N


#--------------------------------------------------------------------------
# Making plot

fig1, ax1  = plt.subplots(1, 1, sharex=True,sharey=True)

ax1.plot(x0,A_normal , color='blue',label='Healthy kidneys, As_normal 100%')

ax1.plot(x0,A_dialysisA, color='red',label=F'Dialysis, As_D/As_N = {"{:0.3f}".format(AC_dialysisA/AC_normal)} with [{"{:0.3f}".format(AC_dialysisL/AC_normal)} , {"{:0.3f}".format(AC_dialysis/AC_normal)}]')
ax1.fill_between(x0, A_dialysisL, A_dialysis,color='red',alpha=0.2)
#ax1.plot(x2,A_D2, color='red')

ax1.plot(x0,A_phys , color='black', linestyle='--',label='Physical decay')

#ax1.axhline(10,color='red')

ax1.set_title('I-131 whole body activity')
ax1.set_ylabel('Relative activity [%]')
ax1.set_xlabel('Tijd [u]')
ax1.legend()

fig2, ax2  = plt.subplots(1, 1, sharex=True,sharey=True)

ax2.plot(xb,bloods , color='blue',label='')
ax2.scatter(xb,bloods , color='blue',label='')

ax2.set_title('I-131 bloed activiteit')
ax2.set_ylabel('Relative activitity [%]')
ax2.set_xlabel('Tijd [u]')
ax2.legend()

fig1.savefig(os.path.join('I-dose resultaten','Whole body dose.pdf'))