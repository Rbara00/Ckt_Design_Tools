import math
# Define Design Parameters
f0=150  # [Hz]
B=200   # [Hz]
H0=1   # Av DC Gain

# Select arbituary values to solve
C=10e-9 # Choose Standard value of C

# Solve for R1, R2, R3
Q=f0/B
R3=Q/((math.pi)*f0*C)
R2=R3/((4*(Q*Q))-(2*H0))
R1=R3/(2*H0)
print(f"\nFor a Deliyannis-Friend Active Bandpass Filter filtering at \n\tf0={f0} Hz with Passband B={B} attenuating DC Gain H0={H0}\n",
      f"\tChoosen C={C} F, with Resulting Q=f0/B={Q}",
      f"\nUse The following Values:\n\tInput Resistor\t\tR1={R1} ohms\n\tShunt Resistor\t\tR2={R2} ohms\n\tFeedback Resistor\tR3={R3} ohms\n")
# Optional Av adjustment if Non-inverting Negative Feedback is Applied
# Select K from chart based on desired n order filter
#
# Order | K
#---------------------------------
#   2   | 1.586
#   4   | 1.152, 2.235
#   6   | 1.068, 1.586, 2.483
#   8   | 1.038, 1.337, 1.18839, 2.610

# Choose a Reasonable Rf2 (Resistor to GND) and solve for Negative Feedback Resistor Rf1
K=1.586
Rf2=10e3 # [Hz]
Rf1=(K-1)*Rf2
Av=(K*K)*(Rf2/(Rf1+Rf2))

print(f"For Av adjustment, using K={K} add Noninverting Negative Feedback and set:",
      f"\n\tNegative Feedback\tRf1={Rf1} ohms\n\tNegative Feedback Shunt\tRf2={Rf2} ohms \nTherefore,\n\tAv=K^2-(Rf2/(Rf1+Rf2))={Av}\n")

