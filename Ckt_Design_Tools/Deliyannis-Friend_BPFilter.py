import math
# Define Design Parameters
f0=1e3  # [Hz]
B=200   # [Hz]
H0=10   # Av DC Gain

# Select arbituary values to solve
C=8.2e-9 # Choose Standard value of C

# Solve for R1, R2, R3
Q=f0/B
R3=Q/((math.pi)*f0*C)
R2=R3/((4*(Q*Q))-(2*H0))
R1=R3/(2*H0)

print(f"\nFor a Deliyannis-Friend Active Bandpass Filter filtering at \n\tf0={f0} Hz with Passband B={B} attenuating DC Gain H0={H0}\n",
      f"\nUse The following Values:\n\tInput Resistor R1={R1} ohms, Shunt Resistor R2={R2} ohms, Feedback Resistor R3={R3} ohms\n",
      f"\tChoosen C={C} F, with Resulting Q=f0/B={Q}\n")