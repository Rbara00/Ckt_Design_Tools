## Installation
Preqs include:
- Python 3.12
- Pip

To Install required libraries:
```
  py -3.12 pip install requirements.txt
```

## The following Repo contains a plethora of miscillenaous design tools for Filters
- filter_sweeper.py is a user-friendly GUI allowing to sweep RLC values and view Bode Plots for 1st order Low Pass/High Pass Filters.
- The GUI also enables the user to export the filter sweep to a spreadsheet.

<img width="2736" height="1738" alt="image" src="https://github.com/user-attachments/assets/81da8608-b693-4fa7-823a-6f921603748d" />

- analog.py contains a plethora of gain calculations and runs within a terminal


## Example circuit for Deliyanni-Friend Active Bandpass Filter
- Deliyannis-Friend_BPFilter.py is an RC/Gain calculator for a narrow active bandpass filter implementation with the following topology. The Deliyannis-Friend Band Pass filter is a modified implemenation of a Sallen-Key Filter that allows the Engineer to create a narrow bandpass fiter with infinite gain of n=2 order, thus allowing for a quick alternative to the typical 2 stage cascading Sallen-Key Bandpass Filter Topology.
- Frequency calculations for the Sallen-Key Filter can more or less be derived using filter_sweeper.py, however I may implement a more focused implementation that helps calculate K values for Gain in the future :)

<img width="478" height="309" alt="image" src="https://github.com/user-attachments/assets/a6790ac9-721f-44dc-b74e-78ab161a6a98" />



