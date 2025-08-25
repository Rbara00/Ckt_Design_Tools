import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Tuple, List, Dict, Any
import warnings
warnings.filterwarnings('ignore')

class AnalogCircuitCalculator:
    """
    A comprehensive analog circuit calculator with user-friendly interface.
    
    Usage:
        calc = AnalogCircuitCalculator()
        calc.amplifier.non_inverting(R1=1000, R2=10000)
        calc.oscillator.wien_bridge(R=1000, C=100e-9)
        calc.plot_summary()
    """
    
    def __init__(self):
        self.amplifier = self.AmplifierCalculator()
        self.oscillator = self.OscillatorCalculator()
        self.power_supply = self.PowerSupplyCalculator()
        self.impedance = self.ImpedanceCalculator()
        self.rf = self.RFCalculator()
        self.control = self.ControlCalculator()
        self.results_history = []
    
    def log_result(self, category: str, function: str, inputs: Dict, outputs: Dict):
        """Log calculation results for history tracking."""
        self.results_history.append({
            'category': category,
            'function': function,
            'inputs': inputs,
            'outputs': outputs
        })
    
    def print_history(self, last_n: int = 5):
        """Print recent calculation history."""
        print(f"\n📊 Recent Calculations (last {min(last_n, len(self.results_history))}):")
        print("=" * 60)
        for result in self.results_history[-last_n:]:
            print(f"{result['category']}.{result['function']}:")
            for key, val in result['inputs'].items():
                if isinstance(val, float) and val > 1000:
                    print(f"  {key}: {val/1000:.1f}k")
                elif isinstance(val, float) and val < 0.001:
                    print(f"  {key}: {val*1e6:.1f}µ")
                else:
                    print(f"  {key}: {val}")
            for key, val in result['outputs'].items():
                if isinstance(val, float):
                    print(f"  → {key}: {val:.3f}")
                else:
                    print(f"  → {key}: {val}")
            print("-" * 40)
    
    class AmplifierCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def non_inverting(self, R1: float, R2: float, verbose: bool = True) -> float:
            """
            Non-inverting op-amp amplifier gain calculation.
            
            Args:
                R1: Resistor to ground (Ω)
                R2: Feedback resistor (Ω)
                verbose: Print formatted result
            
            Returns:
                Voltage gain (V/V)
            
            Example:
                >>> calc.amplifier.non_inverting(R1=1000, R2=10000)
                Non-Inverting Amplifier:
                  R1 (to ground): 1.0 kΩ
                  R2 (feedback): 10.0 kΩ
                  → Gain: 11.0 (20.8 dB)
            """
            gain = 1 + (R2 / R1)
            gain_db = 20 * np.log10(abs(gain))
            
            if verbose:
                print(f"\n🔧 Non-Inverting Amplifier:")
                print(f"  R1 (to ground): {R1/1000:.1f} kΩ")
                print(f"  R2 (feedback): {R2/1000:.1f} kΩ")
                print(f"  → Gain: {gain:.1f} ({gain_db:.1f} dB)")
            
            if self.parent:
                self.parent.log_result('amplifier', 'non_inverting', 
                                     {'R1': R1, 'R2': R2}, 
                                     {'gain': gain, 'gain_dB': gain_db})
            return gain
        
        def inverting(self, R1: float, R2: float, verbose: bool = True) -> float:
            """
            Inverting op-amp amplifier gain calculation.
            
            Args:
                R1: Input resistor (Ω)
                R2: Feedback resistor (Ω)
                verbose: Print formatted result
            
            Returns:
                Voltage gain (V/V) - negative value
            """
            gain = -(R2 / R1)
            gain_db = 20 * np.log10(abs(gain))
            
            if verbose:
                print(f"\n🔧 Inverting Amplifier:")
                print(f"  R1 (input): {R1/1000:.1f} kΩ")
                print(f"  R2 (feedback): {R2/1000:.1f} kΩ")
                print(f"  → Gain: {gain:.1f} ({gain_db:.1f} dB)")
            
            if self.parent:
                self.parent.log_result('amplifier', 'inverting', 
                                     {'R1': R1, 'R2': R2}, 
                                     {'gain': gain, 'gain_dB': gain_db})
            return gain
        
        def differential(self, R1: float, R2: float, R3: float, R4: float, verbose: bool = True) -> Tuple[float, float]:
            """
            Differential amplifier analysis.
            
            Args:
                R1, R3: Input resistors (Ω)
                R2, R4: Feedback resistors (Ω)
            
            Returns:
                (differential_gain, common_mode_gain)
            """
            Ad = R2 / R1
            Ac = (R2/R1) * (R4/(R3+R4)) - (R4/(R3+R4))
            cmrr = 20 * np.log10(abs(Ad/Ac)) if Ac != 0 else float('inf')
            
            if verbose:
                print(f"\n🔧 Differential Amplifier:")
                print(f"  R1: {R1/1000:.1f} kΩ, R2: {R2/1000:.1f} kΩ")
                print(f"  R3: {R3/1000:.1f} kΩ, R4: {R4/1000:.1f} kΩ")
                print(f"  → Differential Gain: {Ad:.1f}")
                print(f"  → Common Mode Gain: {Ac:.4f}")
                print(f"  → CMRR: {cmrr:.1f} dB")
            
            if self.parent:
                self.parent.log_result('amplifier', 'differential', 
                                     {'R1': R1, 'R2': R2, 'R3': R3, 'R4': R4}, 
                                     {'Ad': Ad, 'Ac': Ac, 'CMRR_dB': cmrr})
            return Ad, Ac
    
    class OscillatorCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def wien_bridge(self, R: float, C: float, verbose: bool = True) -> float:
            """
            Wien bridge oscillator frequency calculation.
            
            Args:
                R: Resistance (Ω)
                C: Capacitance (F)
            
            Returns:
                Frequency (Hz)
            """
            frequency = 1 / (2 * np.pi * R * C)
            
            if verbose:
                print(f"\n🌊 Wien Bridge Oscillator:")
                print(f"  R: {R/1000:.1f} kΩ")
                print(f"  C: {C*1e9:.0f} nF")
                print(f"  → Frequency: {frequency:.0f} Hz ({frequency/1000:.2f} kHz)")
            
            if self.parent:
                self.parent.log_result('oscillator', 'wien_bridge', 
                                     {'R': R, 'C': C}, 
                                     {'frequency': frequency})
            return frequency
        
        def timer_555_astable(self, R1: float, R2: float, C: float, verbose: bool = True) -> Tuple[float, float, float]:
            """
            555 timer astable multivibrator calculation.
            
            Args:
                R1: Resistor from Vcc to discharge (Ω)
                R2: Resistor from discharge to threshold (Ω)
                C: Timing capacitor (F)
            
            Returns:
                (frequency_Hz, duty_cycle_percent, period_seconds)
            """
            frequency = 1.44 / ((R1 + 2*R2) * C)
            duty_cycle = (R1 + R2) / (R1 + 2*R2) * 100
            period = 1 / frequency
            
            if verbose:
                print(f"\n⏰ 555 Timer Astable:")
                print(f"  R1: {R1/1000:.1f} kΩ")
                print(f"  R2: {R2/1000:.1f} kΩ")  
                print(f"  C: {C*1e9:.0f} nF")
                print(f"  → Frequency: {frequency:.0f} Hz")
                print(f"  → Duty Cycle: {duty_cycle:.1f}%")
                print(f"  → Period: {period*1000:.2f} ms")
            
            if self.parent:
                self.parent.log_result('oscillator', '555_astable', 
                                     {'R1': R1, 'R2': R2, 'C': C}, 
                                     {'frequency': frequency, 'duty_cycle': duty_cycle, 'period': period})
            return frequency, duty_cycle, period
        
        def timer_555_monostable(self, R: float, C: float, verbose: bool = True) -> float:
            """
            555 timer monostable pulse width calculation.
            
            Args:
                R: Timing resistor (Ω)
                C: Timing capacitor (F)
            
            Returns:
                Pulse width (seconds)
            """
            pulse_width = 1.1 * R * C
            
            if verbose:
                print(f"\n⏰ 555 Timer Monostable:")
                print(f"  R: {R/1000:.1f} kΩ")
                print(f"  C: {C*1e9:.0f} nF")
                print(f"  → Pulse Width: {pulse_width*1000:.2f} ms")
            
            if self.parent:
                self.parent.log_result('oscillator', '555_monostable', 
                                     {'R': R, 'C': C}, 
                                     {'pulse_width': pulse_width})
            return pulse_width
    
    class PowerSupplyCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def voltage_divider(self, Vin: float, R1: float, R2: float, verbose: bool = True) -> float:
            """
            Voltage divider calculation.
            
            Args:
                Vin: Input voltage (V)
                R1: Upper resistor (Ω)
                R2: Lower resistor (Ω)
            
            Returns:
                Output voltage (V)
            """
            Vout = Vin * R2 / (R1 + R2)
            current = Vin / (R1 + R2)
            power = Vin * current
            
            if verbose:
                print(f"\n🔋 Voltage Divider:")
                print(f"  Vin: {Vin:.1f} V")
                print(f"  R1: {R1/1000:.1f} kΩ")
                print(f"  R2: {R2/1000:.1f} kΩ")
                print(f"  → Vout: {Vout:.2f} V")
                print(f"  → Current: {current*1000:.2f} mA")
                print(f"  → Power: {power*1000:.1f} mW")
            
            if self.parent:
                self.parent.log_result('power_supply', 'voltage_divider', 
                                     {'Vin': Vin, 'R1': R1, 'R2': R2}, 
                                     {'Vout': Vout, 'current': current, 'power': power})
            return Vout
        
        def lm317_regulator(self, R1: float, R2: float, Vref: float = 1.25, verbose: bool = True) -> float:
            """
            LM317 adjustable regulator output voltage.
            
            Args:
                R1: Resistor from output to adjust pin (Ω)
                R2: Resistor from adjust pin to ground (Ω)
                Vref: Reference voltage (V)
            
            Returns:
                Output voltage (V)
            """
            Vout = Vref * (1 + R2/R1)
            
            if verbose:
                print(f"\n🔋 LM317 Regulator:")
                print(f"  R1: {R1:.0f} Ω")
                print(f"  R2: {R2/1000:.1f} kΩ")
                print(f"  Vref: {Vref:.2f} V")
                print(f"  → Vout: {Vout:.2f} V")
            
            if self.parent:
                self.parent.log_result('power_supply', 'lm317', 
                                     {'R1': R1, 'R2': R2, 'Vref': Vref}, 
                                     {'Vout': Vout})
            return Vout
        
        def switching_regulator(self, Vin: float, Vout: float, topology: str = 'buck', verbose: bool = True) -> float:
            """
            Switching regulator duty cycle calculation.
            
            Args:
                Vin: Input voltage (V)
                Vout: Output voltage (V)
                topology: 'buck', 'boost', or 'buck-boost'
            
            Returns:
                Duty cycle (0-1)
            """
            if topology.lower() == 'buck':
                duty_cycle = Vout / Vin
            elif topology.lower() == 'boost':
                duty_cycle = 1 - (Vin / Vout)
            elif topology.lower() == 'buck-boost':
                duty_cycle = Vout / (Vin + Vout)
            else:
                raise ValueError("Topology must be 'buck', 'boost', or 'buck-boost'")
            
            if verbose:
                print(f"\n🔋 {topology.title()} Converter:")
                print(f"  Vin: {Vin:.1f} V")
                print(f"  Vout: {Vout:.1f} V")
                print(f"  → Duty Cycle: {duty_cycle:.1%}")
            
            if self.parent:
                self.parent.log_result('power_supply', 'switching_regulator', 
                                     {'Vin': Vin, 'Vout': Vout, 'topology': topology}, 
                                     {'duty_cycle': duty_cycle})
            return duty_cycle
    
    class ImpedanceCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def parallel(self, *Z_values, verbose: bool = True) -> complex:
            """
            Calculate parallel impedance.
            
            Args:
                Z_values: Variable number of impedance values
            
            Returns:
                Total parallel impedance (complex)
            """
            Z_values = [complex(z) for z in Z_values]
            Z_total = 1 / sum(1/z for z in Z_values)
            
            if verbose:
                print(f"\n⚡ Parallel Impedance:")
                for i, z in enumerate(Z_values):
                    print(f"  Z{i+1}: {z:.1f} Ω")
                print(f"  → Total: {Z_total:.1f} Ω")
                print(f"  → Magnitude: {abs(Z_total):.1f} Ω")
                print(f"  → Phase: {np.angle(Z_total, deg=True):.1f}°")
            
            if self.parent:
                self.parent.log_result('impedance', 'parallel', 
                                     {'Z_values': Z_values}, 
                                     {'Z_total': Z_total})
            return Z_total
        
        def series(self, *Z_values, verbose: bool = True) -> complex:
            """
            Calculate series impedance.
            
            Args:
                Z_values: Variable number of impedance values
            
            Returns:
                Total series impedance (complex)
            """
            Z_values = [complex(z) for z in Z_values]
            Z_total = sum(Z_values)
            
            if verbose:
                print(f"\n⚡ Series Impedance:")
                for i, z in enumerate(Z_values):
                    print(f"  Z{i+1}: {z:.1f} Ω")
                print(f"  → Total: {Z_total:.1f} Ω")
                print(f"  → Magnitude: {abs(Z_total):.1f} Ω")
                print(f"  → Phase: {np.angle(Z_total, deg=True):.1f}°")
            
            if self.parent:
                self.parent.log_result('impedance', 'series', 
                                     {'Z_values': Z_values}, 
                                     {'Z_total': Z_total})
            return Z_total
        
        def l_network_matching(self, Rs: float, Rl: float, f: float, verbose: bool = True) -> Tuple[float, float]:
            """
            L-network impedance matching calculation.
            
            Args:
                Rs: Source resistance (Ω)
                Rl: Load resistance (Ω)
                f: Operating frequency (Hz)
            
            Returns:
                (L_henries, C_farads)
            """
            if Rs > Rl:
                Rs, Rl = Rl, Rs
            
            Q = np.sqrt(Rl/Rs - 1)
            omega = 2 * np.pi * f
            
            L = Q * Rs / omega
            C = Q / (omega * Rl)
            
            if verbose:
                print(f"\n⚡ L-Network Matching:")
                print(f"  Rs: {Rs:.0f} Ω")
                print(f"  Rl: {Rl:.0f} Ω")
                print(f"  Frequency: {f/1e6:.1f} MHz")
                print(f"  → L: {L*1e6:.1f} µH")
                print(f"  → C: {C*1e12:.1f} pF")
                print(f"  → Q: {Q:.1f}")
            
            if self.parent:
                self.parent.log_result('impedance', 'l_network', 
                                     {'Rs': Rs, 'Rl': Rl, 'f': f}, 
                                     {'L': L, 'C': C, 'Q': Q})
            return L, C
    
    class RFCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def coax_impedance(self, er: float, D: float, d: float, verbose: bool = True) -> float:
            """
            Coaxial cable characteristic impedance.
            
            Args:
                er: Relative permittivity of dielectric
                D: Inner diameter of outer conductor (mm)
                d: Outer diameter of inner conductor (mm)
            
            Returns:
                Characteristic impedance (Ω)
            """
            Z0 = (138 / np.sqrt(er)) * np.log10(D / d)
            
            if verbose:
                print(f"\n📡 Coaxial Cable:")
                print(f"  εr: {er:.1f}")
                print(f"  D (outer): {D:.1f} mm")
                print(f"  d (inner): {d:.1f} mm")
                print(f"  → Z0: {Z0:.0f} Ω")
            
            if self.parent:
                self.parent.log_result('rf', 'coax_impedance', 
                                     {'er': er, 'D': D, 'd': d}, 
                                     {'Z0': Z0})
            return Z0
        
        def microstrip_impedance(self, er: float, w: float, h: float, verbose: bool = True) -> float:
            """
            Microstrip characteristic impedance.
            
            Args:
                er: Relative permittivity of substrate
                w: Width of trace (mm)
                h: Height of substrate (mm)
            
            Returns:
                Characteristic impedance (Ω)
            """
            ratio = w / h
            
            if ratio >= 1:
                Z0 = (87 / np.sqrt(er + 1.41)) * np.log(5.98 * h / (0.8 * w + h))
            else:
                Z0 = (60 / np.sqrt(er)) * np.log(8 * h / w + 0.25 * w / h)
            
            if verbose:
                print(f"\n📡 Microstrip Trace:")
                print(f"  εr: {er:.1f}")
                print(f"  Width: {w:.2f} mm")
                print(f"  Height: {h:.2f} mm")
                print(f"  w/h ratio: {ratio:.2f}")
                print(f"  → Z0: {Z0:.0f} Ω")
            
            if self.parent:
                self.parent.log_result('rf', 'microstrip_impedance', 
                                     {'er': er, 'w': w, 'h': h}, 
                                     {'Z0': Z0, 'ratio': ratio})
            return Z0
    
    class ControlCalculator:
        def __init__(self):
            self.parent = None
            
        def set_parent(self, parent):
            self.parent = parent
        
        def feedback_gain(self, A_open: float, beta: float, verbose: bool = True) -> float:
            """
            Closed-loop gain with feedback.
            
            Args:
                A_open: Open-loop gain
                beta: Feedback factor (0-1)
            
            Returns:
                Closed-loop gain
            """
            A_closed = A_open / (1 + A_open * beta)
            
            if verbose:
                print(f"\n🔄 Feedback System:")
                print(f"  Open-loop gain: {A_open:.0f}")
                print(f"  Feedback factor: {beta:.3f}")
                print(f"  → Closed-loop gain: {A_closed:.1f}")
                print(f"  → Gain reduction: {A_open/A_closed:.1f}x")
            
            if self.parent:
                self.parent.log_result('control', 'feedback_gain', 
                                     {'A_open': A_open, 'beta': beta}, 
                                     {'A_closed': A_closed})
            return A_closed
    
    def __post_init__(self):
        """Set parent references after initialization."""
        self.amplifier.set_parent(self)
        self.oscillator.set_parent(self)
        self.power_supply.set_parent(self)
        self.impedance.set_parent(self)
        self.rf.set_parent(self)
        self.control.set_parent(self)
    
    def quick_calc(self):
        """Interactive quick calculation mode."""
        print("\n🧮 Quick Calculator Mode")
        print("Available categories:")
        print("1. amplifier    2. oscillator    3. power_supply")
        print("4. impedance    5. rf           6. control")
        print("\nType 'help' for function list, 'quit' to exit")
        
        while True:
            try:
                cmd = input("\n>>> ").strip()
                if cmd.lower() == 'quit':
                    break
                elif cmd.lower() == 'help':
                    self._print_help()
                elif cmd.lower() == 'history':
                    self.print_history()
                else:
                    # Try to execute the command
                    try:
                        result = eval(f"self.{cmd}")
                        if result is not None:
                            print(f"Result: {result}")
                    except Exception as e:
                        print(f"Error: {e}")
                        print("Try: category.function(param1=value1, param2=value2)")
            except KeyboardInterrupt:
                break
    
    def _print_help(self):
        """Print available functions."""
        print("\n📚 Available Functions:")
        print("amplifier.non_inverting(R1=1000, R2=10000)")
        print("amplifier.inverting(R1=1000, R2=10000)")
        print("oscillator.wien_bridge(R=1000, C=100e-9)")
        print("oscillator.timer_555_astable(R1=1000, R2=2200, C=100e-9)")
        print("power_supply.voltage_divider(Vin=12, R1=1000, R2=2000)")
        print("power_supply.lm317_regulator(R1=240, R2=1200)")
        print("impedance.parallel(100, 200, 300)")
        print("rf.coax_impedance(er=2.3, D=5.0, d=1.0)")

# Create a singleton instance for easy access
calc = AnalogCircuitCalculator()
#calc._AnalogCircuitCalculator__post_init__()

# ===================== USAGE EXAMPLES =====================

if __name__ == "__main__":
    print("🔧 Analog Circuit Calculator - Object-Oriented Interface")
    print("=" * 60)
    
    # Example usage
    print("\n📝 Example Calculations:")
    
    # Amplifier calculations
    calc.amplifier.non_inverting(R1=1000, R2=10000)
    calc.amplifier.inverting(R1=2200, R2=22000)
    
    # Oscillator calculations  
    calc.oscillator.wien_bridge(R=1000, C=100e-9)
    calc.oscillator.timer_555_astable(R1=1000, R2=2200, C=100e-9)
    
    # Power supply calculations
    calc.power_supply.voltage_divider(Vin=12, R1=1000, R2=2000)
    calc.power_supply.lm317_regulator(R1=240, R2=1200)
    
    # Impedance calculations
    calc.impedance.parallel(100, 200j, 150-75j)
    calc.impedance.l_network_matching(Rs=50, Rl=200, f=1e6)
    
    # RF calculations
    calc.rf.coax_impedance(er=2.3, D=5.0, d=1.0)
    calc.rf.microstrip_impedance(er=4.2, w=0.5, h=1.6)
    
    # Print calculation history
    calc.print_history()
    
    # Start interactive mode
    print(f"\n🎯 Starting Quick Calculator Mode...")
    print("Try: calc.amplifier.non_inverting(R1=1000, R2=10000)")
    print("Or run: calc.quick_calc() for interactive mode")
    
    # Uncomment the next line to start interactive mode
    # calc.quick_calc()