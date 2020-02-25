import math

######      Constants     ########

G = 6.67e-11 #Gravitational constant
c_s = 299792458 #speed of light (m/s)
Sm = 2.0e30 # Solar mass (kg)
AU = 1.5e11 # Astronomical Unit (m)
pi = math.pi # PI 3.1415...
Sr = 696.34e6 #Solar radius (m)
Gcm =6.67*1e-8 # Gravitational Constant in CGS

####################################





#####     Central Star Parameters     ########

M = 0.6*Sm #Central Star Mass (Kg)
R = 0.01*Sr #Central Star Radius (m)
B = 0 #Central Star magnetic field (Guass) [If 0 Mag drag is off]
inc = 0 # Field inclination to the spin axis (deg)
TEMP = 5000 #Star Temp (kelvin) [If 0 PR drag is off]
omega = 0.0 #Star spin # find units ;)

##############################################
 






#####     Simulation Parameters     ########

dt = 25.1 #time step (s)
n_steps = 10000  #number of time steps
ACC_RAD = 10*R #Collsion radius
EJE_RAD = 10*AU # Ejection radius

integ = 'RK4' #'basic' will use basic integrator 
                #'yoshi' will use Yoshida method !!!WARNING: UNSTABLE!!!
                #'RK4' will use 4th order runge-kata
                #'RK4_cv' will use 4th order convoluted runge-kata !!!WARNING: UNSTABLE!!!
                #TODO 'BS' will use BS method

PLOT_ON = None  #Change to True to get a real time plot... significantly slows down simulation time!
write_files = None #Change to True to get output files, for plots, or other types of analysis 
OUTPUT_int = 100 # Output file and plot update every nth step


Cores = 1 #n=1 - serial, n>1 - will use n parallel processes, cores='max' use maximum available processors.















