""" Lab 4 - Exercise 2 """

import numpy as np
import matplotlib.pyplot as plt
from biopack import integrate, DEFAULT, parse_args
from biopack.plot import save_figure
from SystemParameters import MuscleParameters, MassParameters
from lab4_mass import mass_system
import biolog
from scipy.integrate import odeint
# Import muscule model
import Muscle

DEFAULT["label"] = [r"$\theta$ [rad]", r"$d\theta/dt$ [rad/s]"]

# Global settings for plotting
# You may change as per your requirement
plt.rc('lines', linewidth=2.0)
plt.rc('font', size=12.0)
plt.rc('axes', titlesize=14.0)     # fontsize of the axes title
plt.rc('axes', labelsize=14.0)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14.0)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14.0)    # fontsize of the tick labels


def mass_integration(state, time, *args):
    """ Function to integrate muscle-mass system """
    force = args[0]
    mass_parameters = args[2]
    return mass_system(state[0], state[1], force, mass_parameters)


def muscle_integrate(muscle, deltaLength, activation=0.5, dt=0.001):
    """ Function steps or integrates the muscle model by the specified time_step
    dt.

    Parameters:
    -----
        muscle : <Muscle>
            Instance of Muscle class
        deltaLength : float
            Change in Muscle Tendon Length
        activation : float
            Activation of the muscle
        dt : float
            Time step to integrate (Good value is 0.001)

    Returns:
    --------
        res : dict

        res['l_CE'] :
            Contracticle element length
        res['v_CE'] :
            Contracticle element velocity
        res['l_MTC'] :
            Length of muscle tendon unit
        res['activeForce'] :
            Contracticle element force
        res['passiveForce'] :
            Passive element force
        res['force'] :
            activeForce + passiveForce
        res['tendonForce'] :
            Muscle tendon Force

    Example:
    ========
         >>> res = muscle_integrate(muscle, deltaLength=0.0, activation=0.05,
                                    dt=0.01)
    """
    muscle.stim = activation
    muscle.deltaLength = deltaLength
    muscle.step(dt)
    res = {}
    res['l_CE'] = muscle.l_CE
    res['v_CE'] = muscle.v_CE
    res['l_MTC'] = muscle.l_MTC
    res['activeForce'] = muscle.activeForce
    res['passiveForce'] = muscle.passiveForce ## why passive force is so small?? always 0??
    res['force'] = muscle.force
    res['tendonForce'] = muscle.tendonForce
    return res


def isometric_contraction(muscle, stretch, activation):
    """ This function implements the isometric contraction
    of the muscle.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

    Parameters:
    -----------
        muscle : <Muscle>
            Instance of Muscle class
        stretch : list/array
            A list/array of muscle stretches to be evaluated
        activation : float
            Muscle activation

    Returns:
    -------
    """
    # Time settings
    t_start = 0.0  # Start time
    t_stop = 0.2  # Stop time
    dt = 0.001  # Time step
    time = np.arange(t_start,t_stop,dt)
    stretch = np.array(stretch)
    res_p =  np.ones(len(time))
    res_a =  np.ones(len(time))
    maxpassiveforce = np.ones(len(stretch))
    maxactiveforce = np.ones(len(stretch))
    totallength =  np.ones(len(stretch))
    j = 0
    
    for stretch_ in stretch:
        i = 0
        for time_ in time:
            
            res = muscle_integrate(muscle, stretch_, activation, dt) # best = 0.5
            res_p[i]= res['passiveForce']
            res_a[i]= res['activeForce']
            i += 1
        # the value at t = end
#        maxpassiveforce[j]= res_p[len(time)-1]
#        maxactiveforce[j]= res_a[len(time)-1]
        # OR the maximun value during the hole process (smoother??)
        totallength[j] = muscle.l_MTC
        maxpassiveforce[j]= max(res_p[:])
        maxactiveforce[j]= max(res_a[:])
        j += 1
    return totallength, maxactiveforce, maxpassiveforce


def isotonic_contraction(muscle, load,
                         muscle_parameters=MuscleParameters(),
                         mass_parameters=MassParameters()):
    """ This function implements the isotonic contraction
    of the muscle.

    Parameters:
    -----------
        muscle : <Muscle>
            Instance of Muscle class
        load : list/array
            External load to be applied on the muscle.
            It is the mass suspended by the muscle
        muscle_parameters : MuscleParameters
            Muscle paramters instance
        mass_paramters : MassParameters
            Mass parameters instance


    Since the muscle model is complex and sensitive to integration,
    use the following example as a hint to complete your implementation.


    """
    
    # Time settings
    t_start = 0.0  # Start time
    t_stop = 0.2  # Stop time
    dt = 0.001  # Time step
    time = np.arange(t_start,t_stop,dt)
    load = np.array(load)
    activation = 0.2
    j = 0
    res_ = np.ones(len(time))
    velocity_ce = np.ones(len(load))
    x0 = np.array([-0.0, 0.0])
    for load_ in load:
        muscle.initializeMuscleLength()
        mass_parameters.mass = load_ # load is mass (kg)
        state = np.copy(x0)
        for time_ in time:  # before release, iterate to stable state
            res = muscle_integrate(muscle,state[0],activation, dt)
        i = 0
#        print muscle.force
        for time_ in time:  # release, integrate for the state equations
            mass_res = odeint(mass_integration,state,[time_, time_ + dt],args=(muscle.force, load_, mass_parameters)) # returns the position & velocity at each t
            state[0] = mass_res[-1,0]  # position of the end of the muscle
            state[1] = mass_res[-1,1]  # velocity of the end of the muscle
            res = muscle_integrate(muscle,state[0],activation,dt)
            res_[i] = res['v_CE']
            i += 1
        if(res['l_MTC'] > muscle_parameters.l_opt + muscle_parameters.l_slack):
            velocity_ce[j] = min(res_[:])
        else:
            velocity_ce[j] = max(res_[:])
        j+=1
    return load, velocity_ce,activation

def exercise2a():
    """ Exercise 2a
    The goal of this exercise is to understand the relationship
    between muscle length and tension.
    Here you will re-create the isometric muscle contraction experiment.
    To do so, you will have to keep the muscle at a constant length and
    observe the force while stimulating the muscle at a constant activation."""
    # Defination of muscles
    parameters = MuscleParameters()
    biolog.warning("Loading default muscle parameters")
    biolog.info(parameters.showParameters())

    # Create muscle object
    muscle = Muscle.Muscle(parameters)
    muscle.l_opt = 0.11
    activation = 0.05
    stretch=np.arange(-0.09, 0.09, 0.001)
    res = isometric_contraction(muscle, stretch, activation)
    plt.title('Force-Length when activation is %s' %activation)
    plt.title('Force-Length when Lopt is %s' %muscle.l_opt)
    plt.ylim(0,1.7*max(res[1]))
#    plt.ylim(0,2500)
    plt.plot(res[0],res[1])
    plt.plot(res[0],res[2])
    plt.plot(res[0],res[2]+res[1])
    plt.legend(['active force','passive force','total force'])
    plt.xlabel('Length(m)')
    plt.ylabel('Force(N)')
    plt.grid()
    plt.show()
    
    """ Example for plotting graphs using matplotlib. """
    # plt.figure('fig_name')
    # save_figure('fig_name')

def exercise2b():
    """ Exercise 2b
    Under isotonic conditions external load is kept constant.
    A constant stimulation is applied and then suddenly the muscle
    is allowed contract. The instantaneous velocity at which the muscle
    contracts is of our interest"""

    # Defination of muscles
    muscle_parameters = MuscleParameters()
    print(muscle_parameters.showParameters())

    mass_parameters = MassParameters()
    print(mass_parameters.showParameters())

    # Create muscle object
    muscle = Muscle.Muscle(muscle_parameters)
    
    # Run the isotonic_contraction
    load=np.arange(0, 300, 10)
    
    res = isotonic_contraction(muscle,load,muscle_parameters=MuscleParameters(),mass_parameters=MassParameters())
    plt.plot(res[1],res[0])
    plt.title('force-velocity when activation is %s' %res[2])
    plt.ylim(0,60)
    plt.xlabel('Velocity of the contractile element Vce (m/s)')
    plt.ylabel('Force generated by the muscle(N)')
    plt.show()
    
def exercise2():
    """ Exercise 2 """
#    exercise2a()
    exercise2b()


if __name__ == '__main__':
    exercise2()

