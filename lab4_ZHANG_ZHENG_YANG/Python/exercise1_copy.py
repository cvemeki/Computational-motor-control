""" Lab 4 """

import numpy as np
import matplotlib.pyplot as plt
from biopack import integrate, DEFAULT, parse_args
import biolog
from SystemParameters import PendulumParameters
from lab4_pendulum import pendulum_system

DEFAULT["label"] = [r"$\theta$ [rad]", r"$d\theta/dt$ [rad/s]"]


def pendulum_integration(state, time, *args, **kwargs):
    """ Function for system integration """
#    biolog.warning(
#        "Pendulum equation with spring and damper must be implemented")  # l_S
    return pendulum_system(state[0], state[1], *args, **kwargs)[:, 0]



def exercise1a():
    """ Exercise 1  """
    biolog.info("Executing Lab 4 : Exercise 1");

    def period_(dt,state):
        zero = [0.0,0.0]
        j = 0
        for i in range(1,len(state)):
            if state[i]*state[i-1] <= 0:
                zero[j] = i
                j+=1
            if j == 2:
                break
        period = (zero[1] - zero[0])*dt
        return period
    
    def different_initial_conditions(position,parameters,time):
        for position_ in position:
            x0 = [position_, 0]
            res = integrate(pendulum_integration, x0, time, args=(parameters, ))
            res.plot_state("State")   
            plt.title('State')
            res.plot_phase("Phase")
            plt.title('Phase')
        plt.show()
        return None
    
    def influence_of_k(k,parameters,time,x0,k_test):
        maxposition = np.ones(len(k))
        maxvelocity = np.ones(len(k))
        period = np.ones(len(k))
        
        ''' the change of amplitude and the period in function of k '''
    
        for i in range(0,len(k)): 
            parameters.k1 = k[i]
            parameters.k2 = k[i]
            res = integrate(pendulum_integration, x0, time, args=(parameters, ))
            
            maxposition[i] = max(res.state[:,0])
            maxvelocity[i] = max(res.state[:,1])
            period[i] = period_(dt,res.state[:,1])
            
        plt.plot(k,maxposition)
        plt.plot(k,maxvelocity)
        plt.plot(k,period)
        plt.grid()
        plt.title('The change of the amplitude and the period in function of k')
        plt.legend(['amplitude of position(rad)','amplitude of velocity(rad/s)','period of movement(s)'])
        plt.show()
        
        ''' plot position and velocity seperately for 2 k'''
        for i in range(0,2): 
            for k_ in k_test:
                parameters.k1 = k_
                parameters.k2 = k_
                res = integrate(pendulum_integration, x0, time, args=(parameters, ))
                plt.plot(time,res.state[:,i])
            if i == 0:
                plt.title('Position-time')
            else: plt.title('Velocity-time')
            
            plt.legend(['k1 = %s'%(k_test[0]) , 'k2 = %s'%(k_test[1])])
#            plt.legend(['k = s'])
            plt.grid()
            plt.show()
        return None
    
    def influence_of_theta0(theta0,parameters,time,x0,theta0_test):
        
        maxposition = np.ones(len(theta0))
        minposition = np.ones(len(theta0))
        maxvelocity = np.ones(len(k))
        period = np.ones(len(k))
        
        ''' the change of amplitude, the range of motion, and the period in function of k '''
        for i in range(0,len(theta0)):
            parameters.s_theta_ref1 = theta0[i]
            parameters.s_theta_ref2 = theta0[i]
            res = integrate(pendulum_integration, x0, time, args=(parameters, ))
            
            maxposition[i] = max(res.state[:,0])
            minposition[i] = min(res.state[:,0])
            maxvelocity[i] = max(res.state[:,1])
            period[i] = period_(dt,res.state[:,1])
            
        plt.plot(theta0,maxposition)
        plt.plot(theta0,minposition)
        plt.plot(theta0,maxvelocity)
        plt.plot(theta0,period)
        plt.title('The change of the amplitude, the range of motion, and the period in function of k')
        plt.legend(['amplitude of right position(rad)','amplitude of left position(rad)','amplitude of velocity(rad/s)','period of movement(s)'])
        plt.grid()
        plt.show()
        
        for i in range(0,2): 
            for theta0_ in theta0_test:
                parameters.s_theta_ref1 = theta0_
                parameters.s_theta_ref2 = theta0_
                res = integrate(pendulum_integration, x0, time, args=(parameters, ))
                plt.plot(time,res.state[:,i])
           
            if i == 0:
                plt.title('Position-time')
            else: plt.title('Velocity-time')
        
            plt.legend(['theta01 = %s'%(theta0_test[0]) , 'theta02 = %s'%(theta0_test[1])]) 
            plt.grid()
            plt.show()
         
        return None
     
    def different_k1_k2(k1,k2,parameters,time,x0):
        parameters.k1 = k1
        parameters.k2 = k2
        res = integrate(pendulum_integration, x0, time, args=(parameters, ))
        res.plot_state("State")
        plt.title('State with k1 = %s, k2 = %s' %(k1,k2))
        plt.show()
        return None
    
    def different_theta01_theta02(theta01,theta02,parameters,time,x0):
        parameters.s_theta_ref1 = theta01
        parameters.s_theta_ref2 = theta02
        res = integrate(pendulum_integration, x0, time, args=(parameters, ))
        res.plot_state("State")
        plt.title('State with theta01 = %s, theta02 = %s' %(theta01,theta02))
        plt.show()
        return None
    
    
    # Simulation Parameters
    parameters = PendulumParameters()
    t_start = 0.0
    t_stop = 10.0
    dt = 0.05
    
    biolog.warning("Using large time step dt={}".format(dt))
    time = np.arange(t_start, t_stop, dt)
    
    # no damping
    parameters.b1 = 0
    parameters.b2 = 0
    
    x0 = [np.pi/3, 0]
    
    
    position = [-np.pi/2,0,np.pi/4]
    different_initial_conditions(position,parameters,time)
        
    k = np.linspace(1,60,40)
    k_test = [20,50]
    influence_of_k(k,parameters,time,x0,k_test)
    
    theta0 = np.linspace(-np.pi/2,np.pi/2,40)
    theta0_test = [np.pi/6,np.pi/3]
    influence_of_theta0(theta0,parameters,time,x0,theta0_test)
        
### different k1 and k2
    k1 = 10
    k2 = 30
    different_k1_k2(k1,k2,parameters,time,x0)
    
### different theta01 and theta02
    theta01 = 1.5
    theta02 = 1
    different_theta01_theta02(theta01,theta02,parameters,time,x0)
    
def exercise1b():
    """ Exercise 1  """
    biolog.info("Executing Lab 4 : Exercise 1");
    parameters = PendulumParameters()    
### With damping
    t_start = 0.0
    t_stop = 10.0
    dt = 0.05
    
    biolog.warning("Using large time step dt={}".format(dt))
    time = np.arange(t_start, t_stop, dt)
    
    parameters.b1 = 0.5
    parameters.b2 = 0.5
    x0 = [np.pi/3, 0]
    biolog.info(parameters.showParameters())
    res = integrate(pendulum_integration, x0, time, args=(parameters, ))

    res.plot_state("State")
    
if __name__ == '__main__':
#    exercise1a()
    exercise1b()

