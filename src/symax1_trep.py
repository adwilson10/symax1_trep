#!/usr/bin/env python

# Symax1 Simulator v1
import sys
import math
import roslib; roslib.load_manifest('symax1_trep')
import rospy
import tf

from math import pi as mpi
import numpy as np
import random as rand

import trep
import trep.potentials
import trep.visual as visual
import trep.discopt as discopt

from sensor_msgs.msg import Joy
from std_msgs.msg import Float32MultiArray
from visualization_msgs.msg import Marker

# filter parameters:
meas_cov = np.diag((0.5,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,0.75,0.75)) # measurement covariance
proc_cov = np.diag((0.1,0.1,0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15,0.15)) # process covariance

# initilize joystick global vars
axis = [0.0,0.0,0.0]
buttons=[0,0,0]

# Main Simulation Routine
def QuadMain():

    def generate_desired_trajectory(system, t):
        qd = np.zeros((len(t), system.nQ))
        for i,t in enumerate(t):
            qd[i, 1] = 0
            qd[i, 0] = 0
        return qd

    def make_state_cost(dsys, base, x):
        weight = base*np.ones((dsys.nX,))
        weight[system.get_config('quady').index] = 5
        weight[system.get_config('quadz').index] = 10
        weight[system.get_config('quadx').index] = 5
        weight[system.get_config('quadrx').index] = 10
        weight[system.get_config('quadrz').index] = 10
        weight[system.get_config('quadry').index] = 10

        weight[system.get_config('quady').index+dsys.nX/2] = 1
        weight[system.get_config('quadz').index+dsys.nX/2] = 1
        weight[system.get_config('quadx').index+dsys.nX/2] = 1
        weight[system.get_config('quadrx').index+dsys.nX/2] = 10
        weight[system.get_config('quadrz').index+dsys.nX/2] = 10
        weight[system.get_config('quadry').index+dsys.nX/2] = 10
        return np.diag(weight)

    def make_input_cost(dsys, base, x):
        weight = base*np.ones((dsys.nU,))
       # weight[system.get_input('x-force').index] = x
        return np.diag(weight)

    quad_m = 1.0
    dt = 0.01  # timestep set to 10ms

    # Create a new trep system - 6 generalized coordinates for quadrotor: x,y,z,roll,pitch,yaw
    system = trep.System()
    trep.potentials.Gravity(system, name="Gravity")

    quadz = trep.Frame(system.world_frame, trep.TZ, "quadz")
    quady = trep.Frame(quadz,trep.TY,"quady")
    quadx = trep.Frame(quady,trep.TX,"quadx")

    quadrx = trep.Frame(quadx,trep.RX,"quadrx",kinematic=True)
    quadry = trep.Frame(quadrx,trep.RY,"quadry",kinematic=True)
    quadrz = trep.Frame(quadry,trep.RZ,"quadrz",kinematic=True)
    quadrz.set_mass(quad_m) # central point mass

    # set thrust vector with input u1
    trep.forces.BodyWrench(system,quadry,(0,0,'u1',0,0,0),name='wrench1')

    # set quadrotor initial altitude at 3m
    system.get_config("quadz").q = 0.0

    # Now we'll extract the current configuration into a tuple to use as
    # initial conditions for a variational integrator.
    q0 = system.q

    # Create and initialize the variational integrator
    mvi = trep.MidpointVI(system)
    t = np.arange(0.0, 10.0, 0.01)
    dsys = discopt.DSystem(mvi, t)

    # Generate cost function
    qd = generate_desired_trajectory(system, t)

    Qcost = make_state_cost(dsys, 1, 0.00)
    Rcost = make_input_cost(dsys, 0.01, 0.01)

    u0=np.array([quad_m*9.8])
    u=tuple(u0)

    Uint = np.zeros((len(t)-1,system.nu))
    kin = np.zeros((len(t)-1,system.nQk))

    for i,t in enumerate(t[:-1]):
        Uint[i,:] = u0

    Qk = lambda k: Qcost
    (X,U) = dsys.build_trajectory(u=Uint,rho=kin)
    Kstab = dsys.calc_feedback_controller(X, U, Qk)
    Kstab = Kstab[0]

    mvi.initialize_from_configs(0.0, q0, dt, q0)

    ## These are our simulation variables.
    qref = mvi.q2
    pref = mvi.p2
    vref = (mvi.q2[range(system.nQd,system.nQ)]-mvi.q1[range(system.nQd,system.nQ)])/dt
    
    Xref = np.append(qref,[pref,vref])
    t = np.array([0.0,dt])
    dsys = discopt.DSystem(mvi, t)
    X = Xref

    # start ROS node to broadcast transforms to rviz
    rospy.init_node('state_publisher')
    broadcaster = tf.TransformBroadcaster()
    pub = rospy.Publisher('config',Float32MultiArray)
    markerpub = rospy.Publisher('refmark',Marker)

    configs = Float32MultiArray()
    refmark = Marker()

    # subscribe to joystick topic from joy_node
    rospy.Subscriber("joy",Joy,joycall)

    r = rospy.Rate(100) # simulation rate set to 100Hz

    est_cov = meas_cov


    # Simulator Loop
    while not rospy.is_shutdown():

        # reset simulator if trigger pressed
        if buttons[0]==1:
            X = Xref
            print A

        Xcon = [0,axis[1],axis[0],0,0,0,0,0,0,0,0,0]
        U = np.array([u0[0],0,0,0]-np.dot(Kstab,X-Xref-Xcon)) #+ [rand.gauss(0,0.05) for i in xrange(4)]

        # advance simluation one timestep for model prediction
        dsys.set(X,U,0)  
        X = dsys.f()

        # calculate linearization at current timestep
        A = dsys.fdx()
        P = reduce(np.dot,[A,est_cov,A.T])+proc_cov

        Z = X + [rand.gauss(0,0.05) for i in xrange(12)] #simulated meas noise

        #obtain measurement Z and Kalman gain, K
        y = Z - X
        S = P + meas_cov
        K = np.dot(P,np.linalg.inv(S))

        #update state estimate and cov estimate
        X = X + np.dot(K,y)
        est_cov = np.dot(np.identity(dsys.nX)-K,P)
        
        configs.data = tuple(X)

        refmark.header.frame_id = "/world"
        refmark.header.stamp = rospy.Time.now()
        refmark.type = 2
        refmark.pose.position.z = 3
        refmark.scale.x = 0.15;
        refmark.scale.y = 0.15;
        refmark.scale.z = 0.15;
        refmark.color.r = 0.0;
        refmark.color.g = 1.0;
        refmark.color.b = 0.0;
        refmark.color.a = 1.0;
        refmark.pose.position.y = axis[1]
        refmark.pose.position.x = axis[0]+2

        # send transform data to rviz
        broadcaster.sendTransform((X[2]+2,X[1],X[0]+3),tf.transformations.quaternion_from_euler(X[3], X[4], X[5]),rospy.Time.now(),"quad","world")

        # send transform data to rviz
        broadcaster.sendTransform((Z[2]-2,Z[1],Z[0]+3),tf.transformations.quaternion_from_euler(Z[3], Z[4], Z[5]),rospy.Time.now(),"light_quad","world")
        
        pub.publish(configs)
        markerpub.publish(refmark)
        r.sleep()

# Joystick callback - retrieve current joystick and button data
def joycall(joydata):
    global axis
    global buttons

    axis = joydata.axes
    buttons = joydata.buttons
    axis = [-5.0*axis[0],5.0*axis[1],axis[2]]

# Startup script
if __name__ == '__main__':
    try:
        QuadMain()
    except rospy.ROSInterruptException: pass



